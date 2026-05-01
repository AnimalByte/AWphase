use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::candidate::interface::{CandidateSet, PanelHaplotype};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocalWindowState {
    pub donor_id: String,
    pub hap_id: u8,
    pub group: String,
    pub focal_allele: i8,
    pub local_match_score: f32,
    pub ancestry_weighted_score: f32,
    pub combined_score: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocalPathState {
    pub states: Vec<LocalWindowState>,
    pub path_score: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransitionSummary {
    pub stay_score: f32,
    pub switch_score: f32,
    pub stay_switch_margin: f32,
}

fn normalize_signed(value: f32, denom: f32) -> f32 {
    if denom <= 0.0 {
        0.0
    } else {
        (value / denom).clamp(-1.0, 1.0)
    }
}

pub fn build_local_window_states(
    candidate_set: &CandidateSet,
    panel: &[PanelHaplotype],
    phase_hint: &[i8],
    site_index: usize,
    window_radius: usize,
    top_k_states: usize,
) -> Vec<LocalWindowState> {
    let start = site_index.saturating_sub(window_radius);
    let end = (site_index + window_radius + 1).min(phase_hint.len());

    let mut states: Vec<LocalWindowState> = Vec::new();

    for donor in &candidate_set.donors {
        let maybe_hap = panel
            .iter()
            .find(|h| h.donor_id == donor.donor_id && h.hap_id == donor.hap_id);

        let hap = match maybe_hap {
            Some(h) => h,
            None => continue,
        };

        if site_index >= hap.alleles.len() {
            continue;
        }

        let focal_allele = hap.alleles[site_index];
        if focal_allele == 0 {
            continue;
        }

        let mut local_sum = 0.0_f32;
        let mut local_n = 0.0_f32;

        for j in start..end {
            if j >= hap.alleles.len() {
                break;
            }
            let hint = phase_hint[j];
            let allele = hap.alleles[j];
            if hint == 0 || allele == 0 {
                continue;
            }

            if hint == allele {
                local_sum += 1.0;
            } else {
                local_sum -= 1.0;
            }
            local_n += 1.0;
        }

        let local_match_score = normalize_signed(local_sum, local_n);
        let ancestry_weighted_score = donor.ancestry_weighted_score.max(0.0);
        let combined_score = ancestry_weighted_score * (0.5 + 0.5 * local_match_score.max(0.0));

        states.push(LocalWindowState {
            donor_id: donor.donor_id.clone(),
            hap_id: donor.hap_id,
            group: donor.group.clone(),
            focal_allele,
            local_match_score,
            ancestry_weighted_score,
            combined_score,
        });
    }

    states.sort_by(|a, b| b.combined_score.partial_cmp(&a.combined_score).unwrap());
    states.truncate(top_k_states);
    states
}

fn transition_bonus(prev: &LocalWindowState, curr: &LocalWindowState) -> f32 {
    let same_donor = prev.donor_id == curr.donor_id;
    let same_hap = same_donor && prev.hap_id == curr.hap_id;
    let same_group = prev.group == curr.group;
    let same_allele = prev.focal_allele == curr.focal_allele;

    if same_hap {
        1.0
    } else if same_donor && !same_hap {
        -0.35
    } else if same_group && same_allele {
        0.35
    } else if same_group && !same_allele {
        -0.15
    } else if !same_group && same_allele {
        0.05
    } else {
        -0.30
    }
}

pub fn initialize_path_beam(states: &[LocalWindowState], beam_width: usize) -> Vec<LocalPathState> {
    let mut beam: Vec<LocalPathState> = states
        .iter()
        .map(|s| LocalPathState {
            states: vec![s.clone()],
            path_score: s.combined_score,
        })
        .collect();

    beam.sort_by(|a, b| b.path_score.partial_cmp(&a.path_score).unwrap());
    beam.truncate(beam_width);
    beam
}

pub fn extend_path_beam(
    prev_beam: &[LocalPathState],
    curr_states: &[LocalWindowState],
    lambda_transition: f32,
    beam_width: usize,
) -> Vec<LocalPathState> {
    if prev_beam.is_empty() {
        return initialize_path_beam(curr_states, beam_width);
    }

    let mut candidates = Vec::new();

    for prev_path in prev_beam {
        let prev_state = match prev_path.states.last() {
            Some(s) => s,
            None => continue,
        };

        for curr in curr_states {
            let bonus = transition_bonus(prev_state, curr);
            let transition_term = 1.0 + lambda_transition * bonus;
            let added_score = (curr.combined_score * transition_term).max(0.0);

            let mut new_states = prev_path.states.clone();
            new_states.push(curr.clone());

            candidates.push(LocalPathState {
                states: new_states,
                path_score: prev_path.path_score + added_score,
            });
        }
    }

    let mut best_by_terminal: HashMap<(String, u8), LocalPathState> = HashMap::new();

    for path in candidates {
        let terminal = match path.states.last() {
            Some(s) => (s.donor_id.clone(), s.hap_id),
            None => continue,
        };

        match best_by_terminal.get(&terminal) {
            Some(existing) if existing.path_score >= path.path_score => {}
            _ => {
                best_by_terminal.insert(terminal, path);
            }
        }
    }

    let mut next_beam: Vec<LocalPathState> = best_by_terminal.into_values().collect();
    next_beam.sort_by(|a, b| b.path_score.partial_cmp(&a.path_score).unwrap());
    next_beam.truncate(beam_width);
    next_beam
}

fn normalized_path_weights(beam: &[LocalPathState]) -> Vec<f32> {
    if beam.is_empty() {
        return vec![];
    }

    let max_score = beam
        .iter()
        .map(|p| p.path_score)
        .fold(f32::NEG_INFINITY, f32::max);

    let exps: Vec<f32> = beam
        .iter()
        .map(|p| (p.path_score - max_score).exp())
        .collect();

    let z: f32 = exps.iter().sum();
    if z <= 0.0 {
        vec![0.0; beam.len()]
    } else {
        exps.into_iter().map(|v| v / z).collect()
    }
}

pub fn donor_margin_from_path_beam(beam: &[LocalPathState]) -> f32 {
    let weights = normalized_path_weights(beam);

    let mut pos = 0.0_f32;
    let mut neg = 0.0_f32;

    for (path, w) in beam.iter().zip(weights.iter()) {
        let state = match path.states.last() {
            Some(s) => s,
            None => continue,
        };
        if state.focal_allele > 0 {
            pos += *w;
        } else if state.focal_allele < 0 {
            neg += *w;
        }
    }

    let total = pos + neg;
    if total <= 0.0 {
        0.0
    } else {
        ((pos - neg) / total).clamp(-1.0, 1.0)
    }
}

pub fn summarize_stay_switch(
    prev_beam: &[LocalPathState],
    curr_states: &[LocalWindowState],
    lambda_transition: f32,
) -> TransitionSummary {
    if prev_beam.is_empty() || curr_states.is_empty() {
        return TransitionSummary {
            stay_score: 0.0,
            switch_score: 0.0,
            stay_switch_margin: 0.0,
        };
    }

    let mut best_stay = f32::NEG_INFINITY;
    let mut best_switch = f32::NEG_INFINITY;

    for prev_path in prev_beam {
        let prev_state = match prev_path.states.last() {
            Some(s) => s,
            None => continue,
        };

        for curr in curr_states {
            let bonus = transition_bonus(prev_state, curr);
            let transition_term = 1.0 + lambda_transition * bonus;
            let added_score = (curr.combined_score * transition_term).max(0.0);
            let total_score = prev_path.path_score + added_score;

            let is_stay = prev_state.donor_id == curr.donor_id && prev_state.hap_id == curr.hap_id;
            if is_stay {
                if total_score > best_stay {
                    best_stay = total_score;
                }
            } else if total_score > best_switch {
                best_switch = total_score;
            }
        }
    }

    let stay_score = if best_stay.is_finite() {
        best_stay
    } else {
        0.0
    };
    let switch_score = if best_switch.is_finite() {
        best_switch
    } else {
        0.0
    };

    let denom = stay_score.abs() + switch_score.abs();
    let margin = if denom <= 0.0 {
        0.0
    } else {
        (stay_score - switch_score) / denom
    };

    TransitionSummary {
        stay_score,
        switch_score,
        stay_switch_margin: margin,
    }
}

pub fn current_states_from_beam(beam: &[LocalPathState]) -> Vec<LocalWindowState> {
    let mut out = Vec::new();
    for path in beam {
        if let Some(state) = path.states.last() {
            out.push(state.clone());
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::candidate::interface::{CandidateSet, CandidateState, CompositeState};
    use crate::types::WindowId;

    #[test]
    fn builds_ranked_window_states_and_paths() {
        let candidates = CandidateSet {
            window: WindowId {
                chrom: "chr20".to_string(),
                start: 1,
                end: 1000,
            },
            donors: vec![
                CandidateState {
                    donor_id: "d1".to_string(),
                    hap_id: 0,
                    group: "eur".to_string(),
                    raw_score: 10.0,
                    ancestry_weighted_score: 8.0,
                    matched_sites: 4,
                    total_sites: 4,
                },
                CandidateState {
                    donor_id: "d2".to_string(),
                    hap_id: 0,
                    group: "afr".to_string(),
                    raw_score: 5.0,
                    ancestry_weighted_score: 2.0,
                    matched_sites: 4,
                    total_sites: 4,
                },
            ],
            composites: vec![CompositeState {
                state_id: "top_pair".to_string(),
                members: vec![("d1".to_string(), 0), ("d2".to_string(), 0)],
                score: 5.0,
            }],
        };

        let panel = vec![
            PanelHaplotype {
                donor_id: "d1".to_string(),
                hap_id: 0,
                group: "eur".to_string(),
                alleles: vec![1, 1, 1, 1, 1],
            },
            PanelHaplotype {
                donor_id: "d2".to_string(),
                hap_id: 0,
                group: "afr".to_string(),
                alleles: vec![-1, -1, -1, -1, -1],
            },
        ];

        let phase_hint = vec![1, 1, 1, 1, 1];
        let states = build_local_window_states(&candidates, &panel, &phase_hint, 2, 1, 2);
        assert_eq!(states.len(), 2);

        let beam1 = initialize_path_beam(&states, 2);
        assert_eq!(beam1.len(), 2);

        let beam2 = extend_path_beam(&beam1, &states, 0.15, 2);
        assert!(beam2.len() <= 2);

        let margin = donor_margin_from_path_beam(&beam2);
        assert!(margin > 0.0);

        let ss = summarize_stay_switch(&beam1, &states, 0.15);
        assert!(ss.stay_score >= 0.0);
        assert!(ss.switch_score >= 0.0);
    }
}
