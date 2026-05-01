use serde::{Deserialize, Serialize};

use crate::candidate::interface::{CandidateSet, PanelHaplotype};
use crate::decoder::scoring::score_site_from_reads;
use crate::types::ReadObservation;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CandidateWindowState {
    pub donor_id: String,
    pub hap_id: u8,
    pub group: String,
    pub start_pos: u64,
    pub end_pos: u64,
    pub alleles: Vec<i8>,
    pub read_score: f32,
    pub hint_score: f32,
    pub ancestry_score: f32,
    pub emission_score: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindowViterbiResult {
    pub best_path: Vec<Option<CandidateWindowState>>,
    pub end_scores: Vec<Vec<f32>>,
}

pub fn build_candidate_window_states(
    candidate_set: &CandidateSet,
    panel: &[PanelHaplotype],
    phase_hint: &[i8],
    positions: &[u64],
    reads: &[ReadObservation],
    start_idx: usize,
    end_idx: usize,
    top_k_states: usize,
) -> Vec<CandidateWindowState> {
    if start_idx >= end_idx || end_idx > positions.len() {
        return vec![];
    }

    let start_pos = positions[start_idx];
    let end_pos = positions[end_idx - 1];
    let mut states = Vec::new();

    for donor in &candidate_set.donors {
        let hap = match panel
            .iter()
            .find(|h| h.donor_id == donor.donor_id && h.hap_id == donor.hap_id)
        {
            Some(h) => h,
            None => continue,
        };

        if hap.alleles.len() < end_idx {
            continue;
        }

        let mut alleles = Vec::with_capacity(end_idx - start_idx);
        let mut read_score = 0.0_f32;
        let mut hint_score = 0.0_f32;

        for idx in start_idx..end_idx {
            let allele = hap.alleles[idx];
            alleles.push(allele);

            if allele == 0 {
                continue;
            }

            let (alt_support, ref_support, _) = score_site_from_reads(positions[idx], reads);
            let read_balance = alt_support - ref_support;
            read_score += (allele as f32) * read_balance;

            let hint = phase_hint.get(idx).copied().unwrap_or(0);
            if hint != 0 {
                hint_score += (allele as f32) * (hint as f32);
            }
        }

        let ancestry_score = (1.0 + donor.ancestry_weighted_score.max(0.0)).ln();
        let emission_score = read_score + 0.5 * hint_score + 0.25 * ancestry_score;

        states.push(CandidateWindowState {
            donor_id: donor.donor_id.clone(),
            hap_id: donor.hap_id,
            group: donor.group.clone(),
            start_pos,
            end_pos,
            alleles,
            read_score,
            hint_score,
            ancestry_score,
            emission_score,
        });
    }

    states.sort_by(|a, b| b.emission_score.partial_cmp(&a.emission_score).unwrap());
    states.truncate(top_k_states);
    states
}

pub fn transition_score(prev: &CandidateWindowState, curr: &CandidateWindowState) -> f32 {
    let same_donor = prev.donor_id == curr.donor_id;
    let same_hap = same_donor && prev.hap_id == curr.hap_id;
    let same_group = prev.group == curr.group;

    if same_hap {
        1.5
    } else if same_donor && !same_hap {
        -0.75
    } else if same_group {
        0.25
    } else {
        -0.35
    }
}

pub fn viterbi_decode_windows(
    state_windows: &[Vec<CandidateWindowState>],
    lambda_transition: f32,
) -> WindowViterbiResult {
    if state_windows.is_empty() {
        return WindowViterbiResult {
            best_path: vec![],
            end_scores: vec![],
        };
    }

    let n = state_windows.len();
    let mut dp: Vec<Vec<f32>> = vec![vec![]; n];
    let mut backptr: Vec<Vec<Option<usize>>> = vec![vec![]; n];

    for t in 0..n {
        let curr_states = &state_windows[t];
        if curr_states.is_empty() {
            dp[t] = vec![];
            backptr[t] = vec![];
            continue;
        }

        if t == 0 || state_windows[t - 1].is_empty() || dp[t - 1].is_empty() {
            dp[t] = curr_states.iter().map(|s| s.emission_score).collect();
            backptr[t] = vec![None; curr_states.len()];
            continue;
        }

        let prev_states = &state_windows[t - 1];
        let mut curr_dp = vec![f32::NEG_INFINITY; curr_states.len()];
        let mut curr_bp = vec![None; curr_states.len()];

        for (j, curr) in curr_states.iter().enumerate() {
            for (i, prev) in prev_states.iter().enumerate() {
                let score = dp[t - 1][i]
                    + curr.emission_score
                    + lambda_transition * transition_score(prev, curr);
                if score > curr_dp[j] {
                    curr_dp[j] = score;
                    curr_bp[j] = Some(i);
                }
            }
        }

        dp[t] = curr_dp;
        backptr[t] = curr_bp;
    }

    let mut best_last_t: Option<usize> = None;
    let mut best_last_idx: Option<usize> = None;
    let mut best_last_score = f32::NEG_INFINITY;

    for t in 0..n {
        if let Some((idx, score)) = dp[t]
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, s)| (i, *s))
        {
            if score > best_last_score {
                best_last_score = score;
                best_last_t = Some(t);
                best_last_idx = Some(idx);
            }
        }
    }

    let mut best_path: Vec<Option<CandidateWindowState>> = vec![None; n];

    if let (Some(last_t), Some(last_idx)) = (best_last_t, best_last_idx) {
        let mut path_indices: Vec<Option<usize>> = vec![None; n];
        path_indices[last_t] = Some(last_idx);

        let mut t = last_t;
        while t > 0 {
            let idx = match path_indices[t] {
                Some(v) => v,
                None => break,
            };
            path_indices[t - 1] = backptr[t].get(idx).cloned().unwrap_or(None);
            t -= 1;
        }

        for t in 0..n {
            if let Some(idx) = path_indices[t] {
                best_path[t] = state_windows[t].get(idx).cloned();
            }
        }
    }

    WindowViterbiResult {
        best_path,
        end_scores: dp,
    }
}

pub fn normalized_state_weights(scores: &[f32]) -> Vec<f32> {
    if scores.is_empty() {
        return vec![];
    }

    let max_score = scores.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let exps: Vec<f32> = scores.iter().map(|s| (*s - max_score).exp()).collect();
    let z: f32 = exps.iter().sum();

    if z <= 0.0 {
        vec![0.0; scores.len()]
    } else {
        exps.into_iter().map(|v| v / z).collect()
    }
}

pub fn posterior_site_margins(states: &[CandidateWindowState], weights: &[f32]) -> Vec<f32> {
    if states.is_empty() || weights.is_empty() {
        return vec![];
    }

    let n_sites = states[0].alleles.len();
    let mut margins = vec![0.0_f32; n_sites];

    for (state, w) in states.iter().zip(weights.iter()) {
        for (i, allele) in state.alleles.iter().enumerate() {
            margins[i] += *w * (*allele as f32);
        }
    }

    margins.into_iter().map(|m| m.clamp(-1.0, 1.0)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::candidate::interface::{
        CandidateSet, CandidateState, CompositeState, PanelHaplotype,
    };
    use crate::types::{ReadObservation, WindowId};

    #[test]
    fn viterbi_prefers_consistent_same_hap_path() {
        let s1 = CandidateWindowState {
            donor_id: "d1".to_string(),
            hap_id: 0,
            group: "eur".to_string(),
            start_pos: 101,
            end_pos: 141,
            alleles: vec![1, 1],
            read_score: 3.0,
            hint_score: 2.0,
            ancestry_score: 1.0,
            emission_score: 5.0,
        };
        let s2 = CandidateWindowState {
            donor_id: "d2".to_string(),
            hap_id: 0,
            group: "afr".to_string(),
            start_pos: 101,
            end_pos: 141,
            alleles: vec![-1, -1],
            read_score: 2.0,
            hint_score: 1.0,
            ancestry_score: 1.0,
            emission_score: 4.0,
        };
        let states = vec![
            vec![s1.clone(), s2.clone()],
            vec![s1.clone(), s2.clone()],
            vec![s1.clone(), s2.clone()],
        ];
        let out = viterbi_decode_windows(&states, 0.3);
        assert_eq!(out.best_path.len(), 3);
        assert_eq!(out.best_path[0].as_ref().unwrap().donor_id, "d1");
        assert_eq!(out.best_path[1].as_ref().unwrap().donor_id, "d1");
        assert_eq!(out.best_path[2].as_ref().unwrap().donor_id, "d1");
    }

    #[test]
    fn builds_candidate_window_states() {
        let candidates = CandidateSet {
            window: WindowId {
                chrom: "chr20".to_string(),
                start: 1,
                end: 1000,
            },
            donors: vec![CandidateState {
                donor_id: "d1".to_string(),
                hap_id: 0,
                group: "eur".to_string(),
                raw_score: 1.0,
                ancestry_weighted_score: 2.0,
                matched_sites: 2,
                total_sites: 2,
            }],
            composites: vec![CompositeState {
                state_id: "x".to_string(),
                members: vec![("d1".to_string(), 0)],
                score: 1.0,
            }],
        };

        let panel = vec![PanelHaplotype {
            donor_id: "d1".to_string(),
            hap_id: 0,
            group: "eur".to_string(),
            alleles: vec![1, 1, -1, -1],
        }];

        let phase_hint = vec![1, 1, -1, -1];
        let positions = vec![101, 121, 141, 161];
        let reads: Vec<ReadObservation> = vec![];

        let states = build_candidate_window_states(
            &candidates,
            &panel,
            &phase_hint,
            &positions,
            &reads,
            0,
            2,
            3,
        );

        assert_eq!(states.len(), 1);
        assert_eq!(states[0].alleles.len(), 2);
    }
}
