use serde::{Deserialize, Serialize};

use crate::decoder::window_state::LocalWindowState;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViterbiStepDebug {
    pub pos: u64,
    pub best_state_index: Option<usize>,
    pub best_score: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViterbiPathResult {
    pub best_path: Vec<Option<LocalWindowState>>,
    pub step_debug: Vec<ViterbiStepDebug>,
}

fn transition_score(prev: &LocalWindowState, curr: &LocalWindowState) -> f32 {
    let same_donor = prev.donor_id == curr.donor_id;
    let same_hap = same_donor && prev.hap_id == curr.hap_id;
    let same_group = prev.group == curr.group;
    let same_allele = prev.focal_allele == curr.focal_allele;

    if same_hap {
        1.5
    } else if same_donor && !same_hap {
        -0.5
    } else if same_group && same_allele {
        0.4
    } else if same_group && !same_allele {
        -0.2
    } else if !same_group && same_allele {
        0.05
    } else {
        -0.35
    }
}

pub fn viterbi_decode_local_states(
    positions: &[u64],
    states_per_site: &[Vec<LocalWindowState>],
    lambda_transition: f32,
) -> ViterbiPathResult {
    if states_per_site.is_empty() {
        return ViterbiPathResult {
            best_path: vec![],
            step_debug: vec![],
        };
    }

    let n_sites = states_per_site.len();
    let mut dp: Vec<Vec<f32>> = vec![vec![]; n_sites];
    let mut backptr: Vec<Vec<Option<usize>>> = vec![vec![]; n_sites];

    if states_per_site[0].is_empty() {
        dp[0] = vec![];
        backptr[0] = vec![];
    } else {
        dp[0] = states_per_site[0]
            .iter()
            .map(|s| s.combined_score)
            .collect();
        backptr[0] = vec![None; states_per_site[0].len()];
    }

    for t in 1..n_sites {
        let curr_states = &states_per_site[t];
        let prev_states = &states_per_site[t - 1];

        if curr_states.is_empty() || prev_states.is_empty() || dp[t - 1].is_empty() {
            dp[t] = vec![];
            backptr[t] = vec![];
            continue;
        }

        let mut curr_dp = vec![f32::NEG_INFINITY; curr_states.len()];
        let mut curr_bp = vec![None; curr_states.len()];

        for (j, curr) in curr_states.iter().enumerate() {
            let emission = curr.combined_score;
            for (i, prev) in prev_states.iter().enumerate() {
                let score =
                    dp[t - 1][i] + emission + lambda_transition * transition_score(prev, curr);
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

    for t in 0..n_sites {
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

    let mut best_path: Vec<Option<LocalWindowState>> = vec![None; n_sites];
    let mut step_debug: Vec<ViterbiStepDebug> = positions
        .iter()
        .map(|&pos| ViterbiStepDebug {
            pos,
            best_state_index: None,
            best_score: 0.0,
        })
        .collect();

    if let (Some(last_t), Some(last_idx)) = (best_last_t, best_last_idx) {
        let mut path_indices: Vec<Option<usize>> = vec![None; n_sites];
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

        for t in 0..n_sites {
            if let Some(idx) = path_indices[t] {
                best_path[t] = states_per_site[t].get(idx).cloned();
                step_debug[t].best_state_index = Some(idx);
                step_debug[t].best_score = dp[t].get(idx).copied().unwrap_or(0.0);
            }
        }
    }

    let _ = best_last_score;

    ViterbiPathResult {
        best_path,
        step_debug,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decoder::window_state::LocalWindowState;

    #[test]
    fn viterbi_prefers_consistent_same_hap_path() {
        let s1a = LocalWindowState {
            donor_id: "d1".to_string(),
            hap_id: 0,
            group: "eur".to_string(),
            focal_allele: 1,
            local_match_score: 1.0,
            ancestry_weighted_score: 5.0,
            combined_score: 5.0,
        };
        let s1b = LocalWindowState {
            donor_id: "d2".to_string(),
            hap_id: 0,
            group: "afr".to_string(),
            focal_allele: -1,
            local_match_score: 0.8,
            ancestry_weighted_score: 4.0,
            combined_score: 4.0,
        };
        let positions = vec![101, 121, 141];
        let states = vec![
            vec![s1a.clone(), s1b.clone()],
            vec![s1a.clone(), s1b.clone()],
            vec![s1a.clone(), s1b.clone()],
        ];

        let out = viterbi_decode_local_states(&positions, &states, 0.3);
        assert_eq!(out.best_path.len(), 3);
        assert_eq!(out.best_path[0].as_ref().unwrap().donor_id, "d1");
        assert_eq!(out.best_path[1].as_ref().unwrap().donor_id, "d1");
        assert_eq!(out.best_path[2].as_ref().unwrap().donor_id, "d1");
    }

    #[test]
    fn handles_empty_state_sites() {
        let s = LocalWindowState {
            donor_id: "d1".to_string(),
            hap_id: 0,
            group: "eur".to_string(),
            focal_allele: 1,
            local_match_score: 1.0,
            ancestry_weighted_score: 5.0,
            combined_score: 5.0,
        };
        let positions = vec![101, 121, 141];
        let states = vec![vec![s.clone()], vec![], vec![s.clone()]];
        let out = viterbi_decode_local_states(&positions, &states, 0.3);
        assert_eq!(out.best_path.len(), 3);
    }
}
