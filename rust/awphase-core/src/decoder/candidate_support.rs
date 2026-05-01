// LEGACY: superseded by window_state.rs for current local refinement path.
use crate::candidate::interface::{CandidateSet, PanelHaplotype};

fn normalize_signed(value: f32, denom: f32) -> f32 {
    if denom <= 0.0 {
        0.0
    } else {
        (value / denom).clamp(-1.0, 1.0)
    }
}

pub fn donor_support_for_window(
    candidate_set: &CandidateSet,
    panel: &[PanelHaplotype],
    phase_hint: &[i8],
    site_index: usize,
    window_radius: usize,
) -> f32 {
    let mut weighted_sum = 0.0_f32;
    let mut total_weight = 0.0_f32;

    let start = site_index.saturating_sub(window_radius);
    let end = (site_index + window_radius + 1).min(phase_hint.len());

    for donor in &candidate_set.donors {
        let maybe_hap = panel.iter().find(|h| {
            h.donor_id == donor.donor_id && h.hap_id == donor.hap_id
        });

        let hap = match maybe_hap {
            Some(h) => h,
            None => continue,
        };

        if site_index >= hap.alleles.len() {
            continue;
        }

        let focal_allele = hap.alleles[site_index] as f32;
        if focal_allele == 0.0 {
            continue;
        }

        let mut compat_sum = 0.0_f32;
        let mut compat_n = 0.0_f32;

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
                compat_sum += 1.0;
            } else {
                compat_sum -= 1.0;
            }
            compat_n += 1.0;
        }

        let local_compat = normalize_signed(compat_sum, compat_n);
        let base_weight = donor.ancestry_weighted_score.max(0.0);
        let effective_weight = base_weight * (0.5 + 0.5 * local_compat.max(0.0));

        weighted_sum += effective_weight * focal_allele;
        total_weight += effective_weight;
    }

    normalize_signed(weighted_sum, total_weight)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::candidate::interface::{CandidateSet, CandidateState, CompositeState};
    use crate::types::WindowId;

    #[test]
    fn positive_local_compatibility_yields_positive_support() {
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
        let s = donor_support_for_window(&candidates, &panel, &phase_hint, 2, 1);
        assert!(s > 0.0);
    }

    #[test]
    fn returns_zero_when_no_usable_candidates_exist() {
        let candidates = CandidateSet {
            window: WindowId {
                chrom: "chr20".to_string(),
                start: 1,
                end: 1000,
            },
            donors: vec![],
            composites: vec![],
        };

        let panel = vec![];
        let phase_hint = vec![1, 1, 1];
        let s = donor_support_for_window(&candidates, &panel, &phase_hint, 1, 1);
        assert_eq!(s, 0.0);
    }
}
