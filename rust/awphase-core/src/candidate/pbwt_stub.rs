use anyhow::{Context, Result};

use crate::candidate::interface::{
    CandidateEngine, CandidateSet, CandidateState, CompositeState, PanelHaplotype,
};
use crate::io::panel::load_panel_haplotypes_json;
use crate::types::{AncestryPrior, VariantSite, WindowId};

pub struct PbwtStubEngine {
    pub panel_json_path: String,
    pub top_k: usize,
}

fn ancestry_weight(group: &str, ancestry: &AncestryPrior) -> f32 {
    ancestry
        .groups
        .iter()
        .zip(ancestry.weights.iter())
        .find(|(g, _)| g.as_str() == group)
        .map(|(_, w)| *w)
        .unwrap_or(0.0)
}

fn score_haplotype(hap: &PanelHaplotype, target_phase_guess: &[i8]) -> (f32, usize, usize) {
    let n = target_phase_guess.len().min(hap.alleles.len());
    let mut matched = 0usize;
    let mut raw = 0.0_f32;

    for i in 0..n {
        let t = target_phase_guess[i];
        let h = hap.alleles[i];
        if t == 0 || h == 0 {
            continue;
        }
        if t == h {
            raw += 1.0;
            matched += 1;
        } else {
            raw -= 1.0;
        }
    }

    (raw, matched, n)
}

impl CandidateEngine for PbwtStubEngine {
    fn retrieve(
        &self,
        window: &WindowId,
        _target_variants: &[VariantSite],
        target_phase_guess: &[i8],
        ancestry: &AncestryPrior,
    ) -> Result<CandidateSet> {
        let panel = load_panel_haplotypes_json(&self.panel_json_path).with_context(|| {
            format!(
                "failed to load candidate panel from {}",
                self.panel_json_path
            )
        })?;

        let mut donors: Vec<CandidateState> = panel
            .iter()
            .map(|hap| {
                let (raw_score, matched_sites, total_sites) =
                    score_haplotype(hap, target_phase_guess);
                let weight = ancestry_weight(&hap.group, ancestry);
                let ancestry_weighted_score = raw_score * (0.5 + weight);

                CandidateState {
                    donor_id: hap.donor_id.clone(),
                    hap_id: hap.hap_id,
                    group: hap.group.clone(),
                    raw_score,
                    ancestry_weighted_score,
                    matched_sites,
                    total_sites,
                }
            })
            .collect();

        donors.sort_by(|a, b| {
            b.ancestry_weighted_score
                .partial_cmp(&a.ancestry_weighted_score)
                .unwrap()
        });
        donors.truncate(self.top_k);

        let composites = if donors.len() >= 2 {
            vec![CompositeState {
                state_id: "top_pair".to_string(),
                members: vec![
                    (donors[0].donor_id.clone(), donors[0].hap_id),
                    (donors[1].donor_id.clone(), donors[1].hap_id),
                ],
                score: (donors[0].ancestry_weighted_score + donors[1].ancestry_weighted_score)
                    / 2.0,
            }]
        } else {
            vec![]
        };

        Ok(CandidateSet {
            window: window.clone(),
            donors,
            composites,
        })
    }
}
