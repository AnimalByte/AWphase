use anyhow::{Context, Result};

use crate::candidate::interface::{
    CandidateEngine, CandidateSet, CandidateState, CompositeState, PanelHaplotype,
};
use crate::io::panel::load_panel_haplotypes;
use crate::types::{AncestryPrior, VariantSite, WindowId};

#[derive(Debug, Clone)]
pub struct PanelDonorScoringEngine {
    pub panel_path: String,
    pub top_k: usize,
    pub max_composite_states: usize,
}

pub type PbwtV1Engine = PanelDonorScoringEngine;

fn ancestry_weight_for_group(group: &str, ancestry: &AncestryPrior) -> f32 {
    ancestry
        .groups
        .iter()
        .zip(ancestry.weights.iter())
        .find(|(g, _)| g.as_str() == group)
        .map(|(_, w)| *w)
        .unwrap_or(0.0)
}

fn hint_match_stats(
    hap: &PanelHaplotype,
    phase_hint: &[i8],
    n_sites: usize,
) -> (f32, usize, usize) {
    let usable = n_sites.min(phase_hint.len()).min(hap.alleles.len());
    let mut raw = 0.0f32;
    let mut matched = 0usize;
    let mut total = 0usize;

    for i in 0..usable {
        let h = hap.alleles[i];
        let q = phase_hint[i];
        if h == 0 || q == 0 {
            continue;
        }
        total += 1;
        if h == q {
            raw += 1.0;
            matched += 1;
        } else {
            raw -= 1.0;
        }
    }

    (raw, matched, total)
}

fn donor_score(
    hap: &PanelHaplotype,
    phase_hint: &[i8],
    n_sites: usize,
    ancestry: &AncestryPrior,
) -> CandidateState {
    let (raw_score, matched_sites, total_sites) = hint_match_stats(hap, phase_hint, n_sites);
    let anc = ancestry_weight_for_group(&hap.group, ancestry);
    let match_frac = if total_sites > 0 {
        matched_sites as f32 / total_sites as f32
    } else {
        0.0
    };

    let ancestry_weighted_score = raw_score + 0.75 * match_frac + 0.50 * anc;

    CandidateState {
        donor_id: hap.donor_id.clone(),
        hap_id: hap.hap_id,
        group: hap.group.clone(),
        raw_score,
        ancestry_weighted_score,
        matched_sites,
        total_sites,
    }
}

fn composite_complementarity(a: &PanelHaplotype, b: &PanelHaplotype, n_sites: usize) -> f32 {
    let usable = n_sites.min(a.alleles.len()).min(b.alleles.len());
    let mut denom = 0usize;
    let mut diff = 0usize;

    for i in 0..usable {
        let x = a.alleles[i];
        let y = b.alleles[i];
        if x == 0 || y == 0 {
            continue;
        }
        denom += 1;
        if x != y {
            diff += 1;
        }
    }

    if denom == 0 {
        0.0
    } else {
        diff as f32 / denom as f32
    }
}

impl CandidateEngine for PanelDonorScoringEngine {
    fn retrieve(
        &self,
        window: &WindowId,
        variants: &[VariantSite],
        phase_hint: &[i8],
        ancestry: &AncestryPrior,
    ) -> Result<CandidateSet> {
        let panel = load_panel_haplotypes(&self.panel_path)
            .with_context(|| format!("failed to load candidate panel from {}", self.panel_path))?;

        let n_sites = variants.len();

        let mut all_donors: Vec<CandidateState> = panel
            .iter()
            .map(|hap| donor_score(hap, phase_hint, n_sites, ancestry))
            .collect();

        all_donors.sort_by(|a, b| {
            b.ancestry_weighted_score
                .partial_cmp(&a.ancestry_weighted_score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    b.raw_score
                        .partial_cmp(&a.raw_score)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| b.matched_sites.cmp(&a.matched_sites))
                .then_with(|| a.donor_id.cmp(&b.donor_id))
                .then_with(|| a.hap_id.cmp(&b.hap_id))
        });

        let donors: Vec<CandidateState> = all_donors.iter().take(self.top_k).cloned().collect();

        let pool_size = (self.top_k * 2).max(2);
        let pool: Vec<CandidateState> = all_donors.iter().take(pool_size).cloned().collect();

        let mut composites_scored: Vec<CompositeState> = Vec::new();

        for i in 0..pool.len() {
            for j in (i + 1)..pool.len() {
                let a = &pool[i];
                let b = &pool[j];

                let hap_a = match panel
                    .iter()
                    .find(|h| h.donor_id == a.donor_id && h.hap_id == a.hap_id)
                {
                    Some(h) => h,
                    None => continue,
                };
                let hap_b = match panel
                    .iter()
                    .find(|h| h.donor_id == b.donor_id && h.hap_id == b.hap_id)
                {
                    Some(h) => h,
                    None => continue,
                };

                let comp = composite_complementarity(hap_a, hap_b, n_sites);
                let score = 0.70 * ((a.ancestry_weighted_score + b.ancestry_weighted_score) / 2.0)
                    + 0.30 * comp;

                composites_scored.push(CompositeState {
                    state_id: format!(
                        "pair_{}_{}_{}_{}",
                        a.donor_id, a.hap_id, b.donor_id, b.hap_id
                    ),
                    members: vec![
                        (a.donor_id.clone(), a.hap_id),
                        (b.donor_id.clone(), b.hap_id),
                    ],
                    score,
                });
            }
        }

        composites_scored.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.state_id.cmp(&b.state_id))
        });
        composites_scored.truncate(self.max_composite_states);

        Ok(CandidateSet {
            window: window.clone(),
            donors,
            composites: composites_scored,
        })
    }
}
