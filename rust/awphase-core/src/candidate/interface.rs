use anyhow::Result;
use serde::{Deserialize, Serialize};

use crate::types::{AncestryPrior, VariantSite, WindowId};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelHaplotype {
    pub donor_id: String,
    pub hap_id: u8,
    pub group: String,
    pub alleles: Vec<i8>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CandidateState {
    pub donor_id: String,
    pub hap_id: u8,
    pub group: String,
    pub raw_score: f32,
    pub ancestry_weighted_score: f32,
    pub matched_sites: usize,
    pub total_sites: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompositeState {
    pub state_id: String,
    pub members: Vec<(String, u8)>,
    pub score: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CandidateSet {
    pub window: WindowId,
    pub donors: Vec<CandidateState>,
    pub composites: Vec<CompositeState>,
}

pub trait CandidateEngine {
    fn retrieve(
        &self,
        window: &WindowId,
        target_variants: &[VariantSite],
        target_phase_guess: &[i8],
        ancestry: &AncestryPrior,
    ) -> Result<CandidateSet>;
}
