use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WindowId {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AncestryPrior {
    pub groups: Vec<String>,
    pub weights: Vec<f32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantSite {
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub genotype: (u8, u8),
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ReadObservation {
    pub read_id: String,
    pub site_pos: u64,
    pub allele: i8,
    pub baseq: u8,
    pub mapq: u8,
    #[serde(default)]
    pub raw_baseq: Option<u8>,
    #[serde(default)]
    pub weight_v2: Option<f32>,
    #[serde(default)]
    pub weight_v3: Option<f32>,
    #[serde(default)]
    pub local_ref_score: Option<f32>,
    #[serde(default)]
    pub local_alt_score: Option<f32>,
    #[serde(default)]
    pub allele_score_delta: Option<f32>,
    #[serde(default)]
    pub allele_confidence: Option<f32>,
    #[serde(default)]
    pub is_ambiguous: bool,
    #[serde(default)]
    pub context_bases_compared: u16,
    #[serde(default)]
    pub variant_class: String,
    #[serde(default)]
    pub penalty_multiplier: Option<f32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocalPhaseCall {
    pub pos: u64,
    pub phase_state: i8,
    pub confidence: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocalDecodeResult {
    pub window: WindowId,
    pub calls: Vec<LocalPhaseCall>,
    pub ambiguity_flags: Vec<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SitePriorEntry {
    pub pos: u64,
    pub bias: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SiteDonorPriors {
    pub entries: Vec<SitePriorEntry>,
}
