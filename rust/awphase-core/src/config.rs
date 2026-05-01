use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputsConfig {
    pub target_vcf: String,
    pub target_reads: String,
    pub target_variants: String,
    pub target_site_priors: String,
    pub target_block_summaries: Option<String>,
    pub panel_path: String,
    pub truth_vcf: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AncestryConfig {
    pub groups: Vec<String>,
    pub weights: Vec<f32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CandidateEngineConfig {
    pub top_k_donors: usize,
    pub max_composite_states: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DecoderConfig {
    pub window_size: u64,
    pub overlap_size: u64,
    pub lambda_panel: f32,
    pub lambda_anc: f32,
    pub lambda_switch: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BridgeConfig {
    pub min_informative_sites: usize,
    pub min_confidence: f32,
    pub max_gap_bp: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkConfig {
    pub sample_id: String,
    pub chromosome: String,
    pub region: String,
    pub inputs: InputsConfig,
    pub ancestry: AncestryConfig,
    pub candidate_engine: CandidateEngineConfig,
    pub decoder: DecoderConfig,
    pub bridge: BridgeConfig,
}
