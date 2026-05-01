use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunSummary {
    pub sample_id: String,
    pub region: String,
    pub decoder_mode: String,
    pub number_of_sites: usize,
    pub number_phased: usize,
    pub number_abstained: usize,
    pub informative_boundaries: usize,
    pub joined_boundaries: usize,
    pub abstained_boundaries: usize,
    pub flipped_boundaries: usize,
    pub runtime_seconds: f32,
}
