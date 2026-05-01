use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Decision {
    Keep,
    Flip,
    Join,
    Abstain,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfidenceScore {
    pub value: f32,
    pub informative_sites: usize,
    pub reason: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryDecision {
    pub left_block: String,
    pub right_block: String,
    pub decision: Decision,
    pub confidence: ConfidenceScore,
}
