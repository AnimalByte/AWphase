use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryScoreDebug {
    pub left_block: String,
    pub right_block: String,
    pub same_score: f32,
    pub flip_score: f32,
    pub informative_sites: usize,
    pub decision: String,
    pub confidence: f32,
    pub reason: Option<String>,
}
