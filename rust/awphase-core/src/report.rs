use serde::{Deserialize, Serialize};

use crate::bridge::block_summary::BlockSummary;
use crate::bridge::debug::BoundaryScoreDebug;
use crate::bridge::project::FinalPhaseResult;
use crate::bridge::stitch::StitchedPhaseResult;
use crate::summary::RunSummary;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunReport {
    pub summary: RunSummary,
    pub derived_block_summaries: Vec<BlockSummary>,
    pub derived_boundary_scores: Vec<BoundaryScoreDebug>,
    pub stitched_phase: StitchedPhaseResult,
    pub final_stitched_calls: FinalPhaseResult,
}
