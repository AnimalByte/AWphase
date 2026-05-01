use serde::{Deserialize, Serialize};

use crate::bridge::block_summary::BlockSummary;
use crate::confidence::interface::Decision;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StitchedBlock {
    pub block_id: String,
    pub orientation: Option<i8>,
    pub source: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StitchedPhaseResult {
    pub blocks: Vec<StitchedBlock>,
}

pub fn stitch_blocks(blocks: &[BlockSummary], decisions: &[Decision]) -> StitchedPhaseResult {
    if blocks.is_empty() {
        return StitchedPhaseResult { blocks: vec![] };
    }

    let mut out = Vec::with_capacity(blocks.len());

    // Seed the first block by preserving its local orientation as-is.
    out.push(StitchedBlock {
        block_id: blocks[0].block_id.clone(),
        orientation: Some(1),
        source: "seed".to_string(),
    });

    for (idx, block) in blocks.iter().enumerate().skip(1) {
        let prev = &out[idx - 1];
        let boundary_decision = decisions.get(idx - 1);

        let stitched = match (prev.orientation, boundary_decision) {
            (Some(prev_ori), Some(Decision::Join)) | (Some(prev_ori), Some(Decision::Keep)) => {
                StitchedBlock {
                    block_id: block.block_id.clone(),
                    orientation: Some(prev_ori),
                    source: "join".to_string(),
                }
            }
            (Some(prev_ori), Some(Decision::Flip)) => StitchedBlock {
                block_id: block.block_id.clone(),
                orientation: Some(-prev_ori),
                source: "flip".to_string(),
            },

            // Abstain breaks propagation, but the next block is still kept in its own local orientation.
            (Some(_), Some(Decision::Abstain)) | (Some(_), None) => StitchedBlock {
                block_id: block.block_id.clone(),
                orientation: Some(1),
                source: "reseed_after_abstain".to_string(),
            },

            // Once propagation has broken, just keep reseeding blocks locally.
            (None, Some(Decision::Join))
            | (None, Some(Decision::Keep))
            | (None, Some(Decision::Flip))
            | (None, Some(Decision::Abstain))
            | (None, None) => StitchedBlock {
                block_id: block.block_id.clone(),
                orientation: Some(1),
                source: "reseed_local".to_string(),
            },
        };

        out.push(stitched);
    }

    StitchedPhaseResult { blocks: out }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bridge::block_summary::BlockSummary;

    fn toy_blocks() -> Vec<BlockSummary> {
        vec![
            BlockSummary {
                block_id: "block_1".to_string(),
                orientation_bias: 1.0,
                support_weight: 2.0,
                informative_sites: 4,
                abstained_sites: 0,
                mean_confidence: 0.6,
                signed_confidence_sum: 2.0,
                coherence: 1.0,
            },
            BlockSummary {
                block_id: "block_2".to_string(),
                orientation_bias: 1.0,
                support_weight: 2.0,
                informative_sites: 4,
                abstained_sites: 0,
                mean_confidence: 0.6,
                signed_confidence_sum: 2.0,
                coherence: 1.0,
            },
            BlockSummary {
                block_id: "block_3".to_string(),
                orientation_bias: -1.0,
                support_weight: 2.0,
                informative_sites: 4,
                abstained_sites: 0,
                mean_confidence: 0.6,
                signed_confidence_sum: -2.0,
                coherence: 1.0,
            },
        ]
    }

    #[test]
    fn join_then_flip_propagates_orientation() {
        let blocks = toy_blocks();
        let stitched = stitch_blocks(&blocks, &[Decision::Join, Decision::Flip]);

        assert_eq!(stitched.blocks.len(), 3);
        assert_eq!(stitched.blocks[0].orientation, Some(1));
        assert_eq!(stitched.blocks[1].orientation, Some(1));
        assert_eq!(stitched.blocks[2].orientation, Some(-1));
    }

    #[test]
    fn abstain_reseeds_next_block_locally() {
        let blocks = toy_blocks();
        let stitched = stitch_blocks(&blocks, &[Decision::Join, Decision::Abstain]);

        assert_eq!(stitched.blocks[0].orientation, Some(1));
        assert_eq!(stitched.blocks[1].orientation, Some(1));
        assert_eq!(stitched.blocks[2].orientation, Some(1));
    }

    #[test]
    fn broken_chain_still_reseeds_downstream_blocks() {
        let blocks = toy_blocks();
        let stitched = stitch_blocks(&blocks, &[Decision::Abstain, Decision::Join]);

        assert_eq!(stitched.blocks[0].orientation, Some(1));
        assert_eq!(stitched.blocks[1].orientation, Some(1));
        assert_eq!(stitched.blocks[2].orientation, Some(1));
    }
}
