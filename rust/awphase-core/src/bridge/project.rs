use serde::{Deserialize, Serialize};

use crate::bridge::stitch::StitchedBlock;
use crate::bridge::summarize::BlockRange;
use crate::types::LocalPhaseCall;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FinalPhaseCall {
    pub pos: u64,
    pub local_phase_state: i8,
    pub stitched_phase_state: i8,
    pub confidence: f32,
    pub block_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FinalPhaseResult {
    pub calls: Vec<FinalPhaseCall>,
}

fn block_for_pos<'a>(pos: u64, ranges: &'a [BlockRange]) -> Option<&'a BlockRange> {
    ranges.iter().find(|r| pos >= r.start && pos <= r.end)
}

fn stitched_orientation_for_block<'a>(block_id: &str, stitched: &'a [StitchedBlock]) -> Option<i8> {
    stitched
        .iter()
        .find(|b| b.block_id == block_id)
        .and_then(|b| b.orientation)
}

pub fn project_stitched_phase_to_calls(
    local_calls: &[LocalPhaseCall],
    block_ranges: &[BlockRange],
    stitched_blocks: &[StitchedBlock],
) -> FinalPhaseResult {
    let calls = local_calls
        .iter()
        .map(|call| {
            let block = block_for_pos(call.pos, block_ranges);

            match block {
                Some(block_range) => {
                    let orientation =
                        stitched_orientation_for_block(&block_range.block_id, stitched_blocks);

                    let stitched_phase_state = match orientation {
                        Some(1) => call.phase_state,
                        Some(-1) => {
                            if call.phase_state == 0 {
                                0
                            } else {
                                -call.phase_state
                            }
                        }
                        None => 0,
                        Some(other) => {
                            if other == 0 {
                                0
                            } else {
                                call.phase_state
                            }
                        }
                    };

                    FinalPhaseCall {
                        pos: call.pos,
                        local_phase_state: call.phase_state,
                        stitched_phase_state,
                        confidence: call.confidence,
                        block_id: block_range.block_id.clone(),
                    }
                }
                None => FinalPhaseCall {
                    pos: call.pos,
                    local_phase_state: call.phase_state,
                    stitched_phase_state: 0,
                    confidence: call.confidence,
                    block_id: "unassigned".to_string(),
                },
            }
        })
        .collect();

    FinalPhaseResult { calls }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bridge::stitch::StitchedBlock;
    use crate::bridge::summarize::BlockRange;
    use crate::types::LocalPhaseCall;

    #[test]
    fn projects_join_and_flip_and_abstain() {
        let local_calls = vec![
            LocalPhaseCall {
                pos: 101,
                phase_state: 1,
                confidence: 0.8,
            },
            LocalPhaseCall {
                pos: 201,
                phase_state: 1,
                confidence: 0.7,
            },
            LocalPhaseCall {
                pos: 301,
                phase_state: -1,
                confidence: 0.6,
            },
            LocalPhaseCall {
                pos: 401,
                phase_state: 1,
                confidence: 0.5,
            },
        ];

        let ranges = vec![
            BlockRange {
                block_id: "block_1".to_string(),
                start: 1,
                end: 150,
            },
            BlockRange {
                block_id: "block_2".to_string(),
                start: 151,
                end: 250,
            },
            BlockRange {
                block_id: "block_3".to_string(),
                start: 251,
                end: 350,
            },
            BlockRange {
                block_id: "block_4".to_string(),
                start: 351,
                end: 450,
            },
        ];

        let stitched = vec![
            StitchedBlock {
                block_id: "block_1".to_string(),
                orientation: Some(1),
                source: "seed".to_string(),
            },
            StitchedBlock {
                block_id: "block_2".to_string(),
                orientation: Some(1),
                source: "join".to_string(),
            },
            StitchedBlock {
                block_id: "block_3".to_string(),
                orientation: Some(-1),
                source: "flip".to_string(),
            },
            StitchedBlock {
                block_id: "block_4".to_string(),
                orientation: None,
                source: "abstain".to_string(),
            },
        ];

        let out = project_stitched_phase_to_calls(&local_calls, &ranges, &stitched);
        assert_eq!(out.calls.len(), 4);

        assert_eq!(out.calls[0].stitched_phase_state, 1);
        assert_eq!(out.calls[1].stitched_phase_state, 1);
        assert_eq!(out.calls[2].stitched_phase_state, 1);
        assert_eq!(out.calls[3].stitched_phase_state, 0);
    }
}
