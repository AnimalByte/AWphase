use crate::bridge::summarize::BlockRange;
use crate::types::LocalPhaseCall;

pub fn make_blocks_from_local_calls(calls: &[LocalPhaseCall]) -> Vec<BlockRange> {
    let mut ordered: Vec<&LocalPhaseCall> = calls.iter().collect();
    ordered.sort_by_key(|c| c.pos);

    if ordered.is_empty() {
        return vec![];
    }

    let mut out: Vec<BlockRange> = Vec::new();
    let mut block_idx = 1usize;

    let mut in_block = false;
    let mut current_start = 0u64;
    let mut current_end = 0u64;

    for call in ordered {
        if call.phase_state == 0 {
            if in_block {
                out.push(BlockRange {
                    block_id: format!("block_{}", block_idx),
                    start: current_start,
                    end: current_end,
                });
                block_idx += 1;
                in_block = false;
            }
            continue;
        }

        if !in_block {
            current_start = call.pos;
            current_end = call.pos;
            in_block = true;
        } else {
            current_end = call.pos;
        }
    }

    if in_block {
        out.push(BlockRange {
            block_id: format!("block_{}", block_idx),
            start: current_start,
            end: current_end,
        });
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::LocalPhaseCall;

    #[test]
    fn groups_contiguous_informative_calls_into_blocks_even_if_sign_changes() {
        let calls = vec![
            LocalPhaseCall {
                pos: 101,
                phase_state: 1,
                confidence: 0.8,
            },
            LocalPhaseCall {
                pos: 121,
                phase_state: -1,
                confidence: 0.8,
            },
            LocalPhaseCall {
                pos: 141,
                phase_state: 1,
                confidence: 0.7,
            },
            LocalPhaseCall {
                pos: 201,
                phase_state: 0,
                confidence: 0.5,
            },
            LocalPhaseCall {
                pos: 221,
                phase_state: -1,
                confidence: 0.7,
            },
            LocalPhaseCall {
                pos: 241,
                phase_state: 1,
                confidence: 0.6,
            },
        ];

        let blocks = make_blocks_from_local_calls(&calls);
        assert_eq!(blocks.len(), 2);

        assert_eq!(blocks[0].block_id, "block_1");
        assert_eq!(blocks[0].start, 101);
        assert_eq!(blocks[0].end, 141);

        assert_eq!(blocks[1].block_id, "block_2");
        assert_eq!(blocks[1].start, 221);
        assert_eq!(blocks[1].end, 241);
    }

    #[test]
    fn zero_only_calls_produce_no_blocks() {
        let calls = vec![
            LocalPhaseCall {
                pos: 101,
                phase_state: 0,
                confidence: 0.5,
            },
            LocalPhaseCall {
                pos: 121,
                phase_state: 0,
                confidence: 0.5,
            },
        ];

        let blocks = make_blocks_from_local_calls(&calls);
        assert!(blocks.is_empty());
    }

    #[test]
    fn trailing_informative_block_is_closed() {
        let calls = vec![
            LocalPhaseCall {
                pos: 101,
                phase_state: 0,
                confidence: 0.5,
            },
            LocalPhaseCall {
                pos: 121,
                phase_state: 1,
                confidence: 0.8,
            },
            LocalPhaseCall {
                pos: 141,
                phase_state: -1,
                confidence: 0.7,
            },
        ];

        let blocks = make_blocks_from_local_calls(&calls);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].start, 121);
        assert_eq!(blocks[0].end, 141);
    }
}
