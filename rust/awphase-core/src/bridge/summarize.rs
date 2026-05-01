use crate::bridge::block_summary::BlockSummary;
use crate::types::LocalPhaseCall;

#[derive(Debug, Clone)]
pub struct BlockRange {
    pub block_id: String,
    pub start: u64,
    pub end: u64,
}

pub fn summarize_blocks_from_local_calls(
    calls: &[LocalPhaseCall],
    ranges: &[BlockRange],
) -> Vec<BlockSummary> {
    ranges
        .iter()
        .map(|r| {
            let in_block: Vec<&LocalPhaseCall> = calls
                .iter()
                .filter(|c| c.pos >= r.start && c.pos <= r.end)
                .collect();

            let informative: Vec<&LocalPhaseCall> = in_block
                .iter()
                .copied()
                .filter(|c| c.phase_state != 0)
                .collect();

            let abstained_sites = in_block.iter().filter(|c| c.phase_state == 0).count();
            let informative_sites = informative.len();

            let signed_confidence_sum: f32 = informative
                .iter()
                .map(|c| c.phase_state as f32 * c.confidence)
                .sum();

            let support_weight: f32 = informative.iter().map(|c| c.confidence).sum();

            let orientation_bias = if support_weight > 0.0 {
                (signed_confidence_sum / support_weight).clamp(-1.0, 1.0)
            } else {
                0.0
            };

            let mean_confidence = if informative_sites > 0 {
                support_weight / informative_sites as f32
            } else {
                0.0
            };

            let coherence = if support_weight > 0.0 {
                (signed_confidence_sum.abs() / support_weight).clamp(0.0, 1.0)
            } else {
                0.0
            };

            BlockSummary {
                block_id: r.block_id.clone(),
                orientation_bias,
                support_weight,
                informative_sites,
                abstained_sites,
                mean_confidence,
                signed_confidence_sum,
                coherence,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::LocalPhaseCall;

    #[test]
    fn summarizes_blocks_from_local_calls_correctly() {
        let calls = vec![
            LocalPhaseCall {
                pos: 101,
                phase_state: 1,
                confidence: 0.8,
            },
            LocalPhaseCall {
                pos: 121,
                phase_state: 0,
                confidence: 0.5,
            },
            LocalPhaseCall {
                pos: 251,
                phase_state: -1,
                confidence: 0.7,
            },
            LocalPhaseCall {
                pos: 401,
                phase_state: 0,
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
                end: 300,
            },
            BlockRange {
                block_id: "block_3".to_string(),
                start: 301,
                end: 450,
            },
        ];

        let blocks = summarize_blocks_from_local_calls(&calls, &ranges);
        assert_eq!(blocks.len(), 3);

        assert_eq!(blocks[0].block_id, "block_1");
        assert!(blocks[0].orientation_bias > 0.0);
        assert_eq!(blocks[0].informative_sites, 1);
        assert_eq!(blocks[0].abstained_sites, 1);
        assert!(blocks[0].coherence > 0.0);

        assert_eq!(blocks[1].block_id, "block_2");
        assert!(blocks[1].orientation_bias < 0.0);
        assert_eq!(blocks[1].informative_sites, 1);

        assert_eq!(blocks[2].block_id, "block_3");
        assert_eq!(blocks[2].orientation_bias, 0.0);
        assert_eq!(blocks[2].informative_sites, 0);
        assert_eq!(blocks[2].abstained_sites, 1);
        assert_eq!(blocks[2].coherence, 0.0);
    }
}
