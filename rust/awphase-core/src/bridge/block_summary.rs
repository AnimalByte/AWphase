use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlockSummary {
    pub block_id: String,
    pub orientation_bias: f32,
    pub support_weight: f32,
    pub informative_sites: usize,
    pub abstained_sites: usize,
    pub mean_confidence: f32,
    pub signed_confidence_sum: f32,
    pub coherence: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryScoreInput {
    pub left_block: String,
    pub right_block: String,
    pub same_score: f32,
    pub flip_score: f32,
    pub informative_sites: usize,
}

fn block_reliability(block: &BlockSummary) -> f32 {
    let bias_strength = block.orientation_bias.abs().clamp(0.0, 1.0);
    let support_term = (block.support_weight / 4.0).clamp(0.0, 1.0);
    let info_term = (block.informative_sites as f32 / 4.0).clamp(0.0, 1.0);
    let confidence_term = block.mean_confidence.clamp(0.0, 1.0);
    let coherence_term = block.coherence.clamp(0.0, 1.0);

    let total_sites = (block.informative_sites + block.abstained_sites).max(1) as f32;
    let abstention_rate = block.abstained_sites as f32 / total_sites;
    let abstention_penalty = 1.0 - abstention_rate.clamp(0.0, 1.0);

    (0.22 * bias_strength
        + 0.18 * support_term
        + 0.16 * info_term
        + 0.16 * confidence_term
        + 0.18 * coherence_term
        + 0.10 * abstention_penalty)
        .clamp(0.0, 1.0)
}

pub fn boundary_scores_from_blocks(
    left: &BlockSummary,
    right: &BlockSummary,
) -> BoundaryScoreInput {
    let left_rel = block_reliability(left);
    let right_rel = block_reliability(right);

    let pair_reliability = (left_rel * right_rel).sqrt().clamp(0.0, 1.0);
    let directional_agreement = (left.orientation_bias * right.orientation_bias).clamp(-1.0, 1.0);

    let same_evidence = pair_reliability * directional_agreement.max(0.0);
    let flip_evidence = pair_reliability * (-directional_agreement).max(0.0);

    let uncertainty_mass = (1.0 - pair_reliability).max(0.0);

    let same_score = same_evidence + 0.5 * uncertainty_mass;
    let flip_score = flip_evidence + 0.5 * uncertainty_mass;

    let informative_sites = left.informative_sites.min(right.informative_sites);

    BoundaryScoreInput {
        left_block: left.block_id.clone(),
        right_block: right.block_id.clone(),
        same_score,
        flip_score,
        informative_sites,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positive_bias_pair_favors_same() {
        let a = BlockSummary {
            block_id: "b1".to_string(),
            orientation_bias: 0.9,
            support_weight: 3.2,
            informative_sites: 4,
            abstained_sites: 0,
            mean_confidence: 0.8,
            signed_confidence_sum: 3.2,
            coherence: 1.0,
        };
        let b = BlockSummary {
            block_id: "b2".to_string(),
            orientation_bias: 0.8,
            support_weight: 3.0,
            informative_sites: 4,
            abstained_sites: 0,
            mean_confidence: 0.75,
            signed_confidence_sum: 3.0,
            coherence: 1.0,
        };

        let out = boundary_scores_from_blocks(&a, &b);
        assert!(out.same_score > out.flip_score);
    }

    #[test]
    fn opposite_bias_pair_favors_flip() {
        let a = BlockSummary {
            block_id: "b1".to_string(),
            orientation_bias: 0.9,
            support_weight: 3.2,
            informative_sites: 4,
            abstained_sites: 0,
            mean_confidence: 0.8,
            signed_confidence_sum: 3.2,
            coherence: 1.0,
        };
        let b = BlockSummary {
            block_id: "b2".to_string(),
            orientation_bias: -0.8,
            support_weight: 3.0,
            informative_sites: 4,
            abstained_sites: 0,
            mean_confidence: 0.75,
            signed_confidence_sum: -3.0,
            coherence: 1.0,
        };

        let out = boundary_scores_from_blocks(&a, &b);
        assert!(out.flip_score > out.same_score);
    }

    #[test]
    fn weak_or_uncertain_pair_stays_close() {
        let a = BlockSummary {
            block_id: "b1".to_string(),
            orientation_bias: 0.1,
            support_weight: 0.5,
            informative_sites: 1,
            abstained_sites: 3,
            mean_confidence: 0.5,
            signed_confidence_sum: 0.5,
            coherence: 0.3,
        };
        let b = BlockSummary {
            block_id: "b2".to_string(),
            orientation_bias: 0.05,
            support_weight: 0.5,
            informative_sites: 1,
            abstained_sites: 3,
            mean_confidence: 0.5,
            signed_confidence_sum: 0.5,
            coherence: 0.3,
        };

        let out = boundary_scores_from_blocks(&a, &b);
        assert!((out.same_score - out.flip_score).abs() < 0.2);
    }
}
