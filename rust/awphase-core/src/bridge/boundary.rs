use crate::confidence::interface::{BoundaryDecision, ConfidenceScore, Decision};

pub fn score_boundary(
    left_block_id: &str,
    right_block_id: &str,
    same_score: f32,
    flip_score: f32,
    min_confidence: f32,
    min_info: usize,
    informative_sites: usize,
) -> BoundaryDecision {
    if informative_sites < min_info {
        return BoundaryDecision {
            left_block: left_block_id.to_string(),
            right_block: right_block_id.to_string(),
            decision: Decision::Abstain,
            confidence: ConfidenceScore {
                value: 0.0,
                informative_sites,
                reason: Some("insufficient_informative_sites".to_string()),
            },
        };
    }

    let total = same_score.abs() + flip_score.abs();
    let confidence = if total > 0.0 {
        (same_score - flip_score).abs() / total
    } else {
        0.0
    };

    if confidence < min_confidence {
        return BoundaryDecision {
            left_block: left_block_id.to_string(),
            right_block: right_block_id.to_string(),
            decision: Decision::Abstain,
            confidence: ConfidenceScore {
                value: confidence,
                informative_sites,
                reason: Some("low_confidence".to_string()),
            },
        };
    }

    // No-flip ablation:
    // if flip evidence wins, abstain instead of flipping downstream orientation.
    if same_score >= flip_score {
        BoundaryDecision {
            left_block: left_block_id.to_string(),
            right_block: right_block_id.to_string(),
            decision: Decision::Join,
            confidence: ConfidenceScore {
                value: confidence,
                informative_sites,
                reason: None,
            },
        }
    } else {
        BoundaryDecision {
            left_block: left_block_id.to_string(),
            right_block: right_block_id.to_string(),
            decision: Decision::Abstain,
            confidence: ConfidenceScore {
                value: confidence,
                informative_sites,
                reason: Some("flip_suppressed_no_flip_ablation".to_string()),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::confidence::interface::Decision;

    #[test]
    fn abstains_when_informative_sites_too_low() {
        let d = score_boundary("b1", "b2", 2.0, 1.0, 0.7, 4, 2);
        assert!(matches!(d.decision, Decision::Abstain));
    }

    #[test]
    fn abstains_when_confidence_too_low() {
        let d = score_boundary("b1", "b2", 1.0, 0.9, 0.7, 2, 5);
        assert!(matches!(d.decision, Decision::Abstain));
    }

    #[test]
    fn chooses_join_when_same_wins_confidently() {
        let d = score_boundary("b1", "b2", 5.0, 1.0, 0.5, 2, 5);
        assert!(matches!(d.decision, Decision::Join));
    }

    #[test]
    fn suppresses_flip_when_flip_wins_confidently() {
        let d = score_boundary("b1", "b2", 1.0, 5.0, 0.5, 2, 5);
        assert!(matches!(d.decision, Decision::Abstain));
        assert_eq!(
            d.confidence.reason.as_deref(),
            Some("flip_suppressed_no_flip_ablation")
        );
    }
}
