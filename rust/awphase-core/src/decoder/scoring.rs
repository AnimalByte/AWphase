use crate::types::ReadObservation;

pub fn score_site_from_reads(site_pos: u64, reads: &[ReadObservation]) -> (f32, f32, f32) {
    let mut alt_support = 0.0_f32;
    let mut ref_support = 0.0_f32;

    let mut obs_per_read = std::collections::HashMap::<&str, usize>::new();
    for r in reads.iter() {
        *obs_per_read.entry(r.read_id.as_str()).or_insert(0) += 1;
    }

    for read in reads.iter().filter(|r| r.site_pos == site_pos) {
        let linked = obs_per_read
            .get(read.read_id.as_str())
            .copied()
            .unwrap_or(0)
            >= 2;
        if !linked {
            continue;
        }

        if read.is_ambiguous {
            continue;
        }

        if let Some(delta) = read.allele_score_delta {
            if delta.abs() >= 0.75 {
                let conflict = (read.allele > 0 && delta < 0.0) || (read.allele < 0 && delta > 0.0);
                if conflict {
                    continue;
                }
            }
        }

        let mut weight = read.weight_v2.unwrap_or_else(|| {
            ((read.baseq as f32) / 40.0).clamp(0.1, 1.0)
                * ((read.mapq as f32) / 60.0).clamp(0.1, 1.0)
        });

        if let Some(conf) = read.allele_confidence {
            weight *= 1.0 + 0.20 * conf.clamp(0.0, 1.0);
        }

        if let Some(delta) = read.allele_score_delta {
            let strong_support =
                (read.allele > 0 && delta > 0.75) || (read.allele < 0 && delta < -0.75);
            if strong_support {
                weight *= 1.05;
            }
        }

        match read.allele {
            1 => alt_support += weight,
            -1 => ref_support += weight,
            _ => {}
        }
    }

    let total = alt_support + ref_support;
    let read_conf = if total > 0.0 {
        alt_support.max(ref_support) / total
    } else {
        0.0
    };

    (alt_support, ref_support, read_conf)
}

pub fn combine_read_and_donor(
    alt_support: f32,
    ref_support: f32,
    donor_bias: f32,
    donor_weight: f32,
) -> (i8, f32) {
    let signed_read_score = alt_support - ref_support;
    let combined = signed_read_score + donor_weight * donor_bias;

    let total = alt_support + ref_support;

    if total <= 0.0 {
        if donor_bias > 0.2 {
            return (1, donor_bias.abs().min(1.0) * 0.5);
        } else if donor_bias < -0.2 {
            return (-1, donor_bias.abs().min(1.0) * 0.5);
        } else {
            return (0, 0.0);
        }
    }

    let confidence = ((combined.abs()) / total).clamp(0.0, 1.0);

    if combined > 0.05 {
        (1, confidence.max(0.5))
    } else if combined < -0.05 {
        (-1, confidence.max(0.5))
    } else {
        (0, 0.5)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_scoring_prefers_alt_when_alt_support_is_higher() {
        let reads = vec![
            ReadObservation {
                read_id: "r1".to_string(),
                site_pos: 101,
                allele: 1,
                baseq: 35,
                mapq: 60,
                ..Default::default()
            },
            ReadObservation {
                read_id: "r1".to_string(),
                site_pos: 102,
                allele: 1,
                baseq: 35,
                mapq: 60,
                ..Default::default()
            },
            ReadObservation {
                read_id: "r2".to_string(),
                site_pos: 101,
                allele: 1,
                baseq: 32,
                mapq: 60,
                ..Default::default()
            },
            ReadObservation {
                read_id: "r2".to_string(),
                site_pos: 102,
                allele: 1,
                baseq: 32,
                mapq: 60,
                ..Default::default()
            },
            ReadObservation {
                read_id: "r3".to_string(),
                site_pos: 101,
                allele: -1,
                baseq: 20,
                mapq: 45,
                ..Default::default()
            },
            ReadObservation {
                read_id: "r3".to_string(),
                site_pos: 102,
                allele: -1,
                baseq: 20,
                mapq: 45,
                ..Default::default()
            },
        ];

        let (alt, ref_, conf) = score_site_from_reads(101, &reads);
        assert!(alt > ref_);
        assert!(conf > 0.5);
    }

    #[test]
    fn tie_with_no_donor_bias_abstains() {
        let (phase, conf) = combine_read_and_donor(1.0, 1.0, 0.0, 0.15);
        assert_eq!(phase, 0);
        assert_eq!(conf, 0.5);
    }

    #[test]
    fn donor_bias_breaks_tie_toward_plus_one() {
        let (phase, conf) = combine_read_and_donor(1.0, 1.0, 0.8, 0.15);
        assert_eq!(phase, 1);
        assert!(conf >= 0.5);
    }

    #[test]
    fn donor_bias_breaks_tie_toward_minus_one() {
        let (phase, conf) = combine_read_and_donor(1.0, 1.0, -0.8, 0.15);
        assert_eq!(phase, -1);
        assert!(conf >= 0.5);
    }

    #[test]
    fn no_reads_and_no_donor_bias_abstains() {
        let (phase, conf) = combine_read_and_donor(0.0, 0.0, 0.0, 0.15);
        assert_eq!(phase, 0);
        assert_eq!(conf, 0.0);
    }
}
