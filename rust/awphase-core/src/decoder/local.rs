use anyhow::Result;

use crate::candidate::interface::{CandidateSet, PanelHaplotype};
use crate::decoder::debug::{LocalDecodeDebug, SiteDecodeDebug};
use crate::decoder::read_selection::{select_informative_read_observations, ReadSelectionConfig};
use crate::decoder::scoring::{combine_read_and_donor, score_site_from_reads};
use crate::decoder::window_hmm::{
    build_candidate_window_states, normalized_state_weights, posterior_site_margins,
    viterbi_decode_windows, CandidateWindowState,
};
use crate::types::{
    AncestryPrior, LocalDecodeResult, LocalPhaseCall, ReadObservation, SiteDonorPriors,
    VariantSite, WindowId,
};

pub trait LocalDecoder {
    fn decode(
        &self,
        window: &WindowId,
        variants: &[VariantSite],
        reads: &[ReadObservation],
        candidates: &CandidateSet,
        ancestry: &AncestryPrior,
        site_priors: &SiteDonorPriors,
        panel_haplotypes: &[PanelHaplotype],
        phase_hint: Option<&[i8]>,
    ) -> Result<LocalDecodeResult>;

    fn decode_with_debug(
        &self,
        window: &WindowId,
        variants: &[VariantSite],
        reads: &[ReadObservation],
        candidates: &CandidateSet,
        ancestry: &AncestryPrior,
        site_priors: &SiteDonorPriors,
        panel_haplotypes: &[PanelHaplotype],
        phase_hint: Option<&[i8]>,
    ) -> Result<(LocalDecodeResult, LocalDecodeDebug)>;
}

pub struct ReadDonorLocalDecoder;

fn signf(x: f32) -> i8 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

fn normalized_read_bias(alt_support: f32, ref_support: f32) -> f32 {
    let denom = alt_support.abs() + ref_support.abs();
    if denom <= 1e-6 {
        0.0
    } else {
        ((alt_support - ref_support) / denom).clamp(-1.0, 1.0)
    }
}

impl LocalDecoder for ReadDonorLocalDecoder {
    fn decode(
        &self,
        window: &WindowId,
        variants: &[VariantSite],
        reads: &[ReadObservation],
        candidates: &CandidateSet,
        ancestry: &AncestryPrior,
        site_priors: &SiteDonorPriors,
        panel_haplotypes: &[PanelHaplotype],
        phase_hint: Option<&[i8]>,
    ) -> Result<LocalDecodeResult> {
        let (result, _) = self.decode_with_debug(
            window,
            variants,
            reads,
            candidates,
            ancestry,
            site_priors,
            panel_haplotypes,
            phase_hint,
        )?;
        Ok(result)
    }

    fn decode_with_debug(
        &self,
        window: &WindowId,
        variants: &[VariantSite],
        reads: &[ReadObservation],
        candidates: &CandidateSet,
        _ancestry: &AncestryPrior,
        _site_priors: &SiteDonorPriors,
        panel_haplotypes: &[PanelHaplotype],
        phase_hint: Option<&[i8]>,
    ) -> Result<(LocalDecodeResult, LocalDecodeDebug)> {
        if phase_hint.is_none() {
            let mut calls = Vec::with_capacity(variants.len());
            let mut ambiguity_flags = Vec::new();
            let mut debug_sites = Vec::with_capacity(variants.len());

            for v in variants {
                let (alt_support, ref_support, _) = score_site_from_reads(v.pos, reads);
                let read_bias = normalized_read_bias(alt_support, ref_support);
                let read_sign = signf(read_bias);

                let (phase_state, confidence) =
                    combine_read_and_donor(alt_support, ref_support, 0.0, 0.20);

                if phase_state == 0 {
                    ambiguity_flags.push(v.pos);
                }

                calls.push(LocalPhaseCall {
                    pos: v.pos,
                    phase_state,
                    confidence,
                });

                debug_sites.push(SiteDecodeDebug {
                    pos: v.pos,
                    window_index: 0,
                    candidate_states: vec![],
                    viterbi_best_state: None,
                    alt_support,
                    ref_support,
                    read_bias,
                    donor_margin: 0.0,
                    combined_donor_bias: 0.0,
                    path_sign: 0,
                    read_sign,
                    donor_sign: 0,
                    path_margin_agrees: false,
                    donor_read_conflict: false,
                    weak_read_evidence: read_bias.abs() < 0.20,
                    forced_abstain: false,
                    decision_source: "read_only".to_string(),
                    phase_state,
                    confidence,
                });
            }

            return Ok((
                LocalDecodeResult {
                    window: window.clone(),
                    calls,
                    ambiguity_flags,
                },
                LocalDecodeDebug { sites: debug_sites },
            ));
        }

        let hint = phase_hint.unwrap();

        let (reads, _read_selection_summary) = select_informative_read_observations(
            variants,
            reads,
            &ReadSelectionConfig {
                max_reads: 20000,
                max_per_site: 12,
                min_sites_per_read: 3,
            },
        );

        let donor_weight_default = 0.20_f32;
        let top_k_states = 4usize;
        let lambda_transition = 0.35_f32;
        let sites_per_window = 4usize;

        let positions: Vec<u64> = variants.iter().map(|v| v.pos).collect();

        let mut window_ranges = Vec::new();
        let mut start = 0usize;
        while start < variants.len() {
            let end = (start + sites_per_window).min(variants.len());
            window_ranges.push((start, end));
            start = end;
        }

        let state_windows: Vec<Vec<CandidateWindowState>> = window_ranges
            .iter()
            .map(|(s, e)| {
                build_candidate_window_states(
                    candidates,
                    panel_haplotypes,
                    hint,
                    &positions,
                    &reads,
                    *s,
                    *e,
                    top_k_states,
                )
            })
            .collect();

        let viterbi = viterbi_decode_windows(&state_windows, lambda_transition);

        let mut calls = Vec::with_capacity(variants.len());
        let mut ambiguity_flags = Vec::new();
        let mut debug_sites = Vec::with_capacity(variants.len());

        for (window_index, (start_idx, end_idx)) in window_ranges.iter().enumerate() {
            let states = &state_windows[window_index];
            let scores = &viterbi.end_scores[window_index];
            let weights = normalized_state_weights(scores);
            let site_margins = posterior_site_margins(states, &weights);
            let best_state = viterbi.best_path[window_index].clone();

            for idx in *start_idx..*end_idx {
                let local_i = idx - *start_idx;
                let beam_margin = site_margins.get(local_i).copied().unwrap_or(0.0);

                let path_sign = match &best_state {
                    Some(s) if local_i < s.alleles.len() && s.alleles[local_i] > 0 => 1,
                    Some(s) if local_i < s.alleles.len() && s.alleles[local_i] < 0 => -1,
                    _ => 0,
                };

                let donor_bias = if path_sign > 0 {
                    beam_margin.abs()
                } else if path_sign < 0 {
                    -beam_margin.abs()
                } else {
                    beam_margin
                };

                let donor_sign = signf(donor_bias);
                let beam_sign = signf(beam_margin);
                let path_margin_agrees = path_sign != 0 && path_sign == beam_sign;

                let (alt_support, ref_support, _) =
                    score_site_from_reads(variants[idx].pos, &reads);
                let read_bias = normalized_read_bias(alt_support, ref_support);
                let read_sign = signf(read_bias);

                let strong_donor = donor_bias.abs() >= 0.75;
                let weak_read_evidence = read_bias.abs() < 0.20;
                let donor_read_conflict =
                    donor_sign != 0 && read_sign != 0 && donor_sign != read_sign;

                let combined_donor_bias = if donor_read_conflict {
                    donor_bias * 0.35
                } else if strong_donor && weak_read_evidence {
                    donor_bias * 0.75
                } else {
                    donor_bias
                };

                let donor_weight = if donor_read_conflict {
                    0.08_f32
                } else if strong_donor && weak_read_evidence {
                    0.15_f32
                } else {
                    donor_weight_default
                };

                let (mut phase_state, mut confidence) = combine_read_and_donor(
                    alt_support,
                    ref_support,
                    combined_donor_bias,
                    donor_weight,
                );

                let forced_abstain = donor_read_conflict && weak_read_evidence;

                let decision_source = if forced_abstain {
                    phase_state = 0;
                    confidence = confidence.min(0.30);
                    "forced_abstain_read_conflict_and_weak_read".to_string()
                } else if donor_read_conflict {
                    confidence = confidence.min(0.49);
                    "path_posterior_read_conflict".to_string()
                } else if strong_donor && weak_read_evidence {
                    confidence = confidence.min(0.70);
                    "path_posterior_agree_strong_weak_read".to_string()
                } else if path_sign == 0 {
                    "posterior_only".to_string()
                } else if path_margin_agrees && donor_bias.abs() >= 0.25 {
                    "path_and_posterior_agree_strong".to_string()
                } else if path_margin_agrees {
                    "path_and_posterior_agree_weak".to_string()
                } else {
                    "path_posterior_conflict".to_string()
                };

                if phase_state == 0 {
                    ambiguity_flags.push(variants[idx].pos);
                }

                calls.push(LocalPhaseCall {
                    pos: variants[idx].pos,
                    phase_state,
                    confidence,
                });

                debug_sites.push(SiteDecodeDebug {
                    pos: variants[idx].pos,
                    window_index,
                    candidate_states: states.clone(),
                    viterbi_best_state: best_state.clone(),
                    alt_support,
                    ref_support,
                    read_bias,
                    donor_margin: donor_bias,
                    combined_donor_bias,
                    path_sign,
                    read_sign,
                    donor_sign,
                    path_margin_agrees,
                    donor_read_conflict,
                    weak_read_evidence,
                    forced_abstain,
                    decision_source,
                    phase_state,
                    confidence,
                });
            }
        }

        Ok((
            LocalDecodeResult {
                window: window.clone(),
                calls,
                ambiguity_flags,
            },
            LocalDecodeDebug { sites: debug_sites },
        ))
    }
}
