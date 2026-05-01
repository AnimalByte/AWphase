use anyhow::Result;
use serde::{Deserialize, Serialize};

use crate::candidate::interface::{CandidateEngine, CandidateSet, PanelHaplotype};
use crate::decoder::debug::LocalDecodeDebug;
use crate::decoder::local::LocalDecoder;
use crate::types::{
    AncestryPrior, LocalDecodeResult, ReadObservation, SiteDonorPriors, VariantSite, WindowId,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TwoPassDecodeResult {
    pub initial_guess: LocalDecodeResult,
    pub candidates: CandidateSet,
    pub final_result: LocalDecodeResult,
    pub final_debug: LocalDecodeDebug,
}

pub fn run_two_pass_local_decode<E: CandidateEngine, D: LocalDecoder>(
    engine: &E,
    decoder: &D,
    window: &WindowId,
    variants: &[VariantSite],
    reads: &[ReadObservation],
    ancestry: &AncestryPrior,
    site_priors: &SiteDonorPriors,
    panel_haplotypes: &[PanelHaplotype],
) -> Result<TwoPassDecodeResult> {
    let empty_candidates = CandidateSet {
        window: window.clone(),
        donors: vec![],
        composites: vec![],
    };

    let initial_guess = decoder.decode(
        window,
        variants,
        reads,
        &empty_candidates,
        ancestry,
        site_priors,
        panel_haplotypes,
        None,
    )?;

    let phase_hint: Vec<i8> = initial_guess.calls.iter().map(|c| c.phase_state).collect();

    let candidates = engine.retrieve(window, variants, &phase_hint, ancestry)?;

    let (final_result, final_debug) = decoder.decode_with_debug(
        window,
        variants,
        reads,
        &candidates,
        ancestry,
        site_priors,
        panel_haplotypes,
        Some(&phase_hint),
    )?;

    Ok(TwoPassDecodeResult {
        initial_guess,
        candidates,
        final_result,
        final_debug,
    })
}
