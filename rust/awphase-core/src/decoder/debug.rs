use serde::{Deserialize, Serialize};

use crate::decoder::window_hmm::CandidateWindowState;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SiteDecodeDebug {
    pub pos: u64,
    pub window_index: usize,
    pub candidate_states: Vec<CandidateWindowState>,
    pub viterbi_best_state: Option<CandidateWindowState>,

    pub alt_support: f32,
    pub ref_support: f32,
    pub read_bias: f32,

    pub donor_margin: f32,
    pub combined_donor_bias: f32,

    pub path_sign: i8,
    pub read_sign: i8,
    pub donor_sign: i8,

    pub path_margin_agrees: bool,
    pub donor_read_conflict: bool,
    pub weak_read_evidence: bool,
    pub forced_abstain: bool,

    pub decision_source: String,
    pub phase_state: i8,
    pub confidence: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct LocalDecodeDebug {
    pub sites: Vec<SiteDecodeDebug>,
}
