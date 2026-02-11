use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::{Arc, Mutex};

// ============================================================================
// Core Data Structures (from reference)
// ============================================================================

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FastaRecord {
    pub id: String,
    pub seq: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Oligo {
    pub id: String,
    pub seq: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct OligoSerde {
    pub id: String,
    pub seq: String,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct OligoSet {
    pub forward_primers: Vec<OligoSerde>,
    pub reverse_primers: Vec<OligoSerde>,
    pub probes: Vec<OligoSerde>,
}

#[derive(Clone, Debug, Copy, PartialEq)]
pub enum OligoCategory {
    ForwardPrimer,
    ReversePrimer,
    Probe,
}

#[derive(Clone, Debug)]
pub struct OligoResult {
    pub matched: bool,
    pub orientation: Option<Orientation>,
    pub signature: String,
    pub mismatches: usize,
    pub score: i32,
    pub coverage: f64,
    pub start_pos: Option<usize>,
    pub end_pos: Option<usize>,
}

impl OligoResult {
    pub fn no_match() -> Self {
        Self {
            matched: false,
            orientation: None,
            signature: String::new(),
            mismatches: 0,
            score: 0,
            coverage: 0.0,
            start_pos: None,
            end_pos: None,
        }
    }
}

#[derive(Clone, Debug, Copy, PartialEq)]
pub enum Orientation {
    Sense,
    Antisense,
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Orientation::Sense => write!(f, "fwd"),
            Orientation::Antisense => write!(f, "rev"),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct AmpliconInfo {
    pub found: bool,
    pub forward_oligo_id: Option<String>,
    pub reverse_oligo_id: Option<String>,
    pub start: usize,
    pub end: usize,
    pub size: usize,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct PatternData {
    pub count: usize,
    pub total_mismatches: usize,
    pub matched_fwd: usize,
    pub matched_rev: usize,
    pub matched_probe: usize,
    pub examples: Vec<String>,
    pub amplicon_lengths: Vec<usize>,
}

#[derive(Clone, Debug, Default)]
pub struct OligoStats {
    pub total_matches: usize,
    pub sense_matches: usize,
    pub antisense_matches: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AlignmentSettings {
    pub mode: AlignmentMode,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open_score: i32,
    pub gap_extend_score: i32,
    pub min_fwd_matched: usize,
    pub min_rev_matched: usize,
    pub min_probe_matched: usize,
    pub min_coverage: f64,
    pub max_mismatches_per_oligo: usize,
    pub ambiguity_display: AmbiguityDisplayMode,
    pub min_amplicon_size: Option<usize>,
    pub max_amplicon_size: Option<usize>,
}

impl Default for AlignmentSettings {
    fn default() -> Self {
        Self {
            mode: AlignmentMode::Local,
            match_score: 2,
            mismatch_score: -1,
            gap_open_score: -2,
            gap_extend_score: -1,
            min_fwd_matched: 1,
            min_rev_matched: 1,
            min_probe_matched: 0,
            min_coverage: 0.8,
            max_mismatches_per_oligo: 7,
            ambiguity_display: AmbiguityDisplayMode::ShowDots,
            min_amplicon_size: None,
            max_amplicon_size: None,
        }
    }
}

#[derive(Clone, Debug, Copy, PartialEq, Serialize, Deserialize)]
pub enum AlignmentMode {
    Local,
    Global,
}

#[derive(Clone, Debug, Copy, PartialEq, Serialize, Deserialize)]
pub enum AmbiguityDisplayMode {
    ShowBases,
    ShowDots,
}

#[derive(Clone, Debug)]
pub struct SequenceResult {
    pub fwd_results: HashMap<String, OligoResult>,
    pub rev_results: HashMap<String, OligoResult>,
    pub probe_results: HashMap<String, OligoResult>,
    pub fwd_matched: usize,
    pub rev_matched: usize,
    pub probe_matched: usize,
    pub total_mismatches: usize,
    pub amplicon_info: Option<AmpliconInfo>,
    pub best_fwd_mm: Option<usize>,
    pub best_rev_mm: Option<usize>,
    pub best_probe_mm: Option<usize>,
}

#[derive(Clone, Debug, Default)]
pub struct MismatchDistribution {
    pub zero_mm: usize,
    pub one_mm: usize,
    pub more_mm: usize,
    pub no_match: usize,
}

#[derive(Clone, Debug, Default)]
pub struct AnalysisResults {
    pub alignment_dict: HashMap<String, PatternData>,
    pub oligo_stats: HashMap<String, OligoStats>,
    pub total_sequences: usize,
    pub sequences_with_min_matches: usize,
    pub sequences_with_valid_amplicon: usize,
    pub sequences_failed_amplicon: usize,
    pub fwd_mm_dist: MismatchDistribution,
    pub rev_mm_dist: MismatchDistribution,
    pub probe_mm_dist: MismatchDistribution,
    pub overall_all_perfect: usize,
    pub overall_max_one_mm: usize,
    pub overall_two_plus_mm: usize,
    pub overall_no_match: usize,
    pub output_text: String,
}

#[derive(Clone)]
pub struct ProgressTracker {
    pub current: Arc<AtomicUsize>,
    pub total: Arc<AtomicUsize>,
    pub status: Arc<Mutex<String>>,
    pub running: Arc<AtomicBool>,
}

impl Default for ProgressTracker {
    fn default() -> Self {
        Self {
            current: Arc::new(AtomicUsize::new(0)),
            total: Arc::new(AtomicUsize::new(0)),
            status: Arc::new(Mutex::new("Ready".to_string())),
            running: Arc::new(AtomicBool::new(false)),
        }
    }
}

// ============================================================================
// Frontend-facing payload types
// ============================================================================

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FastaFileInfo {
    pub path: String,
    pub filename: String,
    pub num_sequences: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub avg_length: f64,
    pub first_ids: Vec<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct OligoInfo {
    pub id: String,
    pub seq: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct OligoSetTexts {
    pub forward_text: String,
    pub reverse_text: String,
    pub probe_text: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AnalysisRequest {
    pub forward_text: String,
    pub reverse_text: String,
    pub probe_text: String,
    pub settings: AlignmentSettings,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ProgressPayload {
    pub current: usize,
    pub total: usize,
    pub status: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AnalysisResultsPayload {
    pub output_text: String,
    pub total_sequences: usize,
    pub sequences_with_min_matches: usize,
    pub has_results: bool,
}
