export interface FastaFileInfo {
  path: string;
  filename: string;
  num_sequences: number;
  min_length: number;
  max_length: number;
  avg_length: number;
  first_ids: string[];
}

export interface OligoInfo {
  id: string;
  seq: string;
}

export interface OligoSetTexts {
  forward_text: string;
  reverse_text: string;
  probe_text: string;
}

export interface AlignmentSettings {
  mode: "Local" | "Global";
  match_score: number;
  mismatch_score: number;
  gap_open_score: number;
  gap_extend_score: number;
  min_fwd_matched: number;
  min_rev_matched: number;
  min_probe_matched: number;
  min_coverage: number;
  max_mismatches_per_oligo: number;
  ambiguity_display: "ShowBases" | "ShowDots";
  min_amplicon_size: number | null;
  max_amplicon_size: number | null;
}

export interface AnalysisRequest {
  forward_text: string;
  reverse_text: string;
  probe_text: string;
  settings: AlignmentSettings;
}

export interface ProgressPayload {
  current: number;
  total: number;
  status: string;
}

export interface AnalysisResultsPayload {
  output_text: string;
  total_sequences: number;
  sequences_with_min_matches: number;
  has_results: boolean;
}
