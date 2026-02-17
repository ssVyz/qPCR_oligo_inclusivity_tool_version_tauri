use bio::alignment::pairwise::{Aligner, Scoring};
use bio::alignment::AlignmentOperation;
use bio::io::fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use crate::models::*;

// ============================================================================
// Bioinformatics Functions
// ============================================================================

pub fn iupac_match(a: u8, b: u8) -> bool {
    let a_upper = a.to_ascii_uppercase();
    let b_upper = b.to_ascii_uppercase();

    if a_upper == b'N' || b_upper == b'N' {
        return false;
    }

    if a_upper == b_upper {
        return true;
    }

    let get_bases = |code: u8| -> &'static [u8] {
        match code.to_ascii_uppercase() {
            b'A' => &[b'A'],
            b'T' => &[b'T'],
            b'G' => &[b'G'],
            b'C' => &[b'C'],
            b'U' => &[b'U'],
            b'R' => &[b'A', b'G'],
            b'Y' => &[b'C', b'T'],
            b'S' => &[b'G', b'C'],
            b'W' => &[b'A', b'T'],
            b'K' => &[b'G', b'T'],
            b'M' => &[b'A', b'C'],
            b'B' => &[b'C', b'G', b'T'],
            b'D' => &[b'A', b'G', b'T'],
            b'H' => &[b'A', b'C', b'T'],
            b'V' => &[b'A', b'C', b'G'],
            _ => &[],
        }
    };

    let a_bases = get_bases(a_upper);
    let b_bases = get_bases(b_upper);

    for &base_a in a_bases {
        for &base_b in b_bases {
            let normalized_a = if base_a == b'U' { b'T' } else { base_a };
            let normalized_b = if base_b == b'U' { b'T' } else { base_b };
            if normalized_a == normalized_b {
                return true;
            }
        }
    }

    false
}

pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'U' => 'A',
            'R' => 'Y',
            'Y' => 'R',
            'S' => 'S',
            'W' => 'W',
            'K' => 'M',
            'M' => 'K',
            'B' => 'V',
            'D' => 'H',
            'H' => 'D',
            'V' => 'B',
            'N' => 'N',
            _ => 'N',
        })
        .collect()
}

pub fn parse_fasta(path: &PathBuf) -> Result<Vec<FastaRecord>, String> {
    let file = File::open(path).map_err(|e| format!("Failed to open file: {}", e))?;
    let reader = fasta::Reader::new(BufReader::new(file));

    let mut records = Vec::new();
    for result in reader.records() {
        match result {
            Ok(record) => {
                let id = record.id().to_string();
                let seq = std::str::from_utf8(record.seq())
                    .map_err(|e| format!("Invalid UTF-8 in sequence: {}", e))?
                    .to_uppercase()
                    .chars()
                    .filter(|c| !c.is_whitespace())
                    .collect();
                records.push(FastaRecord { id, seq });
            }
            Err(e) => return Err(format!("Error reading FASTA record: {}", e)),
        }
    }
    Ok(records)
}

pub fn parse_fasta_string(text: &str) -> Result<Vec<Oligo>, String> {
    let cursor = std::io::Cursor::new(text);
    let reader = fasta::Reader::new(cursor);

    let mut oligos = Vec::new();
    for result in reader.records() {
        match result {
            Ok(record) => {
                let id = record.id().to_string();
                let seq = std::str::from_utf8(record.seq())
                    .map_err(|e| format!("Invalid UTF-8 in sequence: {}", e))?
                    .to_uppercase()
                    .chars()
                    .filter(|c| !c.is_whitespace())
                    .collect();
                oligos.push(Oligo { id, seq });
            }
            Err(e) => return Err(format!("Error parsing FASTA: {}", e)),
        }
    }
    Ok(oligos)
}

pub fn oligos_to_fasta_string(oligos: &[Oligo]) -> String {
    oligos
        .iter()
        .map(|o| format!(">{}\n{}", o.id, o.seq))
        .collect::<Vec<_>>()
        .join("\n")
}

pub fn get_best_alignment(
    target_seq: &str,
    oligo_seq: &str,
    settings: &AlignmentSettings,
) -> (
    Option<bio::alignment::Alignment>,
    Option<Orientation>,
    String,
) {
    let scoring = Scoring::new(
        settings.gap_open_score,
        settings.gap_extend_score,
        |a: u8, b: u8| {
            if iupac_match(a, b) {
                settings.match_score
            } else {
                settings.mismatch_score
            }
        },
    );

    let target_bytes = target_seq.as_bytes();
    let oligo_bytes = oligo_seq.as_bytes();
    let oligo_rc = reverse_complement(oligo_seq);
    let oligo_rc_bytes = oligo_rc.as_bytes();

    let mut aligner = Aligner::with_capacity_and_scoring(
        target_bytes.len(),
        oligo_bytes.len(),
        scoring.clone(),
    );

    let alignment_sense = match settings.mode {
        AlignmentMode::Local => aligner.local(target_bytes, oligo_bytes),
        AlignmentMode::Global => aligner.global(target_bytes, oligo_bytes),
    };

    let alignment_antisense = match settings.mode {
        AlignmentMode::Local => aligner.local(target_bytes, oligo_rc_bytes),
        AlignmentMode::Global => aligner.global(target_bytes, oligo_rc_bytes),
    };

    if alignment_sense.score >= alignment_antisense.score {
        (
            Some(alignment_sense),
            Some(Orientation::Sense),
            oligo_seq.to_string(),
        )
    } else {
        (
            Some(alignment_antisense),
            Some(Orientation::Antisense),
            oligo_rc,
        )
    }
}

pub fn is_valid_alignment(
    alignment: &bio::alignment::Alignment,
    oligo_len: usize,
    min_coverage: f64,
) -> bool {
    let aligned_length = alignment
        .operations
        .iter()
        .filter(|op| matches!(op, AlignmentOperation::Match | AlignmentOperation::Subst))
        .count();

    let coverage = aligned_length as f64 / oligo_len as f64;
    coverage >= min_coverage
}

pub fn get_alignment_coverage(alignment: &bio::alignment::Alignment, oligo_len: usize) -> f64 {
    let aligned_length = alignment
        .operations
        .iter()
        .filter(|op| matches!(op, AlignmentOperation::Match | AlignmentOperation::Subst))
        .count();

    aligned_length as f64 / oligo_len as f64
}

pub fn is_exact_match(a: u8, b: u8) -> bool {
    let a_upper = a.to_ascii_uppercase();
    let b_upper = b.to_ascii_uppercase();

    let normalized_a = if a_upper == b'U' { b'T' } else { a_upper };
    let normalized_b = if b_upper == b'U' { b'T' } else { b_upper };

    normalized_a == normalized_b
}

pub fn generate_signature(
    alignment: &bio::alignment::Alignment,
    target_seq: &str,
    oligo_seq: &str,
    orientation: Orientation,
    ambiguity_display: AmbiguityDisplayMode,
) -> (String, usize) {
    let oligo_len = oligo_seq.len();
    let mut sig: Vec<char> = vec!['-'; oligo_len];
    let mut mismatches = 0;

    let target_bytes = target_seq.as_bytes();
    let oligo_bytes = oligo_seq.as_bytes();

    let mut t_pos = alignment.xstart;
    let mut q_pos = alignment.ystart;

    for op in &alignment.operations {
        match op {
            AlignmentOperation::Match => {
                if q_pos < oligo_len && t_pos < target_bytes.len() {
                    let target_base = target_bytes[t_pos];
                    let oligo_base = oligo_bytes[q_pos];

                    if is_exact_match(target_base, oligo_base) {
                        sig[q_pos] = '.';
                    } else {
                        match ambiguity_display {
                            AmbiguityDisplayMode::ShowDots => {
                                sig[q_pos] = '.';
                            }
                            AmbiguityDisplayMode::ShowBases => {
                                sig[q_pos] = target_bytes[t_pos] as char;
                                mismatches += 1;
                            }
                        }
                    }
                }
                t_pos += 1;
                q_pos += 1;
            }
            AlignmentOperation::Subst => {
                if q_pos < oligo_len && t_pos < target_bytes.len() {
                    let target_base = target_bytes[t_pos];
                    let oligo_base = oligo_bytes[q_pos];

                    if iupac_match(target_base, oligo_base) {
                        match ambiguity_display {
                            AmbiguityDisplayMode::ShowDots => {
                                sig[q_pos] = '.';
                            }
                            AmbiguityDisplayMode::ShowBases => {
                                sig[q_pos] = target_base as char;
                                mismatches += 1;
                            }
                        }
                    } else {
                        sig[q_pos] = target_base as char;
                        mismatches += 1;
                    }
                }
                t_pos += 1;
                q_pos += 1;
            }
            AlignmentOperation::Del => {
                t_pos += 1;
            }
            AlignmentOperation::Ins => {
                if q_pos < oligo_len {
                    sig[q_pos] = '-';
                    mismatches += 1;
                }
                q_pos += 1;
            }
            AlignmentOperation::Xclip(_) => {}
            AlignmentOperation::Yclip(_) => {}
        }
    }

    mismatches += sig.iter().filter(|&&c| c == '-').count();

    let mut signature: String = sig.into_iter().collect();

    if orientation == Orientation::Antisense {
        signature = signature
            .chars()
            .map(|c| match c.to_ascii_uppercase() {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                'U' => 'A',
                'R' => 'Y',
                'Y' => 'R',
                'S' => 'S',
                'W' => 'W',
                'K' => 'M',
                'M' => 'K',
                'B' => 'V',
                'D' => 'H',
                'H' => 'D',
                'V' => 'B',
                _ => c,
            })
            .collect();
        signature = signature.chars().rev().collect();
    }

    (signature, mismatches)
}

pub fn align_oligos(
    sequence: &FastaRecord,
    oligos: &[Oligo],
    settings: &AlignmentSettings,
) -> HashMap<String, OligoResult> {
    let mut results = HashMap::new();

    for oligo in oligos {
        let (alignment_opt, orientation_opt, actual_oligo) =
            get_best_alignment(&sequence.seq, &oligo.seq, settings);

        let result =
            if let (Some(alignment), Some(orientation)) = (alignment_opt, orientation_opt) {
                if is_valid_alignment(&alignment, actual_oligo.len(), settings.min_coverage) {
                    let (signature, mismatches) = generate_signature(
                        &alignment,
                        &sequence.seq,
                        &actual_oligo,
                        orientation,
                        settings.ambiguity_display,
                    );

                    if mismatches > settings.max_mismatches_per_oligo {
                        OligoResult::no_match()
                    } else {
                        let coverage = get_alignment_coverage(&alignment, actual_oligo.len());

                        OligoResult {
                            matched: true,
                            orientation: Some(orientation),
                            signature,
                            mismatches,
                            score: alignment.score,
                            coverage,
                            start_pos: Some(alignment.xstart),
                            end_pos: Some(alignment.xend),
                        }
                    }
                } else {
                    OligoResult::no_match()
                }
            } else {
                OligoResult::no_match()
            };

        results.insert(oligo.id.clone(), result);
    }

    results
}

pub fn find_best_amplicon(
    fwd_results: &mut HashMap<String, OligoResult>,
    rev_results: &mut HashMap<String, OligoResult>,
    probe_results: &mut HashMap<String, OligoResult>,
    min_size: Option<usize>,
    max_size: Option<usize>,
) -> AmpliconInfo {
    let forward_hits: Vec<(String, usize, usize)> = fwd_results
        .iter()
        .filter(|(_, r)| r.matched)
        .filter_map(|(id, r)| match (r.start_pos, r.end_pos) {
            (Some(s), Some(e)) => Some((id.clone(), s, e)),
            _ => None,
        })
        .collect();

    let reverse_hits: Vec<(String, usize, usize)> = rev_results
        .iter()
        .filter(|(_, r)| r.matched)
        .filter_map(|(id, r)| match (r.start_pos, r.end_pos) {
            (Some(s), Some(e)) => Some((id.clone(), s, e)),
            _ => None,
        })
        .collect();

    let mut valid_amplicons: Vec<(String, String, usize, usize, usize)> = Vec::new();

    for (fwd_id, fwd_start, _fwd_end) in &forward_hits {
        for (rev_id, _rev_start, rev_end) in &reverse_hits {
            if fwd_start < rev_end {
                let amplicon_size = rev_end - fwd_start + 1;
                let size_ok = match (min_size, max_size) {
                    (Some(min), Some(max)) => amplicon_size >= min && amplicon_size <= max,
                    (Some(min), None) => amplicon_size >= min,
                    (None, Some(max)) => amplicon_size <= max,
                    (None, None) => true,
                };
                if size_ok {
                    valid_amplicons.push((
                        fwd_id.clone(),
                        rev_id.clone(),
                        *fwd_start,
                        *rev_end,
                        amplicon_size,
                    ));
                }
            }
        }
    }

    let mark_all_no_match = |results: &mut HashMap<String, OligoResult>| {
        for result in results.values_mut() {
            *result = OligoResult::no_match();
        }
    };

    if valid_amplicons.is_empty() {
        mark_all_no_match(fwd_results);
        mark_all_no_match(rev_results);
        mark_all_no_match(probe_results);
        return AmpliconInfo::default();
    }

    let best = valid_amplicons
        .iter()
        .max_by_key(|a| a.4)
        .unwrap()
        .clone();

    let (best_fwd_id, best_rev_id, amp_start, amp_end, amp_size) = best;

    let filter_by_bounds = |results: &mut HashMap<String, OligoResult>| {
        for result in results.values_mut() {
            if !result.matched {
                continue;
            }
            if let (Some(start), Some(end)) = (result.start_pos, result.end_pos) {
                if start < amp_start || end > amp_end {
                    *result = OligoResult::no_match();
                }
            } else {
                *result = OligoResult::no_match();
            }
        }
    };

    filter_by_bounds(fwd_results);
    filter_by_bounds(rev_results);
    filter_by_bounds(probe_results);

    AmpliconInfo {
        found: true,
        forward_oligo_id: Some(best_fwd_id),
        reverse_oligo_id: Some(best_rev_id),
        start: amp_start,
        end: amp_end,
        size: amp_size,
    }
}

pub fn analyze_sequence(
    sequence: &FastaRecord,
    fwd_primers: &[Oligo],
    rev_primers: &[Oligo],
    probes: &[Oligo],
    settings: &AlignmentSettings,
) -> SequenceResult {
    let mut fwd_results = align_oligos(sequence, fwd_primers, settings);
    let mut rev_results = align_oligos(sequence, rev_primers, settings);

    let mut probe_results: HashMap<String, OligoResult> = HashMap::new();
    let amplicon_info;

    let has_fwd_match = fwd_results.values().any(|r| r.matched);
    let has_rev_match = rev_results.values().any(|r| r.matched);

    if has_fwd_match && has_rev_match {
        probe_results = align_oligos(sequence, probes, settings);

        amplicon_info = Some(find_best_amplicon(
            &mut fwd_results,
            &mut rev_results,
            &mut probe_results,
            settings.min_amplicon_size,
            settings.max_amplicon_size,
        ));
    } else {
        for result in fwd_results.values_mut() {
            *result = OligoResult::no_match();
        }
        for result in rev_results.values_mut() {
            *result = OligoResult::no_match();
        }
        for probe in probes {
            probe_results.insert(probe.id.clone(), OligoResult::no_match());
        }
        amplicon_info = Some(AmpliconInfo::default());
    }

    let fwd_matched = fwd_results.values().filter(|r| r.matched).count();
    let rev_matched = rev_results.values().filter(|r| r.matched).count();
    let probe_matched = probe_results.values().filter(|r| r.matched).count();

    let best_fwd_mm = fwd_results
        .values()
        .filter(|r| r.matched)
        .map(|r| r.mismatches)
        .min();
    let best_rev_mm = rev_results
        .values()
        .filter(|r| r.matched)
        .map(|r| r.mismatches)
        .min();
    let best_probe_mm = probe_results
        .values()
        .filter(|r| r.matched)
        .map(|r| r.mismatches)
        .min();

    let total_mismatches: usize = best_fwd_mm.unwrap_or(0)
        + best_rev_mm.unwrap_or(0)
        + best_probe_mm.unwrap_or(0);

    SequenceResult {
        fwd_results,
        rev_results,
        probe_results,
        fwd_matched,
        rev_matched,
        probe_matched,
        total_mismatches,
        amplicon_info,
        best_fwd_mm,
        best_rev_mm,
        best_probe_mm,
    }
}

pub fn build_signature_parts(
    oligos: &[Oligo],
    results: &HashMap<String, OligoResult>,
) -> Vec<String> {
    oligos
        .iter()
        .map(|oligo| {
            if let Some(result) = results.get(&oligo.id) {
                if result.matched {
                    let orientation_symbol = match result.orientation {
                        Some(Orientation::Sense) => "(fwd)",
                        Some(Orientation::Antisense) => "(rev)",
                        None => "",
                    };
                    format!("{}{}", result.signature, orientation_symbol)
                } else {
                    "NO_MATCH".to_string()
                }
            } else {
                "NO_MATCH".to_string()
            }
        })
        .collect()
}
