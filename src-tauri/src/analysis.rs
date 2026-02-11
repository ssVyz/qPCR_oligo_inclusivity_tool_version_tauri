use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use tauri::Emitter;

use crate::bio_logic::*;
use crate::models::*;

/// Run the full analysis with parallel processing, emitting progress events via Tauri
pub fn run_analysis(
    sequences: &[FastaRecord],
    fwd_primers: &[Oligo],
    rev_primers: &[Oligo],
    probes: &[Oligo],
    settings: &AlignmentSettings,
    running: &Arc<std::sync::atomic::AtomicBool>,
    app_handle: &tauri::AppHandle,
) -> AnalysisResults {
    let total_sequences = sequences.len();

    let mut initial_stats: HashMap<String, OligoStats> = HashMap::new();
    for oligo in fwd_primers
        .iter()
        .chain(rev_primers.iter())
        .chain(probes.iter())
    {
        initial_stats.insert(oligo.id.clone(), OligoStats::default());
    }
    let oligo_stats = Arc::new(Mutex::new(initial_stats));

    let alignment_dict: Arc<Mutex<HashMap<String, PatternData>>> =
        Arc::new(Mutex::new(HashMap::new()));
    let sequences_with_min_matches = Arc::new(AtomicUsize::new(0));
    let sequences_with_valid_amplicon = Arc::new(AtomicUsize::new(0));
    let sequences_failed_amplicon = Arc::new(AtomicUsize::new(0));

    let fwd_mm_zero = Arc::new(AtomicUsize::new(0));
    let fwd_mm_one = Arc::new(AtomicUsize::new(0));
    let fwd_mm_more = Arc::new(AtomicUsize::new(0));
    let fwd_mm_none = Arc::new(AtomicUsize::new(0));
    let rev_mm_zero = Arc::new(AtomicUsize::new(0));
    let rev_mm_one = Arc::new(AtomicUsize::new(0));
    let rev_mm_more = Arc::new(AtomicUsize::new(0));
    let rev_mm_none = Arc::new(AtomicUsize::new(0));
    let probe_mm_zero = Arc::new(AtomicUsize::new(0));
    let probe_mm_one = Arc::new(AtomicUsize::new(0));
    let probe_mm_more = Arc::new(AtomicUsize::new(0));
    let probe_mm_none = Arc::new(AtomicUsize::new(0));

    let overall_all_perfect = Arc::new(AtomicUsize::new(0));
    let overall_max_one_mm = Arc::new(AtomicUsize::new(0));
    let overall_two_plus_mm = Arc::new(AtomicUsize::new(0));
    let overall_no_match = Arc::new(AtomicUsize::new(0));

    let processed = Arc::new(AtomicUsize::new(0));
    let app_handle_clone = app_handle.clone();

    sequences.par_iter().for_each(|record| {
        if !running.load(Ordering::SeqCst) {
            return;
        }

        let seq_result = analyze_sequence(record, fwd_primers, rev_primers, probes, settings);

        if let Some(ref info) = seq_result.amplicon_info {
            if info.found {
                sequences_with_valid_amplicon.fetch_add(1, Ordering::SeqCst);
            } else {
                sequences_failed_amplicon.fetch_add(1, Ordering::SeqCst);
            }
        }

        match seq_result.best_fwd_mm {
            Some(0) => {
                fwd_mm_zero.fetch_add(1, Ordering::SeqCst);
            }
            Some(1) => {
                fwd_mm_one.fetch_add(1, Ordering::SeqCst);
            }
            Some(_) => {
                fwd_mm_more.fetch_add(1, Ordering::SeqCst);
            }
            None => {
                fwd_mm_none.fetch_add(1, Ordering::SeqCst);
            }
        }
        match seq_result.best_rev_mm {
            Some(0) => {
                rev_mm_zero.fetch_add(1, Ordering::SeqCst);
            }
            Some(1) => {
                rev_mm_one.fetch_add(1, Ordering::SeqCst);
            }
            Some(_) => {
                rev_mm_more.fetch_add(1, Ordering::SeqCst);
            }
            None => {
                rev_mm_none.fetch_add(1, Ordering::SeqCst);
            }
        }
        if !probes.is_empty() {
            match seq_result.best_probe_mm {
                Some(0) => {
                    probe_mm_zero.fetch_add(1, Ordering::SeqCst);
                }
                Some(1) => {
                    probe_mm_one.fetch_add(1, Ordering::SeqCst);
                }
                Some(_) => {
                    probe_mm_more.fetch_add(1, Ordering::SeqCst);
                }
                None => {
                    probe_mm_none.fetch_add(1, Ordering::SeqCst);
                }
            }
        }

        {
            let mut category_bests: Vec<Option<usize>> =
                vec![seq_result.best_fwd_mm, seq_result.best_rev_mm];
            if !probes.is_empty() {
                category_bests.push(seq_result.best_probe_mm);
            }

            if category_bests.iter().any(|b| b.is_none()) {
                overall_no_match.fetch_add(1, Ordering::SeqCst);
            } else {
                let worst_best = category_bests
                    .iter()
                    .filter_map(|b| *b)
                    .max()
                    .unwrap_or(0);
                match worst_best {
                    0 => {
                        overall_all_perfect.fetch_add(1, Ordering::SeqCst);
                    }
                    1 => {
                        overall_max_one_mm.fetch_add(1, Ordering::SeqCst);
                    }
                    _ => {
                        overall_two_plus_mm.fetch_add(1, Ordering::SeqCst);
                    }
                }
            }
        }

        {
            let mut stats = oligo_stats.lock().unwrap();
            let all_results = seq_result
                .fwd_results
                .iter()
                .chain(seq_result.rev_results.iter())
                .chain(seq_result.probe_results.iter());
            for (oligo_id, result) in all_results {
                if result.matched {
                    if let Some(s) = stats.get_mut(oligo_id) {
                        s.total_matches += 1;
                        match result.orientation {
                            Some(Orientation::Sense) => s.sense_matches += 1,
                            Some(Orientation::Antisense) => s.antisense_matches += 1,
                            None => {}
                        }
                    }
                }
            }
        }

        let meets_fwd = seq_result.fwd_matched >= settings.min_fwd_matched;
        let meets_rev = seq_result.rev_matched >= settings.min_rev_matched;
        let meets_probe = seq_result.probe_matched >= settings.min_probe_matched;

        if meets_fwd && meets_rev && meets_probe {
            sequences_with_min_matches.fetch_add(1, Ordering::SeqCst);

            let fwd_parts = build_signature_parts(fwd_primers, &seq_result.fwd_results);
            let probe_parts = build_signature_parts(probes, &seq_result.probe_results);
            let rev_parts = build_signature_parts(rev_primers, &seq_result.rev_results);

            let mut all_parts = Vec::new();
            all_parts.push(fwd_parts.join(" | "));
            if !probes.is_empty() {
                all_parts.push(probe_parts.join(" | "));
            }
            all_parts.push(rev_parts.join(" | "));

            let combined_signature = all_parts.join(" || ");

            {
                let mut dict = alignment_dict.lock().unwrap();
                let entry = dict
                    .entry(combined_signature)
                    .or_insert_with(|| PatternData {
                        count: 0,
                        total_mismatches: seq_result.total_mismatches,
                        matched_fwd: seq_result.fwd_matched,
                        matched_rev: seq_result.rev_matched,
                        matched_probe: seq_result.probe_matched,
                        examples: Vec::new(),
                        amplicon_lengths: Vec::new(),
                    });
                entry.count += 1;
                if entry.examples.len() < 10 {
                    entry.examples.push(record.id.clone());
                }
                if let Some(ref info) = seq_result.amplicon_info {
                    if info.found && info.size > 0 {
                        entry.amplicon_lengths.push(info.size);
                    }
                }
            }
        }

        let current = processed.fetch_add(1, Ordering::SeqCst) + 1;

        if current % 100 == 0 || current == total_sequences {
            let _ = app_handle_clone.emit(
                "analysis-progress",
                ProgressPayload {
                    current,
                    total: total_sequences,
                    status: format!("Processing sequence {}/{}", current, total_sequences),
                },
            );
        }
    });

    let final_oligo_stats = oligo_stats.lock().unwrap().clone();
    let final_alignment_dict = alignment_dict.lock().unwrap().clone();
    let final_sequences_with_min_matches = sequences_with_min_matches.load(Ordering::SeqCst);
    let final_sequences_with_valid_amplicon =
        sequences_with_valid_amplicon.load(Ordering::SeqCst);
    let final_sequences_failed_amplicon = sequences_failed_amplicon.load(Ordering::SeqCst);

    let final_fwd_mm_dist = MismatchDistribution {
        zero_mm: fwd_mm_zero.load(Ordering::SeqCst),
        one_mm: fwd_mm_one.load(Ordering::SeqCst),
        more_mm: fwd_mm_more.load(Ordering::SeqCst),
        no_match: fwd_mm_none.load(Ordering::SeqCst),
    };
    let final_rev_mm_dist = MismatchDistribution {
        zero_mm: rev_mm_zero.load(Ordering::SeqCst),
        one_mm: rev_mm_one.load(Ordering::SeqCst),
        more_mm: rev_mm_more.load(Ordering::SeqCst),
        no_match: rev_mm_none.load(Ordering::SeqCst),
    };
    let final_probe_mm_dist = MismatchDistribution {
        zero_mm: probe_mm_zero.load(Ordering::SeqCst),
        one_mm: probe_mm_one.load(Ordering::SeqCst),
        more_mm: probe_mm_more.load(Ordering::SeqCst),
        no_match: probe_mm_none.load(Ordering::SeqCst),
    };

    let final_overall_all_perfect = overall_all_perfect.load(Ordering::SeqCst);
    let final_overall_max_one_mm = overall_max_one_mm.load(Ordering::SeqCst);
    let final_overall_two_plus_mm = overall_two_plus_mm.load(Ordering::SeqCst);
    let final_overall_no_match = overall_no_match.load(Ordering::SeqCst);

    let output_text = generate_output_text(
        fwd_primers,
        rev_primers,
        probes,
        &final_oligo_stats,
        &final_alignment_dict,
        total_sequences,
        final_sequences_with_min_matches,
        final_sequences_with_valid_amplicon,
        final_sequences_failed_amplicon,
        &final_fwd_mm_dist,
        &final_rev_mm_dist,
        &final_probe_mm_dist,
        final_overall_all_perfect,
        final_overall_max_one_mm,
        final_overall_two_plus_mm,
        final_overall_no_match,
        settings,
    );

    AnalysisResults {
        alignment_dict: final_alignment_dict,
        oligo_stats: final_oligo_stats,
        total_sequences,
        sequences_with_min_matches: final_sequences_with_min_matches,
        sequences_with_valid_amplicon: final_sequences_with_valid_amplicon,
        sequences_failed_amplicon: final_sequences_failed_amplicon,
        fwd_mm_dist: final_fwd_mm_dist,
        rev_mm_dist: final_rev_mm_dist,
        probe_mm_dist: final_probe_mm_dist,
        overall_all_perfect: final_overall_all_perfect,
        overall_max_one_mm: final_overall_max_one_mm,
        overall_two_plus_mm: final_overall_two_plus_mm,
        overall_no_match: final_overall_no_match,
        output_text,
    }
}

pub fn generate_output_text(
    fwd_primers: &[Oligo],
    rev_primers: &[Oligo],
    probes: &[Oligo],
    oligo_stats: &HashMap<String, OligoStats>,
    alignment_dict: &HashMap<String, PatternData>,
    total_sequences: usize,
    sequences_with_min_matches: usize,
    sequences_with_valid_amplicon: usize,
    sequences_failed_amplicon: usize,
    fwd_mm_dist: &MismatchDistribution,
    rev_mm_dist: &MismatchDistribution,
    probe_mm_dist: &MismatchDistribution,
    overall_all_perfect: usize,
    overall_max_one_mm: usize,
    overall_two_plus_mm: usize,
    overall_no_match: usize,
    settings: &AlignmentSettings,
) -> String {
    let mut out = Vec::new();

    out.push("=".repeat(80));
    out.push("qPCR OLIGO INCLUSIVITY ANALYSIS RESULTS".to_string());
    out.push("=".repeat(80));
    out.push(String::new());

    out.push("FORWARD PRIMERS:".to_string());
    for o in fwd_primers {
        out.push(format!("  {} ({})", o.id, o.seq));
    }
    out.push(String::new());

    if !probes.is_empty() {
        out.push("PROBES:".to_string());
        for o in probes {
            out.push(format!("  {} ({})", o.id, o.seq));
        }
        out.push(String::new());
    }

    out.push("REVERSE PRIMERS:".to_string());
    for o in rev_primers {
        out.push(format!("  {} ({})", o.id, o.seq));
    }
    out.push(String::new());

    out.push(format!(
        "Analysis settings: Min fwd matched = {}, Min rev matched = {}, Min probes matched = {}, Min coverage = {}, Max mismatches/oligo = {}",
        settings.min_fwd_matched, settings.min_rev_matched, settings.min_probe_matched, settings.min_coverage, settings.max_mismatches_per_oligo
    ));
    if settings.min_amplicon_size.is_some() || settings.max_amplicon_size.is_some() {
        match (settings.min_amplicon_size, settings.max_amplicon_size) {
            (Some(min), Some(max)) => {
                out.push(format!("Amplicon size constraint: {} - {} bp", min, max));
            }
            (Some(min), None) => {
                out.push(format!("Amplicon size constraint: >= {} bp", min));
            }
            (None, Some(max)) => {
                out.push(format!("Amplicon size constraint: <= {} bp", max));
            }
            (None, None) => {}
        }
    }
    out.push(String::new());

    out.push("SEQUENCE SIGNATURE PATTERNS:".to_string());
    out.push("Column order: [Forward Primers] || [Probes] || [Reverse Primers]".to_string());
    out.push("-".repeat(50));

    if !alignment_dict.is_empty() {
        let mut sorted_patterns: Vec<_> = alignment_dict.iter().collect();
        sorted_patterns.sort_by(|a, b| b.1.count.cmp(&a.1.count));

        for (signature, data) in sorted_patterns {
            let mut examples_str = data.examples.iter().take(3).cloned().collect::<Vec<_>>();
            if data.examples.len() > 3 {
                examples_str.push(format!("... (+{} more)", data.examples.len() - 3));
            }

            out.push(format!("Pattern: {}", signature));
            out.push(format!(
                "  Count: {}, Mismatches: {}, Fwd matched: {}, Rev matched: {}, Probes matched: {}",
                data.count, data.total_mismatches, data.matched_fwd, data.matched_rev, data.matched_probe
            ));
            if !data.amplicon_lengths.is_empty() {
                let avg_len: f64 = data.amplicon_lengths.iter().sum::<usize>() as f64
                    / data.amplicon_lengths.len() as f64;
                out.push(format!("  Amplicon length: ~{:.0} bp", avg_len));
            }
            out.push(format!("  Examples: {}", examples_str.join(", ")));
            out.push(String::new());
        }
    } else {
        out.push("No sequences met the minimum matching criteria.".to_string());
        out.push(String::new());
    }

    out.push("PER-OLIGO STATISTICS:".to_string());
    out.push("-".repeat(30));

    out.push("  Forward Primers:".to_string());
    for oligo in fwd_primers {
        if let Some(stats) = oligo_stats.get(&oligo.id) {
            let percentage = if total_sequences > 0 {
                (stats.total_matches as f64 / total_sequences as f64) * 100.0
            } else {
                0.0
            };
            out.push(format!(
                "    {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                oligo.id,
                stats.total_matches,
                total_sequences,
                percentage,
                stats.sense_matches,
                stats.antisense_matches
            ));
        }
    }

    if !probes.is_empty() {
        out.push("  Probes:".to_string());
        for oligo in probes {
            if let Some(stats) = oligo_stats.get(&oligo.id) {
                let percentage = if total_sequences > 0 {
                    (stats.total_matches as f64 / total_sequences as f64) * 100.0
                } else {
                    0.0
                };
                out.push(format!(
                    "    {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                    oligo.id,
                    stats.total_matches,
                    total_sequences,
                    percentage,
                    stats.sense_matches,
                    stats.antisense_matches
                ));
            }
        }
    }

    out.push("  Reverse Primers:".to_string());
    for oligo in rev_primers {
        if let Some(stats) = oligo_stats.get(&oligo.id) {
            let percentage = if total_sequences > 0 {
                (stats.total_matches as f64 / total_sequences as f64) * 100.0
            } else {
                0.0
            };
            out.push(format!(
                "    {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                oligo.id,
                stats.total_matches,
                total_sequences,
                percentage,
                stats.sense_matches,
                stats.antisense_matches
            ));
        }
    }

    out.push(String::new());
    out.push("SUMMARY:".to_string());
    out.push("-".repeat(20));
    out.push(format!("Total sequences analyzed: {}", total_sequences));
    let percentage = if total_sequences > 0 {
        (sequences_with_min_matches as f64 / total_sequences as f64) * 100.0
    } else {
        0.0
    };
    out.push(format!(
        "Sequences meeting all thresholds: {} ({:.1}%)",
        sequences_with_min_matches, percentage
    ));

    out.push(String::new());
    out.push("AMPLICON STATISTICS:".to_string());
    out.push("-".repeat(20));
    let amp_percentage = if total_sequences > 0 {
        (sequences_with_valid_amplicon as f64 / total_sequences as f64) * 100.0
    } else {
        0.0
    };
    out.push(format!(
        "Sequences with valid amplicon (fwd+rev pair): {} ({:.1}%)",
        sequences_with_valid_amplicon, amp_percentage
    ));
    out.push(format!(
        "Sequences without valid amplicon: {}",
        sequences_failed_amplicon
    ));
    if settings.min_amplicon_size.is_some() || settings.max_amplicon_size.is_some() {
        let constraint_desc = match (settings.min_amplicon_size, settings.max_amplicon_size) {
            (Some(min), Some(max)) => format!("between {} and {} bp", min, max),
            (Some(min), None) => format!(">= {} bp", min),
            (None, Some(max)) => format!("<= {} bp", max),
            (None, None) => "".to_string(),
        };
        out.push(format!("(Amplicon size constraint: {})", constraint_desc));
    }

    out.push(String::new());
    out.push("MISMATCH DISTRIBUTION (best oligo per category per sequence):".to_string());
    out.push("-".repeat(50));

    let fmt_pct = |count: usize| -> String {
        if total_sequences > 0 {
            format!(
                "{} ({:.1}%)",
                count,
                count as f64 / total_sequences as f64 * 100.0
            )
        } else {
            format!("{}", count)
        }
    };

    out.push("  Forward Primers:".to_string());
    out.push(format!(
        "    0 mismatches: {}",
        fmt_pct(fwd_mm_dist.zero_mm)
    ));
    out.push(format!(
        "    1 mismatch:   {}",
        fmt_pct(fwd_mm_dist.one_mm)
    ));
    out.push(format!(
        "    >1 mismatches: {}",
        fmt_pct(fwd_mm_dist.more_mm)
    ));
    out.push(format!(
        "    No match:     {}",
        fmt_pct(fwd_mm_dist.no_match)
    ));

    if !probes.is_empty() {
        out.push("  Probes:".to_string());
        out.push(format!(
            "    0 mismatches: {}",
            fmt_pct(probe_mm_dist.zero_mm)
        ));
        out.push(format!(
            "    1 mismatch:   {}",
            fmt_pct(probe_mm_dist.one_mm)
        ));
        out.push(format!(
            "    >1 mismatches: {}",
            fmt_pct(probe_mm_dist.more_mm)
        ));
        out.push(format!(
            "    No match:     {}",
            fmt_pct(probe_mm_dist.no_match)
        ));
    }

    out.push("  Reverse Primers:".to_string());
    out.push(format!(
        "    0 mismatches: {}",
        fmt_pct(rev_mm_dist.zero_mm)
    ));
    out.push(format!(
        "    1 mismatch:   {}",
        fmt_pct(rev_mm_dist.one_mm)
    ));
    out.push(format!(
        "    >1 mismatches: {}",
        fmt_pct(rev_mm_dist.more_mm)
    ));
    out.push(format!(
        "    No match:     {}",
        fmt_pct(rev_mm_dist.no_match)
    ));

    out.push(String::new());
    out.push("  Overall Pattern (worst best-match across all categories):".to_string());
    out.push(format!(
        "    All categories 0 mismatches:  {}",
        fmt_pct(overall_all_perfect)
    ));
    out.push(format!(
        "    All categories ≤1 mismatch:   {}",
        fmt_pct(overall_max_one_mm)
    ));
    out.push(format!(
        "    ≥2 mismatches in any category: {}",
        fmt_pct(overall_two_plus_mm)
    ));
    out.push(format!(
        "    No match in any category:     {}",
        fmt_pct(overall_no_match)
    ));

    out.push(String::new());
    out.push("LEGEND:".to_string());
    out.push("'(fwd)' = sense orientation, '(rev)' = antisense orientation".to_string());
    out.push("'.' = match, letter = mismatch, '-' = gap/unaligned".to_string());
    out.push(
        "'||' separates oligo categories: Forward Primers || Probes || Reverse Primers".to_string(),
    );
    out.push("'|' separates individual oligos within a category".to_string());
    out.push(
        "Probes are only matched within the amplicon region defined by the best fwd+rev pair"
            .to_string(),
    );
    out.push("=".repeat(80));

    out.join("\n")
}
