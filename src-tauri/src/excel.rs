use rust_xlsxwriter::{Format, Workbook};
use std::path::PathBuf;

use crate::models::*;

pub fn write_excel(
    path: &PathBuf,
    results: &AnalysisResults,
    fwd_primers: &[Oligo],
    rev_primers: &[Oligo],
    probes: &[Oligo],
    settings: &AlignmentSettings,
) -> Result<(), String> {
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();

    let header_format = Format::new().set_bold();
    let title_format = Format::new().set_bold().set_font_size(14);
    let category_format = Format::new().set_bold().set_font_size(11);

    let mut row: u32 = 0;

    worksheet
        .write_string_with_format(
            row,
            0,
            "qPCR Oligo Inclusivity Analysis Results",
            &title_format,
        )
        .map_err(|e| e.to_string())?;
    row += 2;

    worksheet
        .write_string(
            row,
            0,
            &format!(
                "Settings: Min fwd = {}, Min rev = {}, Min probes = {}, Min coverage = {}, Max mm/oligo = {}",
                settings.min_fwd_matched,
                settings.min_rev_matched,
                settings.min_probe_matched,
                settings.min_coverage,
                settings.max_mismatches_per_oligo
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;

    if settings.min_amplicon_size.is_some() || settings.max_amplicon_size.is_some() {
        let constraint_text = match (settings.min_amplicon_size, settings.max_amplicon_size) {
            (Some(min), Some(max)) => {
                format!("Amplicon size constraint: {} - {} bp", min, max)
            }
            (Some(min), None) => format!("Amplicon size constraint: >= {} bp", min),
            (None, Some(max)) => format!("Amplicon size constraint: <= {} bp", max),
            (None, None) => "".to_string(),
        };
        worksheet
            .write_string(row, 0, &constraint_text)
            .map_err(|e| e.to_string())?;
        row += 1;
    }
    row += 2;

    // Category labels row
    let mut col: u16 = 1;
    if !fwd_primers.is_empty() {
        worksheet
            .write_string_with_format(row, col, "--- Forward Primers ---", &category_format)
            .map_err(|e| e.to_string())?;
    }
    col += fwd_primers.len() as u16;
    if !probes.is_empty() {
        worksheet
            .write_string_with_format(row, col, "--- Probes ---", &category_format)
            .map_err(|e| e.to_string())?;
    }
    col += probes.len() as u16;
    if !rev_primers.is_empty() {
        worksheet
            .write_string_with_format(row, col, "--- Reverse Primers ---", &category_format)
            .map_err(|e| e.to_string())?;
    }
    row += 1;

    // Headers row
    col = 0;
    worksheet
        .write_string_with_format(row, col, "Pattern #", &header_format)
        .map_err(|e| e.to_string())?;
    col += 1;

    let ordered_oligos: Vec<&Oligo> = fwd_primers
        .iter()
        .chain(probes.iter())
        .chain(rev_primers.iter())
        .collect();

    for oligo in &ordered_oligos {
        worksheet
            .write_string_with_format(
                row,
                col,
                &format!("{} Pattern", oligo.id),
                &header_format,
            )
            .map_err(|e| e.to_string())?;
        col += 1;
    }

    worksheet
        .write_string_with_format(row, col, "Count", &header_format)
        .map_err(|e| e.to_string())?;
    col += 1;
    worksheet
        .write_string_with_format(row, col, "Percentage", &header_format)
        .map_err(|e| e.to_string())?;
    col += 1;
    worksheet
        .write_string_with_format(row, col, "Total Mismatches", &header_format)
        .map_err(|e| e.to_string())?;
    col += 1;
    worksheet
        .write_string_with_format(row, col, "Fwd Matched", &header_format)
        .map_err(|e| e.to_string())?;
    col += 1;
    worksheet
        .write_string_with_format(row, col, "Rev Matched", &header_format)
        .map_err(|e| e.to_string())?;
    col += 1;
    worksheet
        .write_string_with_format(row, col, "Probes Matched", &header_format)
        .map_err(|e| e.to_string())?;
    col += 1;
    worksheet
        .write_string_with_format(row, col, "Amplicon Length", &header_format)
        .map_err(|e| e.to_string())?;
    col += 1;
    worksheet
        .write_string_with_format(row, col, "Example Sequences", &header_format)
        .map_err(|e| e.to_string())?;
    row += 1;

    // Row with full oligo sequences
    col = 1;
    for oligo in &ordered_oligos {
        worksheet
            .write_string(row, col, &oligo.seq)
            .map_err(|e| e.to_string())?;
        col += 1;
    }
    row += 1;

    // Data rows
    let mut sorted_patterns: Vec<_> = results.alignment_dict.iter().collect();
    sorted_patterns.sort_by(|a, b| b.1.count.cmp(&a.1.count));

    let num_fwd = fwd_primers.len();
    let num_probes = probes.len();
    let num_rev = rev_primers.len();

    let mut pattern_num = 1u32;
    for (signature, data) in sorted_patterns {
        col = 0;
        worksheet
            .write_number(row, col, pattern_num as f64)
            .map_err(|e| e.to_string())?;
        col += 1;

        let category_sections: Vec<&str> = signature.split(" || ").collect();

        let mut all_patterns: Vec<String> = Vec::new();

        let fwd_section = category_sections.first().unwrap_or(&"");
        let fwd_patterns: Vec<&str> = if !fwd_section.is_empty() {
            fwd_section.split(" | ").collect()
        } else {
            Vec::new()
        };
        for i in 0..num_fwd {
            all_patterns.push(
                fwd_patterns
                    .get(i)
                    .unwrap_or(&"NO_MATCH")
                    .to_string(),
            );
        }

        if num_probes > 0 {
            let probe_section = if category_sections.len() >= 3 {
                category_sections.get(1).unwrap_or(&"")
            } else {
                &""
            };
            let probe_patterns: Vec<&str> = if !probe_section.is_empty() {
                probe_section.split(" | ").collect()
            } else {
                Vec::new()
            };
            for i in 0..num_probes {
                all_patterns.push(
                    probe_patterns
                        .get(i)
                        .unwrap_or(&"NO_MATCH")
                        .to_string(),
                );
            }
        }

        let rev_section_idx = if num_probes > 0 { 2 } else { 1 };
        let rev_section = category_sections.get(rev_section_idx).unwrap_or(&"");
        let rev_patterns: Vec<&str> = if !rev_section.is_empty() {
            rev_section.split(" | ").collect()
        } else {
            Vec::new()
        };
        for i in 0..num_rev {
            all_patterns.push(
                rev_patterns
                    .get(i)
                    .unwrap_or(&"NO_MATCH")
                    .to_string(),
            );
        }

        for pattern in &all_patterns {
            worksheet
                .write_string(row, col, pattern)
                .map_err(|e| e.to_string())?;
            col += 1;
        }

        worksheet
            .write_number(row, col, data.count as f64)
            .map_err(|e| e.to_string())?;
        col += 1;

        let percentage = if results.total_sequences > 0 {
            (data.count as f64 / results.total_sequences as f64) * 100.0
        } else {
            0.0
        };
        worksheet
            .write_number(row, col, percentage)
            .map_err(|e| e.to_string())?;
        col += 1;

        worksheet
            .write_number(row, col, data.total_mismatches as f64)
            .map_err(|e| e.to_string())?;
        col += 1;

        worksheet
            .write_number(row, col, data.matched_fwd as f64)
            .map_err(|e| e.to_string())?;
        col += 1;

        worksheet
            .write_number(row, col, data.matched_rev as f64)
            .map_err(|e| e.to_string())?;
        col += 1;

        worksheet
            .write_number(row, col, data.matched_probe as f64)
            .map_err(|e| e.to_string())?;
        col += 1;

        if !data.amplicon_lengths.is_empty() {
            use std::collections::HashMap;
            let mut length_counts: HashMap<usize, usize> = HashMap::new();
            for &len in &data.amplicon_lengths {
                *length_counts.entry(len).or_insert(0) += 1;
            }
            let amplicon_length = length_counts
                .iter()
                .max_by_key(|(_, &count)| count)
                .map(|(&len, _)| len)
                .unwrap_or(0);
            worksheet
                .write_number(row, col, amplicon_length as f64)
                .map_err(|e| e.to_string())?;
        } else {
            worksheet
                .write_string(row, col, "")
                .map_err(|e| e.to_string())?;
        }
        col += 1;

        let examples: Vec<_> = data.member_ids.iter().take(3).collect();
        let mut examples_str = examples
            .iter()
            .map(|s| s.as_str())
            .collect::<Vec<_>>()
            .join(", ");
        if data.member_ids.len() > 3 {
            examples_str.push_str(&format!(" (+{} more)", data.member_ids.len() - 3));
        }
        worksheet
            .write_string(row, col, &examples_str)
            .map_err(|e| e.to_string())?;

        row += 1;
        pattern_num += 1;
    }

    // Per-oligo Statistics
    row += 2;
    worksheet
        .write_string_with_format(row, 0, "PER-OLIGO STATISTICS:", &header_format)
        .map_err(|e| e.to_string())?;
    row += 1;

    worksheet
        .write_string_with_format(row, 0, "Forward Primers:", &category_format)
        .map_err(|e| e.to_string())?;
    row += 1;
    for oligo in fwd_primers {
        if let Some(stats) = results.oligo_stats.get(&oligo.id) {
            let percentage = if results.total_sequences > 0 {
                (stats.total_matches as f64 / results.total_sequences as f64) * 100.0
            } else {
                0.0
            };
            worksheet
                .write_string(
                    row,
                    0,
                    &format!(
                        "  {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                        oligo.id,
                        stats.total_matches,
                        results.total_sequences,
                        percentage,
                        stats.sense_matches,
                        stats.antisense_matches
                    ),
                )
                .map_err(|e| e.to_string())?;
            row += 1;
        }
    }

    if !probes.is_empty() {
        worksheet
            .write_string_with_format(row, 0, "Probes:", &category_format)
            .map_err(|e| e.to_string())?;
        row += 1;
        for oligo in probes {
            if let Some(stats) = results.oligo_stats.get(&oligo.id) {
                let percentage = if results.total_sequences > 0 {
                    (stats.total_matches as f64 / results.total_sequences as f64) * 100.0
                } else {
                    0.0
                };
                worksheet
                    .write_string(
                        row,
                        0,
                        &format!(
                            "  {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                            oligo.id,
                            stats.total_matches,
                            results.total_sequences,
                            percentage,
                            stats.sense_matches,
                            stats.antisense_matches
                        ),
                    )
                    .map_err(|e| e.to_string())?;
                row += 1;
            }
        }
    }

    worksheet
        .write_string_with_format(row, 0, "Reverse Primers:", &category_format)
        .map_err(|e| e.to_string())?;
    row += 1;
    for oligo in rev_primers {
        if let Some(stats) = results.oligo_stats.get(&oligo.id) {
            let percentage = if results.total_sequences > 0 {
                (stats.total_matches as f64 / results.total_sequences as f64) * 100.0
            } else {
                0.0
            };
            worksheet
                .write_string(
                    row,
                    0,
                    &format!(
                        "  {}: {}/{} matches ({:.1}%) - Sense: {}, Antisense: {}",
                        oligo.id,
                        stats.total_matches,
                        results.total_sequences,
                        percentage,
                        stats.sense_matches,
                        stats.antisense_matches
                    ),
                )
                .map_err(|e| e.to_string())?;
            row += 1;
        }
    }

    // Summary
    row += 1;
    worksheet
        .write_string_with_format(row, 0, "SUMMARY:", &header_format)
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!("Total sequences analyzed: {}", results.total_sequences),
        )
        .map_err(|e| e.to_string())?;
    row += 1;

    let percentage = if results.total_sequences > 0 {
        (results.sequences_with_min_matches as f64 / results.total_sequences as f64) * 100.0
    } else {
        0.0
    };
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "Sequences meeting all thresholds: {} ({:.1}%)",
                results.sequences_with_min_matches, percentage
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;

    // Amplicon statistics
    row += 1;
    worksheet
        .write_string_with_format(row, 0, "AMPLICON STATISTICS:", &header_format)
        .map_err(|e| e.to_string())?;
    row += 1;

    let amp_percentage = if results.total_sequences > 0 {
        (results.sequences_with_valid_amplicon as f64 / results.total_sequences as f64) * 100.0
    } else {
        0.0
    };
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "Sequences with valid amplicon: {} ({:.1}%)",
                results.sequences_with_valid_amplicon, amp_percentage
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;

    worksheet
        .write_string(
            row,
            0,
            &format!(
                "Sequences without valid amplicon: {}",
                results.sequences_failed_amplicon
            ),
        )
        .map_err(|e| e.to_string())?;

    // Mismatch distribution
    row += 2;
    worksheet
        .write_string_with_format(
            row,
            0,
            "MISMATCH DISTRIBUTION (best oligo per category per sequence):",
            &header_format,
        )
        .map_err(|e| e.to_string())?;
    row += 1;

    let fmt_pct_xl = |count: usize| -> String {
        if results.total_sequences > 0 {
            format!(
                "{} ({:.1}%)",
                count,
                count as f64 / results.total_sequences as f64 * 100.0
            )
        } else {
            format!("{}", count)
        }
    };

    worksheet
        .write_string_with_format(row, 0, "Forward Primers:", &category_format)
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "  0 mismatches: {}",
                fmt_pct_xl(results.fwd_mm_dist.zero_mm)
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!("  1 mismatch: {}", fmt_pct_xl(results.fwd_mm_dist.one_mm)),
        )
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "  >1 mismatches: {}",
                fmt_pct_xl(results.fwd_mm_dist.more_mm)
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "  No match: {}",
                fmt_pct_xl(results.fwd_mm_dist.no_match)
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;

    if !probes.is_empty() {
        worksheet
            .write_string_with_format(row, 0, "Probes:", &category_format)
            .map_err(|e| e.to_string())?;
        row += 1;
        worksheet
            .write_string(
                row,
                0,
                &format!(
                    "  0 mismatches: {}",
                    fmt_pct_xl(results.probe_mm_dist.zero_mm)
                ),
            )
            .map_err(|e| e.to_string())?;
        row += 1;
        worksheet
            .write_string(
                row,
                0,
                &format!(
                    "  1 mismatch: {}",
                    fmt_pct_xl(results.probe_mm_dist.one_mm)
                ),
            )
            .map_err(|e| e.to_string())?;
        row += 1;
        worksheet
            .write_string(
                row,
                0,
                &format!(
                    "  >1 mismatches: {}",
                    fmt_pct_xl(results.probe_mm_dist.more_mm)
                ),
            )
            .map_err(|e| e.to_string())?;
        row += 1;
        worksheet
            .write_string(
                row,
                0,
                &format!(
                    "  No match: {}",
                    fmt_pct_xl(results.probe_mm_dist.no_match)
                ),
            )
            .map_err(|e| e.to_string())?;
        row += 1;
    }

    worksheet
        .write_string_with_format(row, 0, "Reverse Primers:", &category_format)
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "  0 mismatches: {}",
                fmt_pct_xl(results.rev_mm_dist.zero_mm)
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!("  1 mismatch: {}", fmt_pct_xl(results.rev_mm_dist.one_mm)),
        )
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "  >1 mismatches: {}",
                fmt_pct_xl(results.rev_mm_dist.more_mm)
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "  No match: {}",
                fmt_pct_xl(results.rev_mm_dist.no_match)
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;

    // Overall pattern distribution
    row += 1;
    worksheet
        .write_string_with_format(
            row,
            0,
            "Overall Pattern (worst best-match across all categories):",
            &category_format,
        )
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "  All categories 0 mismatches: {}",
                fmt_pct_xl(results.overall_all_perfect)
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "  All categories ≤1 mismatch: {}",
                fmt_pct_xl(results.overall_max_one_mm)
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "  ≥2 mismatches in any category: {}",
                fmt_pct_xl(results.overall_two_plus_mm)
            ),
        )
        .map_err(|e| e.to_string())?;
    row += 1;
    worksheet
        .write_string(
            row,
            0,
            &format!(
                "  No match in any category: {}",
                fmt_pct_xl(results.overall_no_match)
            ),
        )
        .map_err(|e| e.to_string())?;

    workbook.save(path).map_err(|e| e.to_string())?;

    Ok(())
}
