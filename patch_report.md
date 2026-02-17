# Patch Report: Total Mismatches Calculation Change

## Date: 2026-02-17

## Intended Behavior

### Previous behavior
"Total Mismatches" per pattern was calculated as the **sum of mismatches from all matched oligos** across all categories (forward primers, reverse primers, probes).

Example with 3 fwd primers (5mm, 2mm, 1mm), 1 rev primer (0mm), 1 probe (1mm):
- Old result: 5 + 2 + 1 + 0 + 1 = **9 total mismatches**

### New behavior
"Total Mismatches" per pattern is now calculated as the **sum of the minimum mismatch count per oligo category**. For each category (forward, reverse, probe), the oligo with the fewest mismatches is selected, and those minimum values are summed.

Same example:
- min(fwd) = min(5, 2, 1) = 1
- min(rev) = min(0) = 0
- min(probe) = min(1) = 1
- New result: 1 + 0 + 1 = **2 total mismatches**

### Rationale
The minimum mismatch per category better represents the actual diagnostic relevance of a pattern. When multiple oligos exist per category (e.g., multiple forward primers), only the best-matching one matters for assay performance. The old sum-of-all approach inflated mismatch counts for assays with many oligos, making patterns appear worse than they functionally are.

### Scope
This change affects the "Total Mismatches" / "Mismatches" value reported **per pattern** in:
- The text report (displayed in the application UI and saved as .txt)
- The Excel export (.xlsx)

The "Overall Pattern" section at the end of the report is **unchanged** — it already used per-category minimum logic independently.

## Technical Implementation

### File changed: `src-tauri/src/bio_logic.rs`

**Function:** `analyze_sequence()` (around line 540)

**What was changed:**
1. The `best_fwd_mm`, `best_rev_mm`, and `best_probe_mm` calculations (per-category minimums) were moved **above** the `total_mismatches` calculation so they can be referenced.
2. The `total_mismatches` formula was replaced.

**Before:**
```rust
let total_mismatches: usize = fwd_results
    .values()
    .chain(rev_results.values())
    .chain(probe_results.values())
    .filter(|r| r.matched)
    .map(|r| r.mismatches)
    .sum();
```

**After:**
```rust
let total_mismatches: usize = best_fwd_mm.unwrap_or(0)
    + best_rev_mm.unwrap_or(0)
    + best_probe_mm.unwrap_or(0);
```

**Why `unwrap_or(0)`:** If no oligo in a category matched (the minimum is `None`), that category contributes 0 mismatches to the total. This is consistent with how the field was used before — unmatched oligos were already excluded from the sum via `.filter(|r| r.matched)`.

### Files NOT changed (but propagate the new value automatically)
- `src-tauri/src/analysis.rs` — `PatternData.total_mismatches` is populated from `seq_result.total_mismatches` (line 195), unchanged.
- `src-tauri/src/analysis.rs` — text report output references `data.total_mismatches` (line 379-380), unchanged.
- `src-tauri/src/excel.rs` — Excel output references `data.total_mismatches` (line 249), unchanged.
- `src-tauri/src/models.rs` — `PatternData` and `SequenceResult` structs unchanged.

### Data flow
```
bio_logic::analyze_sequence()
  -> SequenceResult.total_mismatches  (NEW: sum of per-category minimums)
     -> analysis::run_analysis()
        -> PatternData.total_mismatches  (stored per pattern signature)
           -> generate_output_text()  (text report, line 379)
           -> excel::write_excel()    (Excel, line 249)
```

## Testing Checklist

### Test 1: Basic verification with multiple oligos per category
- **Setup:** Use an assay with 2+ forward primers, 1+ reverse primer, 1+ probe
- **Input:** A FASTA file where sequences have varying mismatch counts across oligos
- **Verify:** "Mismatches" in the text report and "Total Mismatches" in Excel show the sum of per-category minimums, not the sum of all oligos

### Test 2: Single oligo per category
- **Setup:** Use an assay with exactly 1 forward, 1 reverse, 1 probe
- **Expected:** Behavior is identical to old behavior (min of a single value = that value)

### Test 3: Category with no matches
- **Setup:** Use sequences where one category (e.g., probe) has no matches for any oligo
- **Expected:** That category contributes 0 to the total (same as before — unmatched oligos were filtered out)

### Test 4: All perfect matches
- **Setup:** Use sequences where all oligos match perfectly (0 mismatches)
- **Expected:** Total mismatches = 0

### Test 5: Consistency across output formats
- **Verify:** The "Mismatches" value in the text report matches the "Total Mismatches" value in the Excel export for every pattern

### Test 6: Overall Pattern section unchanged
- **Verify:** The "Overall Pattern (worst best-match across all categories)" section at the end of both text and Excel reports is unchanged — it uses its own independent logic based on `best_fwd_mm`, `best_rev_mm`, `best_probe_mm` directly
