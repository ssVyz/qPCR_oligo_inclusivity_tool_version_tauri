use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use tauri::{Emitter, Manager};

use crate::analysis;
use crate::bio_logic;
use crate::excel;
use crate::models::*;

pub struct AppState {
    pub sequences: Mutex<Vec<FastaRecord>>,
    pub sequence_file: Mutex<Option<PathBuf>>,
    pub results: Mutex<Option<AnalysisResults>>,
    pub settings_snapshot: Mutex<AlignmentSettings>,
    pub fwd_primers_snapshot: Mutex<Vec<Oligo>>,
    pub rev_primers_snapshot: Mutex<Vec<Oligo>>,
    pub probes_snapshot: Mutex<Vec<Oligo>>,
    pub running: Arc<AtomicBool>,
}

#[tauri::command]
pub fn load_fasta_file(
    path: String,
    state: tauri::State<'_, AppState>,
) -> Result<FastaFileInfo, String> {
    let path_buf = PathBuf::from(&path);
    let records = bio_logic::parse_fasta(&path_buf)?;

    let num_sequences = records.len();
    let seq_lengths: Vec<usize> = records.iter().map(|r| r.seq.len()).collect();
    let min_length = seq_lengths.iter().copied().min().unwrap_or(0);
    let max_length = seq_lengths.iter().copied().max().unwrap_or(0);
    let avg_length = if !seq_lengths.is_empty() {
        seq_lengths.iter().sum::<usize>() as f64 / seq_lengths.len() as f64
    } else {
        0.0
    };
    let first_ids: Vec<String> = records.iter().take(5).map(|r| r.id.clone()).collect();

    let filename = path_buf
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_default();

    *state.sequences.lock().unwrap() = records;
    *state.sequence_file.lock().unwrap() = Some(path_buf);

    Ok(FastaFileInfo {
        path,
        filename,
        num_sequences,
        min_length,
        max_length,
        avg_length,
        first_ids,
    })
}

#[tauri::command]
pub fn parse_oligos(text: String) -> Result<Vec<OligoInfo>, String> {
    let oligos = bio_logic::parse_fasta_string(&text)?;
    Ok(oligos
        .into_iter()
        .map(|o| OligoInfo {
            id: o.id,
            seq: o.seq,
        })
        .collect())
}

#[tauri::command]
pub fn run_analysis(
    request: AnalysisRequest,
    state: tauri::State<'_, AppState>,
    app: tauri::AppHandle,
) -> Result<(), String> {
    let sequences = state.sequences.lock().unwrap().clone();
    if sequences.is_empty() {
        return Err("No sequences loaded. Please select a FASTA file first.".to_string());
    }

    let fwd_primers = if request.forward_text.trim().is_empty() {
        Vec::new()
    } else {
        bio_logic::parse_fasta_string(&request.forward_text)?
    };

    if fwd_primers.is_empty() {
        return Err("Please enter at least one forward primer.".to_string());
    }

    let rev_primers = if request.reverse_text.trim().is_empty() {
        Vec::new()
    } else {
        bio_logic::parse_fasta_string(&request.reverse_text)?
    };

    if rev_primers.is_empty() {
        return Err("Please enter at least one reverse primer.".to_string());
    }

    let probes = if request.probe_text.trim().is_empty() {
        Vec::new()
    } else {
        bio_logic::parse_fasta_string(&request.probe_text)
            .unwrap_or_default()
    };

    let settings = request.settings;

    // Store snapshots for later Excel export
    *state.settings_snapshot.lock().unwrap() = settings.clone();
    *state.fwd_primers_snapshot.lock().unwrap() = fwd_primers.clone();
    *state.rev_primers_snapshot.lock().unwrap() = rev_primers.clone();
    *state.probes_snapshot.lock().unwrap() = probes.clone();

    // Clear previous results
    *state.results.lock().unwrap() = None;

    let running = state.running.clone();
    running.store(true, Ordering::SeqCst);

    // Get a reference to the results mutex via a raw pointer approach
    // We need to move ownership into the thread
    let results_holder: Arc<Mutex<Option<AnalysisResults>>> = {
        // We'll use the app state indirectly via a channel
        Arc::new(Mutex::new(None))
    };
    let results_for_thread = results_holder.clone();

    let app_clone = app.clone();
    let running_clone = running.clone();

    std::thread::spawn(move || {
        let results = analysis::run_analysis(
            &sequences,
            &fwd_primers,
            &rev_primers,
            &probes,
            &settings,
            &running_clone,
            &app_clone,
        );

        running_clone.store(false, Ordering::SeqCst);

        *results_for_thread.lock().unwrap() = Some(results.clone());

        let payload = AnalysisResultsPayload {
            output_text: results.output_text.clone(),
            total_sequences: results.total_sequences,
            sequences_with_min_matches: results.sequences_with_min_matches,
            has_results: true,
        };

        let _ = app_clone.emit("analysis-complete", payload);
    });

    // We need a way to store results back into AppState after the thread completes.
    // Since we can't hold the State across threads, we'll use a watcher thread.
    let app_for_watcher = app.clone();
    std::thread::spawn(move || {
        // Poll until results are available
        loop {
            std::thread::sleep(std::time::Duration::from_millis(100));
            let maybe_results = results_holder.lock().unwrap().take();
            if let Some(results) = maybe_results {
                // Store into app state
                let state: tauri::State<'_, AppState> = app_for_watcher.state();
                *state.results.lock().unwrap() = Some(results);
                break;
            }
            if !running.load(Ordering::SeqCst) {
                // Analysis was cancelled or completed
                let maybe = results_holder.lock().unwrap().take();
                if let Some(results) = maybe {
                    let state: tauri::State<'_, AppState> = app_for_watcher.state();
                    *state.results.lock().unwrap() = Some(results);
                }
                break;
            }
        }
    });

    Ok(())
}

#[tauri::command]
pub fn cancel_analysis(state: tauri::State<'_, AppState>) -> Result<(), String> {
    state.running.store(false, Ordering::SeqCst);
    Ok(())
}

#[tauri::command]
pub fn save_oligos_json(texts: OligoSetTexts, path: String) -> Result<(), String> {
    let fwd = if texts.forward_text.trim().is_empty() {
        Vec::new()
    } else {
        bio_logic::parse_fasta_string(&texts.forward_text).unwrap_or_default()
    };
    let rev = if texts.reverse_text.trim().is_empty() {
        Vec::new()
    } else {
        bio_logic::parse_fasta_string(&texts.reverse_text).unwrap_or_default()
    };
    let probes = if texts.probe_text.trim().is_empty() {
        Vec::new()
    } else {
        bio_logic::parse_fasta_string(&texts.probe_text).unwrap_or_default()
    };

    let oligo_set = OligoSet {
        forward_primers: fwd
            .iter()
            .map(|o| OligoSerde {
                id: o.id.clone(),
                seq: o.seq.clone(),
            })
            .collect(),
        reverse_primers: rev
            .iter()
            .map(|o| OligoSerde {
                id: o.id.clone(),
                seq: o.seq.clone(),
            })
            .collect(),
        probes: probes
            .iter()
            .map(|o| OligoSerde {
                id: o.id.clone(),
                seq: o.seq.clone(),
            })
            .collect(),
    };

    let json =
        serde_json::to_string_pretty(&oligo_set).map_err(|e| format!("Serialization error: {}", e))?;
    std::fs::write(&path, json).map_err(|e| format!("Error writing file: {}", e))?;
    Ok(())
}

#[tauri::command]
pub fn load_oligos_json(path: String) -> Result<OligoSetTexts, String> {
    let content =
        std::fs::read_to_string(&path).map_err(|e| format!("Error reading file: {}", e))?;
    let oligo_set: OligoSet =
        serde_json::from_str(&content).map_err(|e| format!("Error parsing JSON: {}", e))?;

    let fwd_oligos: Vec<Oligo> = oligo_set
        .forward_primers
        .iter()
        .map(|o| Oligo {
            id: o.id.clone(),
            seq: o.seq.clone(),
        })
        .collect();
    let rev_oligos: Vec<Oligo> = oligo_set
        .reverse_primers
        .iter()
        .map(|o| Oligo {
            id: o.id.clone(),
            seq: o.seq.clone(),
        })
        .collect();
    let probe_oligos: Vec<Oligo> = oligo_set
        .probes
        .iter()
        .map(|o| Oligo {
            id: o.id.clone(),
            seq: o.seq.clone(),
        })
        .collect();

    Ok(OligoSetTexts {
        forward_text: bio_logic::oligos_to_fasta_string(&fwd_oligos),
        reverse_text: bio_logic::oligos_to_fasta_string(&rev_oligos),
        probe_text: bio_logic::oligos_to_fasta_string(&probe_oligos),
    })
}

#[tauri::command]
pub fn export_excel(path: String, state: tauri::State<'_, AppState>) -> Result<(), String> {
    let results = state.results.lock().unwrap();
    let results = results
        .as_ref()
        .ok_or("No results to export. Run analysis first.")?;

    let fwd = state.fwd_primers_snapshot.lock().unwrap().clone();
    let rev = state.rev_primers_snapshot.lock().unwrap().clone();
    let probes = state.probes_snapshot.lock().unwrap().clone();
    let settings = state.settings_snapshot.lock().unwrap().clone();

    let path_buf = PathBuf::from(&path);
    excel::write_excel(&path_buf, results, &fwd, &rev, &probes, &settings)
}

#[tauri::command]
pub fn save_results_text(path: String, state: tauri::State<'_, AppState>) -> Result<(), String> {
    let results = state.results.lock().unwrap();
    let results = results
        .as_ref()
        .ok_or("No results to save. Run analysis first.")?;

    std::fs::write(&path, &results.output_text).map_err(|e| format!("Error writing file: {}", e))
}

#[tauri::command]
pub fn get_fasta_info(state: tauri::State<'_, AppState>) -> Result<String, String> {
    let sequences = state.sequences.lock().unwrap();
    let sequence_file = state.sequence_file.lock().unwrap();

    if sequences.is_empty() {
        return Err("No sequences loaded. Please select a FASTA file first.".to_string());
    }

    let mut info = format!("FASTA File Information:\n{}\n\n", "=".repeat(30));
    info.push_str(&format!("Total sequences: {}\n", sequences.len()));

    if let Some(ref path) = *sequence_file {
        info.push_str(&format!(
            "File: {}\n\n",
            path.file_name()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_default()
        ));
    }

    let seq_lengths: Vec<usize> = sequences.iter().map(|r| r.seq.len()).collect();
    let min_len = seq_lengths.iter().min().unwrap_or(&0);
    let max_len = seq_lengths.iter().max().unwrap_or(&0);
    let avg_len: f64 = seq_lengths.iter().sum::<usize>() as f64 / seq_lengths.len() as f64;

    info.push_str("Sequence length statistics:\n");
    info.push_str(&format!("  Min: {} bp\n", min_len));
    info.push_str(&format!("  Max: {} bp\n", max_len));
    info.push_str(&format!("  Average: {:.1} bp\n\n", avg_len));

    info.push_str("First 5 sequence IDs:\n");
    for (i, record) in sequences.iter().take(5).enumerate() {
        info.push_str(&format!("  {}. {}\n", i + 1, record.id));
    }

    Ok(info)
}

#[tauri::command]
pub fn export_results_json(path: String, state: tauri::State<'_, AppState>) -> Result<(), String> {
    let results = state.results.lock().unwrap();
    let results = results
        .as_ref()
        .ok_or("No results to export. Run analysis first.")?;

    let json = serde_json::to_string_pretty(results)
        .map_err(|e| format!("Serialization error: {}", e))?;
    std::fs::write(&path, json).map_err(|e| format!("Error writing file: {}", e))?;
    Ok(())
}

#[tauri::command]
pub fn dump_fasta(
    path: String,
    mode: String,
    threshold: Option<usize>,
    range_start: Option<usize>,
    range_end: Option<usize>,
    state: tauri::State<'_, AppState>,
) -> Result<String, String> {
    let results = state.results.lock().unwrap();
    let results = results
        .as_ref()
        .ok_or("No results available. Run analysis first.")?;

    let sequences = state.sequences.lock().unwrap();
    if sequences.is_empty() {
        return Err("No sequences loaded.".to_string());
    }

    // Collect IDs from patterns matching the criteria
    let mut matched_ids: HashSet<String> = HashSet::new();
    let mut include_unmatched = false;

    match mode.as_str() {
        "less_than" => {
            let thresh = threshold.ok_or("Threshold required for 'less than' mode.")?;
            for (_sig, data) in &results.alignment_dict {
                if data.total_mismatches < thresh {
                    for id in &data.member_ids {
                        matched_ids.insert(id.clone());
                    }
                }
            }
        }
        "more_than" => {
            let thresh = threshold.ok_or("Threshold required for 'more than' mode.")?;
            for (_sig, data) in &results.alignment_dict {
                if data.total_mismatches > thresh {
                    for id in &data.member_ids {
                        matched_ids.insert(id.clone());
                    }
                }
            }
            include_unmatched = true;
        }
        "unmatched" => {
            include_unmatched = true;
        }
        "range" => {
            let start = range_start.ok_or("Range start required.")?;
            let end = range_end.ok_or("Range end required.")?;
            for (_sig, data) in &results.alignment_dict {
                if data.total_mismatches >= start && data.total_mismatches <= end {
                    for id in &data.member_ids {
                        matched_ids.insert(id.clone());
                    }
                }
            }
        }
        _ => return Err(format!("Unknown dump mode: {}", mode)),
    }

    // If we need unmatched sequences, find all IDs that are NOT in any pattern
    if include_unmatched {
        let mut all_pattern_ids: HashSet<String> = HashSet::new();
        for data in results.alignment_dict.values() {
            for id in &data.member_ids {
                all_pattern_ids.insert(id.clone());
            }
        }
        for record in sequences.iter() {
            if !all_pattern_ids.contains(&record.id) {
                matched_ids.insert(record.id.clone());
            }
        }
    }

    // Write matching sequences to FASTA
    let mut output = String::new();
    let mut count = 0;
    for record in sequences.iter() {
        if matched_ids.contains(&record.id) {
            output.push('>');
            output.push_str(&record.id);
            output.push('\n');
            output.push_str(&record.seq);
            output.push('\n');
            count += 1;
        }
    }

    std::fs::write(&path, &output).map_err(|e| format!("Error writing file: {}", e))?;
    Ok(format!("{} sequences written", count))
}
