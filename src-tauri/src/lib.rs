mod analysis;
mod bio_logic;
mod commands;
mod excel;
mod models;

use commands::AppState;
use std::sync::atomic::AtomicBool;
use std::sync::{Arc, Mutex};

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .plugin(
            tauri_plugin_log::Builder::default()
                .level(log::LevelFilter::Info)
                .build(),
        )
        .manage(AppState {
            sequences: Mutex::new(Vec::new()),
            sequence_file: Mutex::new(None),
            results: Mutex::new(None),
            settings_snapshot: Mutex::new(models::AlignmentSettings::default()),
            fwd_primers_snapshot: Mutex::new(Vec::new()),
            rev_primers_snapshot: Mutex::new(Vec::new()),
            probes_snapshot: Mutex::new(Vec::new()),
            running: Arc::new(AtomicBool::new(false)),
        })
        .invoke_handler(tauri::generate_handler![
            commands::load_fasta_file,
            commands::parse_oligos,
            commands::run_analysis,
            commands::cancel_analysis,
            commands::save_oligos_json,
            commands::load_oligos_json,
            commands::export_excel,
            commands::save_results_text,
            commands::get_fasta_info,
            commands::export_results_json,
            commands::dump_fasta,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
