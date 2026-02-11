# qPCR Oligo Inclusivity Tool v2.0

A desktop application for evaluating how well a set of qPCR oligonucleotides (forward primers, reverse primers, and probes) matches against a collection of target sequences in FASTA format. It reports per-sequence alignment signatures, mismatch distributions, amplicon statistics, and per-oligo match rates.

## Usage

1. **Load sequences** - Click "Browse" to select a FASTA file containing your target sequences.
2. **Enter oligos** - Use the Forward / Reverse / Probes tabs to paste oligo sequences in FASTA format. Forward and reverse primers are required; probes are optional.
3. **Configure settings** - Adjust minimum matched oligos per category, minimum alignment coverage, and maximum mismatches per oligo. Optionally enable amplicon size constraints.
4. **Run analysis** - Click "Run Analysis". Progress is shown in real time. Results appear in a modal when complete.
5. **Export** - Save results as plain text or Excel (.xlsx). Oligo sets can be saved/loaded as JSON via the File menu.

## How It Works

### Alignment

Each oligo is aligned against every target sequence using pairwise semi-global alignment (rust-bio). The alignment is performed in both sense and antisense orientations. IUPAC ambiguity codes are supported in the scoring function. An oligo is considered matched if it meets the minimum coverage threshold and does not exceed the maximum allowed mismatches.

### Amplicon Detection

For each sequence, the tool finds the best-matching forward and reverse primer pair and checks whether they produce a convergent amplicon (forward on sense strand upstream of reverse on antisense strand, or vice versa). If amplicon size constraints are set, the pair must also fall within the specified range.

### Signature Patterns

For sequences that meet all matching thresholds, a mismatch signature is generated for each oligo alignment. Signatures use `.` for matches, the target base letter for mismatches, and `-` for gaps. Sequences with identical signature patterns across all oligos are grouped and counted.

### Parallelism

Sequence analysis is parallelized across CPU cores using rayon. Progress events are emitted to the frontend every 100 sequences.

## Project Structure

```
src/                        # Frontend (TypeScript + Vite)
  main.ts                   # UI rendering, event wiring, state management
  types.ts                  # TypeScript interfaces matching Rust types
  api.ts                    # Tauri invoke() wrappers
  events.ts                 # Tauri event listeners
  style.css                 # Application styles

src-tauri/src/              # Backend (Rust)
  lib.rs                    # Tauri app setup, plugin registration, state init
  commands.rs               # Tauri command handlers and AppState definition
  models.rs                 # All data structures and payload types
  bio_logic.rs              # FASTA parsing, alignment, IUPAC matching, signatures
  analysis.rs               # Parallel analysis engine, output text generation
  excel.rs                  # Excel export (rust_xlsxwriter)
  main.rs                   # Entry point
```

Communication between frontend and backend uses Tauri commands (request/response via `invoke()`) and Tauri events (backend-to-frontend push for progress updates and results).

## Building from Source

### Prerequisites

- **Rust** >= 1.77.2 (install via [rustup](https://rustup.rs/))
- **Node.js** >= 18 (includes npm)
- **Tauri v2 system dependencies** - On Windows, this means a working MSVC toolchain (Visual Studio Build Tools with the "Desktop development with C++" workload) and WebView2 (pre-installed on Windows 10/11). See [Tauri prerequisites](https://v2.tauri.app/start/prerequisites/) for other platforms.

### Steps

1. Clone the repository:
   ```
   git clone <repo-url>
   cd Aligner_tauriport
   ```

2. Install npm dependencies:
   ```
   npm install
   ```

3. Run in development mode (hot-reloading frontend, debug Rust build):
   ```
   npx tauri dev
   ```

4. Build a release binary and installers:
   ```
   npx tauri build
   ```
   Output will be in `src-tauri/target/release/`. On Windows, this also produces MSI and NSIS installers under `src-tauri/target/release/bundle/`.

### Rust Dependencies

Pulled automatically by Cargo on first build:

- `bio` - Pairwise sequence alignment
- `rayon` - Parallel iteration
- `rust_xlsxwriter` - Excel file generation
- `tauri-plugin-dialog` - Native file open/save dialogs
- `serde` / `serde_json` - Serialization
