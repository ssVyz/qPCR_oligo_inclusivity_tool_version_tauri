import "./style.css";
import { open, save } from "@tauri-apps/plugin-dialog";
import * as api from "./api";
import { onAnalysisProgress, onAnalysisComplete, onAnalysisError } from "./events";
import type { AlignmentSettings, AnalysisResultsPayload } from "./types";

// ============================================================================
// Application State
// ============================================================================

interface AppState {
  fastaLoaded: boolean;
  fastaFilename: string;
  fastaCount: number;
  activeCategory: "forward" | "reverse" | "probes";
  forwardText: string;
  reverseText: string;
  probeText: string;
  settings: AlignmentSettings;
  ampliconEnabled: boolean;
  ampliconMinInput: string;
  ampliconMaxInput: string;
  isRunning: boolean;
  progressCurrent: number;
  progressTotal: number;
  statusMessage: string;
  hasResults: boolean;
  resultsText: string;
}

const state: AppState = {
  fastaLoaded: false,
  fastaFilename: "",
  fastaCount: 0,
  activeCategory: "forward",
  forwardText: "",
  reverseText: "",
  probeText: "",
  settings: {
    mode: "Local",
    match_score: 2,
    mismatch_score: -1,
    gap_open_score: -2,
    gap_extend_score: -1,
    min_fwd_matched: 1,
    min_rev_matched: 1,
    min_probe_matched: 0,
    min_coverage: 0.8,
    max_mismatches_per_oligo: 7,
    ambiguity_display: "ShowDots",
    min_amplicon_size: null,
    max_amplicon_size: null,
  },
  ampliconEnabled: false,
  ampliconMinInput: "",
  ampliconMaxInput: "1000",
  isRunning: false,
  progressCurrent: 0,
  progressTotal: 0,
  statusMessage: "Ready",
  hasResults: false,
  resultsText: "",
};

// ============================================================================
// Render HTML
// ============================================================================

function renderApp(): string {
  return `
    <div class="app-container">
      <header class="menu-bar">
        <div class="menu-item" id="menu-file">
          File
          <div class="menu-dropdown" id="dropdown-file">
            <button id="btn-select-file">Select Sequence File...</button>
            <div class="menu-separator"></div>
            <button id="btn-export-excel">Export to Excel...</button>
          </div>
        </div>
        <div class="menu-item" id="menu-tools">
          Tools
          <div class="menu-dropdown" id="dropdown-tools">
            <button id="btn-fasta-info-menu">FASTA Info</button>
            <button id="btn-advanced-settings-menu">Advanced Settings...</button>
          </div>
        </div>
        <div class="menu-item" id="menu-help">
          Help
          <div class="menu-dropdown" id="dropdown-help">
            <button id="btn-about-menu">About</button>
          </div>
        </div>
      </header>

      <main class="content">
        <h1>qPCR Oligo Inclusivity Tool</h1>
        <p class="subtitle">Categorized oligo analysis with amplicon detection</p>

        <!-- Sequence File Section -->
        <section class="panel">
          <h2>Sequence File</h2>
          <div class="row">
            <button id="btn-browse" class="btn">Browse...</button>
            <span id="file-status" class="file-status">No file selected</span>
          </div>
        </section>

        <!-- Oligo Sequences Section -->
        <section class="panel">
          <h2>Oligo Sequences</h2>
          <div class="row oligo-controls">
            <div class="tab-group">
              <button class="tab-btn active" data-category="forward" id="tab-forward">Forward Primers</button>
              <button class="tab-btn" data-category="reverse" id="tab-reverse">Reverse Primers</button>
              <button class="tab-btn" data-category="probes" id="tab-probes">Probes</button>
            </div>
            <div class="oligo-actions">
              <button id="btn-save-json" class="btn btn-small">Save JSON</button>
              <button id="btn-load-json" class="btn btn-small">Load JSON</button>
            </div>
          </div>
          <div class="category-label" id="category-label">Forward Primers (FASTA format)</div>
          <textarea id="oligo-textarea" class="oligo-textarea" spellcheck="false"
            placeholder=">Primer1&#10;ATGCGTACGTAGC&#10;>Primer2&#10;GCTAGCTAGCTA"></textarea>
          <div id="oligo-count" class="oligo-count"></div>
        </section>

        <!-- Quick Settings Section -->
        <section class="panel">
          <h2>Quick Settings</h2>
          <div class="settings-grid">
            <div class="setting-row">
              <label>Min Fwd Matched:</label>
              <input type="number" id="setting-min-fwd" min="0" max="100" value="1" class="setting-input">
              <label>Min Rev Matched:</label>
              <input type="number" id="setting-min-rev" min="0" max="100" value="1" class="setting-input">
              <label>Min Probes Matched:</label>
              <input type="number" id="setting-min-probe" min="0" max="100" value="0" class="setting-input">
            </div>
            <div class="setting-row">
              <label>Min Coverage:</label>
              <input type="number" id="setting-min-coverage" min="0" max="1" step="0.01" value="0.8" class="setting-input">
              <label>Max Mismatches/Oligo:</label>
              <input type="number" id="setting-max-mm" min="0" max="50" value="7" class="setting-input">
              <label>Ambiguity Display:</label>
              <select id="setting-ambiguity" class="setting-input">
                <option value="ShowDots">Show Dots</option>
                <option value="ShowBases">Show Base Letters</option>
              </select>
            </div>
          </div>
        </section>

        <!-- Amplicon Constraints Section -->
        <section class="panel">
          <h2>Amplicon Size Constraints</h2>
          <div class="row">
            <label class="checkbox-label">
              <input type="checkbox" id="amplicon-enabled">
              Enable amplicon size limits
            </label>
          </div>
          <div id="amplicon-inputs" class="amplicon-inputs hidden">
            <div class="setting-row">
              <label>Min amplicon size (bp):</label>
              <input type="text" id="amplicon-min" value="" class="setting-input" placeholder="(none)">
              <label>Max amplicon size (bp):</label>
              <input type="text" id="amplicon-max" value="1000" class="setting-input">
            </div>
          </div>
          <p id="amplicon-hint" class="hint">Amplicon detection is always active (fwd+rev pairing). This only constrains the allowed size range.</p>
        </section>

        <!-- Action Buttons -->
        <section class="actions-row">
          <button id="btn-run" class="btn btn-primary" disabled>Run Analysis</button>
          <button id="btn-stop" class="btn btn-danger hidden">Stop</button>
          <button id="btn-fasta-info" class="btn">FASTA Info</button>
          <button id="btn-advanced-settings" class="btn">Advanced Settings</button>
        </section>

        <!-- Progress Section -->
        <section class="panel">
          <h2>Progress</h2>
          <div class="progress-container">
            <div class="progress-bar">
              <div class="progress-fill" id="progress-fill" style="width: 0%"></div>
            </div>
            <div id="progress-text" class="progress-text">0%</div>
          </div>
          <div id="progress-detail" class="progress-detail"></div>
        </section>
      </main>

      <footer class="status-bar">
        <span id="status-message">Ready</span>
      </footer>
    </div>

    <!-- Results Modal -->
    <div class="modal-overlay hidden" id="modal-results">
      <div class="modal modal-large">
        <div class="modal-header">
          <h2>Analysis Results</h2>
          <div class="modal-actions">
            <button id="btn-save-text" class="btn btn-small">Save Text</button>
            <button id="btn-export-excel-results" class="btn btn-small">Export Excel</button>
            <button id="btn-close-results" class="btn btn-small">Close</button>
          </div>
        </div>
        <textarea id="results-textarea" class="results-textarea" readonly></textarea>
      </div>
    </div>

    <!-- Advanced Settings Modal -->
    <div class="modal-overlay hidden" id="modal-advanced">
      <div class="modal">
        <h2>Advanced Settings</h2>
        <div class="adv-section">
          <h3>Alignment Parameters</h3>
          <div class="adv-grid">
            <label>Alignment Mode:</label>
            <select id="adv-mode" class="setting-input">
              <option value="Local">Local</option>
              <option value="Global">Global</option>
            </select>
            <label>Match Score:</label>
            <input type="number" id="adv-match" min="-10" max="10" value="2" class="setting-input">
            <label>Mismatch Score:</label>
            <input type="number" id="adv-mismatch" min="-10" max="10" value="-1" class="setting-input">
            <label>Gap Open Score:</label>
            <input type="number" id="adv-gap-open" min="-20" max="0" value="-2" class="setting-input">
            <label>Gap Extend Score:</label>
            <input type="number" id="adv-gap-extend" min="-10" max="0" value="-1" class="setting-input">
            <label>Min Coverage:</label>
            <input type="number" id="adv-coverage" min="0" max="1" step="0.01" value="0.8" class="setting-input">
            <label>Max Mismatches/Oligo:</label>
            <input type="number" id="adv-max-mm" min="0" max="50" value="7" class="setting-input">
            <label>Ambiguity Display:</label>
            <select id="adv-ambiguity" class="setting-input">
              <option value="ShowDots">Show Dots</option>
              <option value="ShowBases">Show Base Letters</option>
            </select>
          </div>
        </div>
        <div class="adv-section">
          <h3>Per-Category Match Thresholds</h3>
          <div class="adv-grid">
            <label>Min Forward Primers Matched:</label>
            <input type="number" id="adv-min-fwd" min="0" max="100" value="1" class="setting-input">
            <label>Min Reverse Primers Matched:</label>
            <input type="number" id="adv-min-rev" min="0" max="100" value="1" class="setting-input">
            <label>Min Probes Matched:</label>
            <input type="number" id="adv-min-probe" min="0" max="100" value="0" class="setting-input">
          </div>
        </div>
        <div class="adv-section">
          <h3>Amplicon Size Constraints</h3>
          <label class="checkbox-label">
            <input type="checkbox" id="adv-amplicon-enabled">
            Enable amplicon size limits
          </label>
          <div id="adv-amplicon-inputs" class="setting-row" style="margin-top:8px">
            <label>Min (bp):</label>
            <input type="text" id="adv-amplicon-min" class="setting-input" placeholder="(none)">
            <label>Max (bp):</label>
            <input type="text" id="adv-amplicon-max" value="1000" class="setting-input">
          </div>
        </div>
        <div class="modal-footer">
          <button id="btn-reset-defaults" class="btn">Reset to Defaults</button>
          <button id="btn-close-advanced" class="btn btn-primary">OK</button>
        </div>
      </div>
    </div>

    <!-- FASTA Info Modal -->
    <div class="modal-overlay hidden" id="modal-fasta-info">
      <div class="modal">
        <h2>FASTA Information</h2>
        <textarea id="fasta-info-text" class="info-textarea" readonly></textarea>
        <div class="modal-footer">
          <button id="btn-close-fasta-info" class="btn">Close</button>
        </div>
      </div>
    </div>

    <!-- About Modal -->
    <div class="modal-overlay hidden" id="modal-about">
      <div class="modal modal-small">
        <h2>qPCR Oligo Inclusivity Tool</h2>
        <p>Version 2.0.0</p>
        <p>A high-performance oligo inclusivity analysis tool built with Rust and rust-bio.</p>
        <h3>Features:</h3>
        <ul>
          <li>Categorized oligo input (Forward/Reverse/Probe)</li>
          <li>Amplicon-aware analysis with fwd+rev pairing</li>
          <li>Probe matching within amplicon regions</li>
          <li>Bidirectional primer alignment (sense/antisense)</li>
          <li>Per-category match thresholds</li>
          <li>Parallel processing with Rayon</li>
          <li>Excel export with categorized columns</li>
          <li>JSON save/load for oligo sets</li>
        </ul>
        <div class="modal-footer">
          <button id="btn-close-about" class="btn btn-primary">OK</button>
        </div>
      </div>
    </div>
  `;
}

// ============================================================================
// Initialization
// ============================================================================

function $(id: string): HTMLElement {
  return document.getElementById(id)!;
}

function init() {
  document.querySelector<HTMLDivElement>("#app")!.innerHTML = renderApp();

  setupMenus();
  setupFileActions();
  setupOligoTabs();
  setupOligoActions();
  setupSettingsSync();
  setupAmpliconToggle();
  setupAnalysisControls();
  setupResultsModal();
  setupAdvancedSettingsModal();
  setupFastaInfoModal();
  setupAboutModal();
  setupBackendEvents();

  updateRunButton();
}

// ============================================================================
// Menu Dropdowns
// ============================================================================

function setupMenus() {
  const menuItems = document.querySelectorAll(".menu-item");
  menuItems.forEach((item) => {
    item.addEventListener("click", (e) => {
      e.stopPropagation();
      const wasOpen = item.classList.contains("open");
      closeAllMenus();
      if (!wasOpen) {
        item.classList.add("open");
      }
    });
  });

  document.addEventListener("click", () => {
    closeAllMenus();
  });
}

function closeAllMenus() {
  document.querySelectorAll(".menu-item").forEach((m) => m.classList.remove("open"));
}

// ============================================================================
// File Actions
// ============================================================================

function setupFileActions() {
  $("btn-browse").addEventListener("click", browseFile);
  $("btn-select-file").addEventListener("click", () => {
    closeAllMenus();
    browseFile();
  });
}

async function browseFile() {
  const path = await open({
    filters: [
      { name: "FASTA files", extensions: ["fasta", "fas", "fa", "fna"] },
      { name: "All files", extensions: ["*"] },
    ],
    title: "Select Sequence File",
  });

  if (path) {
    try {
      const info = await api.loadFastaFile(path as string);
      state.fastaLoaded = true;
      state.fastaFilename = info.filename;
      state.fastaCount = info.num_sequences;
      $("file-status").textContent = `${info.filename} (${info.num_sequences} sequences)`;
      $("file-status").classList.add("loaded");
      setStatus(`Loaded ${info.num_sequences} sequences`);
      updateRunButton();
    } catch (e) {
      setStatus(`Error loading sequences: ${e}`);
    }
  }
}

// ============================================================================
// Oligo Tab Management
// ============================================================================

function setupOligoTabs() {
  const tabs = document.querySelectorAll(".tab-btn");
  tabs.forEach((tab) => {
    tab.addEventListener("click", () => {
      const category = (tab as HTMLElement).dataset.category as "forward" | "reverse" | "probes";
      switchCategory(category);
    });
  });

  const textarea = $("oligo-textarea") as HTMLTextAreaElement;
  textarea.addEventListener("input", () => {
    saveCurrentCategoryText();
    updateOligoCount();
    updateTabIndicators();
    updateRunButton();
  });
}

function switchCategory(category: "forward" | "reverse" | "probes") {
  saveCurrentCategoryText();
  state.activeCategory = category;

  document.querySelectorAll(".tab-btn").forEach((t) => t.classList.remove("active"));
  $(`tab-${category}`).classList.add("active");

  const textarea = $("oligo-textarea") as HTMLTextAreaElement;
  const labels: Record<string, string> = {
    forward: "Forward Primers (FASTA format)",
    reverse: "Reverse Primers (FASTA format)",
    probes: "Probes (FASTA format)",
  };
  $("category-label").textContent = labels[category];

  const texts: Record<string, string> = {
    forward: state.forwardText,
    reverse: state.reverseText,
    probes: state.probeText,
  };
  textarea.value = texts[category];
  updateOligoCount();
}

function saveCurrentCategoryText() {
  const textarea = $("oligo-textarea") as HTMLTextAreaElement;
  switch (state.activeCategory) {
    case "forward":
      state.forwardText = textarea.value;
      break;
    case "reverse":
      state.reverseText = textarea.value;
      break;
    case "probes":
      state.probeText = textarea.value;
      break;
  }
}

function updateOligoCount() {
  const textarea = $("oligo-textarea") as HTMLTextAreaElement;
  const text = textarea.value;
  const count = text.split("\n").filter((line) => line.startsWith(">")).length;
  const el = $("oligo-count");
  if (count > 0) {
    el.textContent = `${count} sequence(s) defined`;
    el.classList.add("has-oligos");
  } else {
    el.textContent = "";
    el.classList.remove("has-oligos");
  }
}

function updateTabIndicators() {
  const countSeqs = (text: string) => text.split("\n").filter((l) => l.startsWith(">")).length;

  const fwdOk = countSeqs(state.forwardText) > 0;
  const revOk = countSeqs(state.reverseText) > 0;
  const probeOk = countSeqs(state.probeText) > 0;

  $("tab-forward").classList.toggle("valid", fwdOk);
  $("tab-forward").classList.toggle("invalid", !fwdOk);
  $("tab-reverse").classList.toggle("valid", revOk);
  $("tab-reverse").classList.toggle("invalid", !revOk);
  $("tab-probes").classList.toggle("valid", probeOk);
  // Probes are optional, so don't mark as invalid
  if (!probeOk) {
    $("tab-probes").classList.remove("invalid");
    $("tab-probes").classList.remove("valid");
  }
}

// ============================================================================
// Oligo JSON Save/Load
// ============================================================================

function setupOligoActions() {
  $("btn-save-json").addEventListener("click", saveOligosJson);
  $("btn-load-json").addEventListener("click", loadOligosJson);
}

async function saveOligosJson() {
  saveCurrentCategoryText();
  const path = await save({
    filters: [{ name: "JSON files", extensions: ["json"] }],
    title: "Save Oligo Set",
  });
  if (path) {
    try {
      await api.saveOligosJson(
        {
          forward_text: state.forwardText,
          reverse_text: state.reverseText,
          probe_text: state.probeText,
        },
        path
      );
      setStatus("Oligo set saved");
    } catch (e) {
      setStatus(`Error saving oligo set: ${e}`);
    }
  }
}

async function loadOligosJson() {
  const path = await open({
    filters: [{ name: "JSON files", extensions: ["json"] }],
    title: "Load Oligo Set",
  });
  if (path) {
    try {
      const texts = await api.loadOligosJson(path as string);
      state.forwardText = texts.forward_text;
      state.reverseText = texts.reverse_text;
      state.probeText = texts.probe_text;

      const textarea = $("oligo-textarea") as HTMLTextAreaElement;
      const textMap: Record<string, string> = {
        forward: state.forwardText,
        reverse: state.reverseText,
        probes: state.probeText,
      };
      textarea.value = textMap[state.activeCategory];
      updateOligoCount();
      updateTabIndicators();
      updateRunButton();

      const fwdCount = state.forwardText.split("\n").filter((l) => l.startsWith(">")).length;
      const revCount = state.reverseText.split("\n").filter((l) => l.startsWith(">")).length;
      const probeCount = state.probeText.split("\n").filter((l) => l.startsWith(">")).length;
      setStatus(`Loaded ${fwdCount} fwd, ${revCount} rev, ${probeCount} probes`);
    } catch (e) {
      setStatus(`Error loading oligo set: ${e}`);
    }
  }
}

// ============================================================================
// Settings Sync
// ============================================================================

function setupSettingsSync() {
  const syncNum = (id: string, key: keyof AlignmentSettings) => {
    $(id).addEventListener("change", () => {
      const val = ($(id) as HTMLInputElement).value;
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      (state.settings as any)[key] = Number(val);
    });
  };
  const syncFloat = (id: string, key: keyof AlignmentSettings) => {
    $(id).addEventListener("change", () => {
      const val = ($(id) as HTMLInputElement).value;
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      (state.settings as any)[key] = parseFloat(val);
    });
  };
  const syncStr = (id: string, key: keyof AlignmentSettings) => {
    $(id).addEventListener("change", () => {
      const val = ($(id) as HTMLInputElement | HTMLSelectElement).value;
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      (state.settings as any)[key] = val;
    });
  };

  syncNum("setting-min-fwd", "min_fwd_matched");
  syncNum("setting-min-rev", "min_rev_matched");
  syncNum("setting-min-probe", "min_probe_matched");
  syncFloat("setting-min-coverage", "min_coverage");
  syncNum("setting-max-mm", "max_mismatches_per_oligo");
  syncStr("setting-ambiguity", "ambiguity_display");
}

// ============================================================================
// Amplicon Toggle
// ============================================================================

function setupAmpliconToggle() {
  $("amplicon-enabled").addEventListener("change", () => {
    const checked = ($("amplicon-enabled") as HTMLInputElement).checked;
    state.ampliconEnabled = checked;
    $("amplicon-inputs").classList.toggle("hidden", !checked);
    $("amplicon-hint").textContent = checked
      ? "Limits the amplicon size for valid fwd+rev primer pairs. Without size limits, any convergent pair is accepted."
      : "Amplicon detection is always active (fwd+rev pairing). This only constrains the allowed size range.";
  });
}

// ============================================================================
// Analysis Controls
// ============================================================================

function setupAnalysisControls() {
  $("btn-run").addEventListener("click", startAnalysis);
  $("btn-stop").addEventListener("click", stopAnalysis);
  $("btn-fasta-info").addEventListener("click", showFastaInfo);
  $("btn-fasta-info-menu").addEventListener("click", () => {
    closeAllMenus();
    showFastaInfo();
  });
  $("btn-advanced-settings").addEventListener("click", () => openAdvancedSettings());
  $("btn-advanced-settings-menu").addEventListener("click", () => {
    closeAllMenus();
    openAdvancedSettings();
  });
}

function updateRunButton() {
  const fwdOk = state.forwardText.split("\n").some((l) => l.startsWith(">"));
  const revOk = state.reverseText.split("\n").some((l) => l.startsWith(">"));
  const canRun = state.fastaLoaded && fwdOk && revOk && !state.isRunning;

  const btn = $("btn-run") as HTMLButtonElement;
  btn.disabled = !canRun;
  btn.classList.toggle("ready", canRun);
}

function collectSettings(): AlignmentSettings {
  const s = { ...state.settings };

  s.min_fwd_matched = Number(($("setting-min-fwd") as HTMLInputElement).value);
  s.min_rev_matched = Number(($("setting-min-rev") as HTMLInputElement).value);
  s.min_probe_matched = Number(($("setting-min-probe") as HTMLInputElement).value);
  s.min_coverage = parseFloat(($("setting-min-coverage") as HTMLInputElement).value);
  s.max_mismatches_per_oligo = Number(($("setting-max-mm") as HTMLInputElement).value);
  s.ambiguity_display = ($("setting-ambiguity") as HTMLSelectElement).value as "ShowDots" | "ShowBases";

  if (state.ampliconEnabled) {
    const minVal = ($("amplicon-min") as HTMLInputElement).value.trim();
    const maxVal = ($("amplicon-max") as HTMLInputElement).value.trim();
    s.min_amplicon_size = minVal ? Number(minVal) || null : null;
    s.max_amplicon_size = maxVal ? Number(maxVal) || null : null;
  } else {
    s.min_amplicon_size = null;
    s.max_amplicon_size = null;
  }

  return s;
}

async function startAnalysis() {
  saveCurrentCategoryText();
  const settings = collectSettings();
  state.settings = settings;

  state.isRunning = true;
  state.progressCurrent = 0;
  state.progressTotal = 0;
  updateRunButton();
  $("btn-stop").classList.remove("hidden");
  updateProgress(0, 1);
  setStatus("Starting analysis...");

  try {
    await api.runAnalysis({
      forward_text: state.forwardText,
      reverse_text: state.reverseText,
      probe_text: state.probeText,
      settings,
    });
  } catch (e) {
    state.isRunning = false;
    $("btn-stop").classList.add("hidden");
    updateRunButton();
    setStatus(`Error: ${e}`);
  }
}

async function stopAnalysis() {
  try {
    await api.cancelAnalysis();
  } catch (_e) {
    // ignore
  }
  state.isRunning = false;
  $("btn-stop").classList.add("hidden");
  updateRunButton();
  setStatus("Analysis cancelled");
}

function updateProgress(current: number, total: number) {
  const pct = total > 0 ? (current / total) * 100 : 0;
  $("progress-fill").style.width = `${pct}%`;
  $("progress-text").textContent = `${Math.round(pct)}%`;
  if (total > 0) {
    $("progress-detail").textContent = `${current} / ${total} sequences processed`;
  }
}

// ============================================================================
// Backend Events
// ============================================================================

function setupBackendEvents() {
  onAnalysisProgress((payload) => {
    state.progressCurrent = payload.current;
    state.progressTotal = payload.total;
    updateProgress(payload.current, payload.total);
  });

  onAnalysisComplete((payload: AnalysisResultsPayload) => {
    state.isRunning = false;
    state.hasResults = true;
    state.resultsText = payload.output_text;
    $("btn-stop").classList.add("hidden");
    updateRunButton();
    updateProgress(state.progressTotal, state.progressTotal);
    setStatus("Analysis complete!");
    showResults();
  });

  onAnalysisError((message: string) => {
    state.isRunning = false;
    $("btn-stop").classList.add("hidden");
    updateRunButton();
    setStatus(`Analysis error: ${message}`);
  });
}

// ============================================================================
// Results Modal
// ============================================================================

function setupResultsModal() {
  $("btn-close-results").addEventListener("click", () => {
    $("modal-results").classList.add("hidden");
  });

  $("btn-save-text").addEventListener("click", async () => {
    const path = await save({
      filters: [{ name: "Text files", extensions: ["txt"] }],
      title: "Save Results",
    });
    if (path) {
      try {
        await api.saveResultsText(path);
        setStatus("Results saved");
      } catch (e) {
        setStatus(`Error saving: ${e}`);
      }
    }
  });

  $("btn-export-excel-results").addEventListener("click", doExportExcel);
  $("btn-export-excel").addEventListener("click", () => {
    closeAllMenus();
    doExportExcel();
  });
}

function showResults() {
  ($("results-textarea") as HTMLTextAreaElement).value = state.resultsText;
  $("modal-results").classList.remove("hidden");
}

async function doExportExcel() {
  if (!state.hasResults) {
    setStatus("No results to export. Run analysis first.");
    return;
  }
  const path = await save({
    filters: [{ name: "Excel files", extensions: ["xlsx"] }],
    title: "Save Excel Results As",
  });
  if (path) {
    try {
      await api.exportExcel(path);
      setStatus("Excel export complete");
    } catch (e) {
      setStatus(`Error exporting Excel: ${e}`);
    }
  }
}

// ============================================================================
// Advanced Settings Modal
// ============================================================================

function setupAdvancedSettingsModal() {
  $("btn-close-advanced").addEventListener("click", () => {
    applyAdvancedSettings();
    $("modal-advanced").classList.add("hidden");
  });

  $("btn-reset-defaults").addEventListener("click", resetDefaults);

  $("modal-advanced").addEventListener("click", (e) => {
    if (e.target === $("modal-advanced")) {
      applyAdvancedSettings();
      $("modal-advanced").classList.add("hidden");
    }
  });
}

function openAdvancedSettings() {
  const s = state.settings;
  ($("adv-mode") as HTMLSelectElement).value = s.mode;
  ($("adv-match") as HTMLInputElement).value = String(s.match_score);
  ($("adv-mismatch") as HTMLInputElement).value = String(s.mismatch_score);
  ($("adv-gap-open") as HTMLInputElement).value = String(s.gap_open_score);
  ($("adv-gap-extend") as HTMLInputElement).value = String(s.gap_extend_score);
  ($("adv-coverage") as HTMLInputElement).value = String(s.min_coverage);
  ($("adv-max-mm") as HTMLInputElement).value = String(s.max_mismatches_per_oligo);
  ($("adv-ambiguity") as HTMLSelectElement).value = s.ambiguity_display;
  ($("adv-min-fwd") as HTMLInputElement).value = String(s.min_fwd_matched);
  ($("adv-min-rev") as HTMLInputElement).value = String(s.min_rev_matched);
  ($("adv-min-probe") as HTMLInputElement).value = String(s.min_probe_matched);
  ($("adv-amplicon-enabled") as HTMLInputElement).checked = state.ampliconEnabled;
  ($("adv-amplicon-min") as HTMLInputElement).value = state.ampliconMinInput;
  ($("adv-amplicon-max") as HTMLInputElement).value = state.ampliconMaxInput;

  $("modal-advanced").classList.remove("hidden");
}

function applyAdvancedSettings() {
  state.settings.mode = ($("adv-mode") as HTMLSelectElement).value as "Local" | "Global";
  state.settings.match_score = Number(($("adv-match") as HTMLInputElement).value);
  state.settings.mismatch_score = Number(($("adv-mismatch") as HTMLInputElement).value);
  state.settings.gap_open_score = Number(($("adv-gap-open") as HTMLInputElement).value);
  state.settings.gap_extend_score = Number(($("adv-gap-extend") as HTMLInputElement).value);
  state.settings.min_coverage = parseFloat(($("adv-coverage") as HTMLInputElement).value);
  state.settings.max_mismatches_per_oligo = Number(($("adv-max-mm") as HTMLInputElement).value);
  state.settings.ambiguity_display = ($("adv-ambiguity") as HTMLSelectElement).value as "ShowDots" | "ShowBases";
  state.settings.min_fwd_matched = Number(($("adv-min-fwd") as HTMLInputElement).value);
  state.settings.min_rev_matched = Number(($("adv-min-rev") as HTMLInputElement).value);
  state.settings.min_probe_matched = Number(($("adv-min-probe") as HTMLInputElement).value);

  state.ampliconEnabled = ($("adv-amplicon-enabled") as HTMLInputElement).checked;
  state.ampliconMinInput = ($("adv-amplicon-min") as HTMLInputElement).value;
  state.ampliconMaxInput = ($("adv-amplicon-max") as HTMLInputElement).value;

  // Sync back to quick settings UI
  ($("setting-min-fwd") as HTMLInputElement).value = String(state.settings.min_fwd_matched);
  ($("setting-min-rev") as HTMLInputElement).value = String(state.settings.min_rev_matched);
  ($("setting-min-probe") as HTMLInputElement).value = String(state.settings.min_probe_matched);
  ($("setting-min-coverage") as HTMLInputElement).value = String(state.settings.min_coverage);
  ($("setting-max-mm") as HTMLInputElement).value = String(state.settings.max_mismatches_per_oligo);
  ($("setting-ambiguity") as HTMLSelectElement).value = state.settings.ambiguity_display;
  ($("amplicon-enabled") as HTMLInputElement).checked = state.ampliconEnabled;
  $("amplicon-inputs").classList.toggle("hidden", !state.ampliconEnabled);
  ($("amplicon-min") as HTMLInputElement).value = state.ampliconMinInput;
  ($("amplicon-max") as HTMLInputElement).value = state.ampliconMaxInput;
}

function resetDefaults() {
  state.settings = {
    mode: "Local",
    match_score: 2,
    mismatch_score: -1,
    gap_open_score: -2,
    gap_extend_score: -1,
    min_fwd_matched: 1,
    min_rev_matched: 1,
    min_probe_matched: 0,
    min_coverage: 0.8,
    max_mismatches_per_oligo: 7,
    ambiguity_display: "ShowDots",
    min_amplicon_size: null,
    max_amplicon_size: null,
  };
  state.ampliconEnabled = false;
  state.ampliconMinInput = "";
  state.ampliconMaxInput = "1000";
  openAdvancedSettings();
}

// ============================================================================
// FASTA Info Modal
// ============================================================================

function setupFastaInfoModal() {
  $("btn-close-fasta-info").addEventListener("click", () => {
    $("modal-fasta-info").classList.add("hidden");
  });
}

async function showFastaInfo() {
  try {
    const info = await api.getFastaInfo();
    ($("fasta-info-text") as HTMLTextAreaElement).value = info;
    $("modal-fasta-info").classList.remove("hidden");
  } catch (e) {
    setStatus(`${e}`);
  }
}

// ============================================================================
// About Modal
// ============================================================================

function setupAboutModal() {
  $("btn-about-menu").addEventListener("click", () => {
    closeAllMenus();
    $("modal-about").classList.remove("hidden");
  });

  $("btn-close-about").addEventListener("click", () => {
    $("modal-about").classList.add("hidden");
  });
}

// ============================================================================
// Status Bar
// ============================================================================

function setStatus(msg: string) {
  state.statusMessage = msg;
  $("status-message").textContent = msg;
}

// ============================================================================
// Start
// ============================================================================

init();
