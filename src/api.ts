import { invoke } from "@tauri-apps/api/core";
import type { FastaFileInfo, OligoInfo, OligoSetTexts, AnalysisRequest } from "./types";

export async function loadFastaFile(path: string): Promise<FastaFileInfo> {
  return invoke<FastaFileInfo>("load_fasta_file", { path });
}

export async function parseOligos(text: string): Promise<OligoInfo[]> {
  return invoke<OligoInfo[]>("parse_oligos", { text });
}

export async function runAnalysis(request: AnalysisRequest): Promise<void> {
  return invoke<void>("run_analysis", { request });
}

export async function cancelAnalysis(): Promise<void> {
  return invoke<void>("cancel_analysis");
}

export async function saveOligosJson(texts: OligoSetTexts, path: string): Promise<void> {
  return invoke<void>("save_oligos_json", { texts, path });
}

export async function loadOligosJson(path: string): Promise<OligoSetTexts> {
  return invoke<OligoSetTexts>("load_oligos_json", { path });
}

export async function exportExcel(path: string): Promise<void> {
  return invoke<void>("export_excel", { path });
}

export async function saveResultsText(path: string): Promise<void> {
  return invoke<void>("save_results_text", { path });
}

export async function getFastaInfo(): Promise<string> {
  return invoke<string>("get_fasta_info");
}

export async function exportResultsJson(path: string): Promise<void> {
  return invoke<void>("export_results_json", { path });
}

export async function dumpFasta(
  path: string,
  mode: string,
  threshold: number | null,
  rangeStart: number | null,
  rangeEnd: number | null
): Promise<string> {
  return invoke<string>("dump_fasta", {
    path,
    mode,
    threshold,
    rangeStart,
    rangeEnd,
  });
}
