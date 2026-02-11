import { listen } from "@tauri-apps/api/event";
import type { ProgressPayload, AnalysisResultsPayload } from "./types";

export function onAnalysisProgress(callback: (payload: ProgressPayload) => void) {
  return listen<ProgressPayload>("analysis-progress", (event) => {
    callback(event.payload);
  });
}

export function onAnalysisComplete(callback: (payload: AnalysisResultsPayload) => void) {
  return listen<AnalysisResultsPayload>("analysis-complete", (event) => {
    callback(event.payload);
  });
}

export function onAnalysisError(callback: (message: string) => void) {
  return listen<string>("analysis-error", (event) => {
    callback(event.payload);
  });
}
