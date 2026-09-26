// Application state: one plain object, replaced (never mutated) on every change.
import { isFinal } from "./decisions.js";
import { sameRegions } from "./regions.js";

/**
 * report: a validated findings.json; manifest: a validated manifest or null.
 * Decisions from the manifest whose finding is not in this report are kept aside as
 * orphans, shown to the user and never exported silently.
 */
export function createState(report, manifest = null) {
  const ids = new Set(report.findings.map((f) => f.id));
  const decisions = {};
  const orphans = [];
  for (const decision of manifest ? manifest.decisions : []) {
    if (ids.has(decision.finding_id)) decisions[decision.finding_id] = decision;
    else orphans.push(decision);
  }
  return {
    report,
    baseManifest: manifest,
    decisions,
    orphans,
    regions: manifest ? manifest.regions : report.audit_config.regions,
    selectedId: report.findings[0] ? report.findings[0].id : null,
    filter: { severities: null, family: null, undecidedOnly: false },
  };
}

export const select = (state, findingId) => ({ ...state, selectedId: findingId });
export const setFilter = (state, patch) => ({ ...state, filter: { ...state.filter, ...patch } });

export function recordDecision(state, decision) {
  return { ...state, decisions: { ...state.decisions, [decision.finding_id]: decision } };
}

export function removeDecision(state, findingId) {
  const { [findingId]: _removed, ...rest } = state.decisions;
  return { ...state, decisions: rest };
}

export const setRegions = (state, regions) => ({ ...state, regions });
export const dropOrphans = (state) => ({ ...state, orphans: [] });

/** Severities in the report were computed for other regions than the ones now set. */
export function severityIsStale(state) {
  return !sameRegions(state.regions, state.report.audit_config.regions);
}

/** Findings with a final decision (an "expert review" decision leaves a finding open). */
export function decidedIds(state) {
  return new Set(
    state.report.findings
      .filter((f) => state.decisions[f.id] && isFinal(f, state.decisions[f.id].option_id))
      .map((f) => f.id),
  );
}

/** Restore a saved draft (decisions + regions) for the same input file. */
export function applyDraft(state, draft) {
  const ids = new Set(state.report.findings.map((f) => f.id));
  const decisions = Object.fromEntries(Object.entries(draft.decisions).filter(([id]) => ids.has(id)));
  return { ...state, decisions, regions: draft.regions };
}
