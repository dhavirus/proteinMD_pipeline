// Application state: one plain object, replaced (never mutated) on every change.
import { isFinal } from "./decisions.js";
import { sameRegions } from "./regions.js";

/**
 * report: a validated findings.json; manifest: a validated manifest or null.
 * The manifest's variant_build findings (variant_snapshot) are reviewed alongside the
 * audit findings. Decisions whose finding is in neither are kept aside as orphans,
 * shown to the user and never exported silently.
 */
export function createState(report, manifest = null) {
  const variantFindings = manifest?.variant_snapshot?.findings || [];
  const ids = new Set([...report.findings, ...variantFindings].map((f) => f.id));
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
    variantFindings,
    variants: manifest?.variants || [],
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
export const setVariants = (state, variants) => ({ ...state, variants });

/** Audit findings followed by the variant_build findings of the loaded manifest. */
export const allFindings = (state) => [...state.report.findings, ...state.variantFindings];

/** The variants were edited since the loaded variant_snapshot was computed. */
export function variantsAreStale(state) {
  const loaded = state.baseManifest?.variants || [];
  return JSON.stringify(loaded) !== JSON.stringify(state.variants);
}
export const dropOrphans = (state) => ({ ...state, orphans: [] });

/** Severities in the report were computed for other regions than the ones now set. */
export function severityIsStale(state) {
  return !sameRegions(state.regions, state.report.audit_config.regions);
}

/** Findings with a final decision (an "expert review" decision leaves a finding open). */
export function decidedIds(state) {
  return new Set(
    allFindings(state)
      .filter((f) => state.decisions[f.id] && isFinal(f, state.decisions[f.id].option_id))
      .map((f) => f.id),
  );
}

/** Restore a saved draft (decisions, regions, variants) for the same input file. */
export function applyDraft(state, draft) {
  const ids = new Set(allFindings(state).map((f) => f.id));
  const decisions = Object.fromEntries(Object.entries(draft.decisions).filter(([id]) => ids.has(id)));
  return { ...state, decisions, regions: draft.regions, variants: draft.variants || state.variants };
}
