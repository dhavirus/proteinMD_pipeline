// Decisions as recorded in manifest.schema.json (decision), validated before they are saved.
import { altlocChoices } from "./findings.js";

/** Options that need parameters, and how to read them. */
export const OPTION_PARAMETERS = { specific_altloc: "altloc" };

/**
 * Validate a draft {option_id, rationale, decided_by, confirmed_explicit, parameters}
 * against its finding; returns a list of problems (empty when it can be saved).
 */
export function decisionProblems(finding, draft) {
  const problems = [];
  const option = finding.options.find((o) => o.id === draft.option_id);
  if (!option) return [`Choose one of the options for ${finding.id}.`];
  if (!(draft.rationale || "").trim()) problems.push("Rationale: say why this option fits the study.");
  if (!(draft.decided_by || "").trim()) problems.push("Decided by: enter your name.");
  if (option.requires_explicit_choice && !draft.confirmed_explicit) {
    problems.push(`"${option.label}" is never a default: confirm that you choose it deliberately.`);
  }
  problems.push(...parameterProblems(finding, draft));
  return problems;
}

function parameterProblems(finding, draft) {
  const key = OPTION_PARAMETERS[draft.option_id];
  if (!key) return [];
  const value = ((draft.parameters || {})[key] || "").trim();
  const choices = altlocChoices(finding);
  if (!value) return ["Choose which altloc to keep."];
  if (choices.length && !choices.includes(value)) return [`Altloc "${value}" is not one of ${choices.join(", ")}.`];
  if (value.length !== 1) return ["An altloc identifier is a single character."];
  return [];
}

/** The decision object to store; call only when decisionProblems() is empty. */
export function makeDecision(finding, draft, timestamp) {
  const decision = {
    finding_id: finding.id,
    option_id: draft.option_id,
    rationale: draft.rationale.trim(),
    decided_by: draft.decided_by.trim(),
    timestamp,
  };
  const key = OPTION_PARAMETERS[draft.option_id];
  if (key) decision.parameters = { [key]: draft.parameters[key].trim() };
  return decision;
}

export function isoSeconds(date) {
  return date.toISOString().replace(/\.\d{3}Z$/, "Z");
}

/** A decision is final unless its option only defers the decision (e.g. expert_review). */
export function isFinal(finding, optionId) {
  const option = finding.options.find((o) => o.id === optionId);
  return Boolean(option) && option.prep_action !== "unresolved";
}

/** The recommended option if it is final, else the first final option. */
export function finalOption(finding) {
  if (isFinal(finding, finding.recommended_option)) return finding.recommended_option;
  const option = finding.options.find((o) => isFinal(finding, o.id));
  return option ? option.id : null;
}
