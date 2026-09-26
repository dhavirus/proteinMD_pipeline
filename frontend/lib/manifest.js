// Assemble manifest.json (manifest.schema.json) from the loaded documents and the state.
import { sha256Canonical } from "./canonical.js";

export const SCHEMA_VERSION = "0.1.0";

/**
 * The manifest to export: the loaded manifest (or a new one from findings.json
 * provenance), with the current regions, decisions and a snapshot of these findings.
 * Orphaned decisions are not exported; the UI lists them before export.
 */
export async function buildManifest(state, createdAt) {
  const { report, baseManifest } = state;
  const base = baseManifest || newManifest(report, createdAt);
  const order = new Map(report.findings.map((f, index) => [f.id, index]));
  const decisions = Object.values(state.decisions).sort((a, b) => order.get(a.finding_id) - order.get(b.finding_id));
  return {
    ...base,
    regions: state.regions,
    findings_snapshot: {
      generated_at: report.generated_at,
      findings_sha256: await sha256Canonical(report.findings),
      findings: report.findings,
    },
    decisions,
  };
}

function newManifest(report, createdAt) {
  return {
    schema_version: SCHEMA_VERSION,
    created_at: createdAt,
    input: report.input,
    knowledge_base: report.knowledge_base,
    simprep: report.simprep,
    regions: [],
    severity_escalation: report.audit_config.severity_escalation,
    findings_snapshot: null,
    decisions: [],
  };
}

/** Why a manifest cannot be used with this findings.json (empty list when it can). */
export function manifestMismatch(manifest, report) {
  const problems = [];
  if (manifest.input.sha256 !== report.input.sha256) {
    problems.push(
      `The manifest was made for input SHA-256 ${manifest.input.sha256}, but findings.json describes ${report.input.sha256}.`,
    );
  }
  return problems;
}
