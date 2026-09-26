import assert from "node:assert/strict";
import { test } from "node:test";
import { buildManifest, manifestMismatch } from "../lib/manifest.js";
import { applyDraft, createState, decidedIds, recordDecision, removeDecision, setRegions, severityIsStale } from "../lib/state.js";
import { finalOption, isFinal } from "../lib/decisions.js";
import { exampleManifest, miniReport } from "./fixtures.mjs";

const decision = (id, option) => ({ finding_id: id, option_id: option, rationale: "r", decided_by: "b", timestamp: "2026-01-01T00:00:00Z" });
const region = { name: "zn", description: "zinc", residues: [{ chain: "A", seq_num: 901, ins_code: "" }] };

test("state starts from the report and is never mutated", () => {
  const report = miniReport();
  const state = createState(report);
  const next = recordDecision(state, decision("metals/A:901", "bonded_mcpb"));
  assert.deepEqual(state.decisions, {});
  assert.equal(Object.keys(next.decisions).length, 1);
  assert.deepEqual(removeDecision(next, "metals/A:901").decisions, {});
});

test("editing regions marks severities stale", () => {
  const state = createState(miniReport());
  assert.ok(!severityIsStale(state));
  assert.ok(severityIsStale(setRegions(state, [region])));
});

test("manifest decisions for unknown findings become orphans", () => {
  const manifest = exampleManifest();
  manifest.decisions.push(decision("metals/Z:1", "bonded_mcpb"));
  const state = createState(miniReport(), manifest);
  assert.deepEqual(state.orphans.map((d) => d.finding_id), ["metals/Z:1"]);
  assert.ok(state.decisions["metals/A:901"]);
});

test("a new manifest carries provenance, regions, decisions and a verified snapshot", async () => {
  const report = miniReport();
  let state = createState(report);
  state = recordDecision(state, decision("unrecognized/group/A:902", "exclude"));
  state = recordDecision(state, decision("metals/A:901", "bonded_mcpb"));
  state = setRegions(state, [region]);
  const manifest = await buildManifest(state, "2026-01-01T00:00:00Z");
  assert.equal(manifest.input.sha256, report.input.sha256);
  assert.equal(manifest.created_at, "2026-01-01T00:00:00Z");
  assert.deepEqual(manifest.decisions.map((d) => d.finding_id), ["metals/A:901", "unrecognized/group/A:902"]);
  assert.deepEqual(manifest.regions, [region]);
  assert.deepEqual(manifest.severity_escalation, report.audit_config.severity_escalation);
  assert.equal(manifest.findings_snapshot.findings.length, report.findings.length);
});

test("a loaded manifest keeps its creation time and thresholds", async () => {
  const manifest = exampleManifest();
  const state = createState(miniReport(), manifest);
  const exported = await buildManifest(state, "2030-01-01T00:00:00Z");
  assert.equal(exported.created_at, manifest.created_at);
  assert.deepEqual(exported.severity_escalation, manifest.severity_escalation);
});

test("a manifest for another input is refused", () => {
  const manifest = exampleManifest();
  manifest.input.sha256 = "f".repeat(64);
  assert.equal(manifestMismatch(manifest, miniReport()).length, 1);
});

test("drafts restore only decisions for findings in this report", () => {
  const state = applyDraft(createState(miniReport()), {
    decisions: { "metals/A:901": decision("metals/A:901", "restrained"), "metals/Z:9": decision("metals/Z:9", "x") },
    regions: [region],
  });
  assert.deepEqual(Object.keys(state.decisions), ["metals/A:901"]);
  assert.deepEqual(state.regions, [region]);
});

test("an expert-review decision leaves the finding undecided", () => {
  const report = miniReport();
  const group = report.findings.find((f) => f.id === "unrecognized/group/A:902");
  group.options = group.options.map((o) => ({ ...o, prep_action: o.id === "expert_review" || o.id === "add_rule" ? "unresolved" : "apply" }));
  assert.ok(!isFinal(group, "expert_review"));
  assert.equal(finalOption(group), "exclude");
  const state = recordDecision(createState(report), decision(group.id, "expert_review"));
  assert.ok(!decidedIds(state).has(group.id));
  assert.ok(decidedIds(recordDecision(state, decision(group.id, "exclude"))).has(group.id));
});
