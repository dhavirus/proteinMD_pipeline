import assert from "node:assert/strict";
import { test } from "node:test";
import { decisionProblems, isoSeconds, makeDecision } from "../lib/decisions.js";
import { altlocChoices, formatEvidenceValue, groupFindings, residuesToShow, severityCounts } from "../lib/findings.js";
import { miniReport } from "./fixtures.mjs";

const findings = miniReport().findings;
const byId = (id) => findings.find((f) => f.id === id);
const draft = (patch) => ({ rationale: "because", decided_by: "Bem", ...patch });

test("a complete draft has no problems and becomes a decision", () => {
  const metal = byId("metals/A:901");
  assert.deepEqual(decisionProblems(metal, draft({ option_id: "bonded_mcpb" })), []);
  const decision = makeDecision(metal, draft({ option_id: "bonded_mcpb" }), "2026-01-01T00:00:00Z");
  assert.deepEqual(Object.keys(decision).sort(), ["decided_by", "finding_id", "option_id", "rationale", "timestamp"]);
});

test("rationale, name and a valid option are required", () => {
  const metal = byId("metals/A:901");
  assert.equal(decisionProblems(metal, { option_id: "nope" }).length, 1);
  assert.equal(decisionProblems(metal, { option_id: "bonded_mcpb", rationale: " ", decided_by: "" }).length, 2);
});

test("explicit-choice options need confirmation", () => {
  const residue = byId("nonstandard_residues/A:3");
  const problems = decisionProblems(residue, draft({ option_id: "revert_to_parent" }));
  assert.match(problems[0], /deliberately/);
  assert.deepEqual(decisionProblems(residue, draft({ option_id: "revert_to_parent", confirmed_explicit: true })), []);
});

test("specific_altloc needs a single-character altloc", () => {
  const water = byId("altlocs/water");
  assert.deepEqual(altlocChoices(water), []);
  assert.match(decisionProblems(water, draft({ option_id: "specific_altloc" }))[0], /Choose which altloc/);
  assert.match(decisionProblems(water, draft({ option_id: "specific_altloc", parameters: { altloc: "AB" } }))[0], /single character/);
  const ok = draft({ option_id: "specific_altloc", parameters: { altloc: "B" } });
  assert.deepEqual(decisionProblems(water, ok), []);
  assert.deepEqual(makeDecision(water, ok, "t").parameters, { altloc: "B" });
});

test("altloc choices come from the evidence", () => {
  const finding = { evidence: [{ key: "altloc_ids", type: "label", value: "A,B" }], options: [{ id: "specific_altloc" }], id: "altlocs/A:1" };
  assert.deepEqual(altlocChoices(finding), ["A", "B"]);
  assert.match(decisionProblems(finding, draft({ option_id: "specific_altloc", parameters: { altloc: "C" } }))[0], /not one of A, B/);
});

test("grouping, filtering and counts", () => {
  assert.deepEqual(severityCounts(findings), { blocking: 1, warn: 3, info: 2 });
  const groups = groupFindings(findings, { severities: null, family: null, undecidedOnly: false, decidedIds: new Set() });
  assert.deepEqual(groups.map((g) => g.family), ["metals", "nonstandard_residues", "altlocs", "missing_residues", "unrecognized"]);
  const blocking = groupFindings(findings, { severities: ["blocking"], family: null, undecidedOnly: false, decidedIds: new Set() });
  assert.deepEqual(blocking.flatMap((g) => g.items.map((f) => f.id)), ["unrecognized/group/A:902"]);
  const undecided = groupFindings(findings, { severities: null, family: null, undecidedOnly: true, decidedIds: new Set(findings.map((f) => f.id)) });
  assert.equal(undecided.length, 0);
});

test("evidence formatting with units", () => {
  assert.equal(formatEvidenceValue({ type: "distance", value: 2.05, unit: "angstrom" }), "2.05 Å");
  assert.equal(formatEvidenceValue({ type: "flag", value: false }), "no");
  assert.equal(formatEvidenceValue({ type: "source_record", category: "c", fields: { id: "x", d: null } }), "id: x");
});

test("timestamps are whole seconds", () => {
  assert.equal(isoSeconds(new Date("2026-01-02T03:04:05.678Z")), "2026-01-02T03:04:05Z");
});

test("contact findings show both atoms of the pair as context", () => {
  const atom = (chain, seq) => ({ chain, seq_num: seq, ins_code: "", res_name: "X", atom_name: "C1", altloc: "" });
  const finding = {
    anchor_residues: [{ chain: "B", seq_num: 1, ins_code: "" }, { chain: "H", seq_num: 314, ins_code: "" }],
    locus: { extent: [] },
    evidence: [{ key: "contact_distance", type: "distance", value: 2.6, unit: "angstrom", atoms: [atom("B", 1), atom("H", 314)] }],
  };
  const { focus, context } = residuesToShow(finding);
  assert.equal(focus.length, 2);
  assert.deepEqual(context.map((r) => `${r.chain}:${r.seq_num}`), ["B:1", "H:314"]);
});
