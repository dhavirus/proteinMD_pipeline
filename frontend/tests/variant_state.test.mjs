import assert from "node:assert/strict";
import { test } from "node:test";
import { decisionProblems, makeDecision } from "../lib/decisions.js";
import { formatEvidenceValue, rotamerChoices } from "../lib/findings.js";
import { buildManifest } from "../lib/manifest.js";
import { allFindings, createState, decidedIds, setVariants, variantsAreStale } from "../lib/state.js";
import { miniReport, variantManifest } from "./fixtures.mjs";

const opened = () => createState(miniReport(), variantManifest());
const variantFinding = (state) => state.variantFindings[0];
const draft = (patch) => ({ option_id: "choose_rotamer", rationale: "r", decided_by: "Bem", ...patch });

test("variant_build findings from the manifest are reviewed with the audit findings", () => {
  const state = opened();
  assert.equal(allFindings(state).length, miniReport().findings.length + 1);
  assert.deepEqual(state.orphans, []);
  assert.ok(decidedIds(state).has("variant_build/R468Q/A:468"));
});

test("rotamer choices come from the candidate evidence; clashing ones are allowed (relaxation follows)", () => {
  const finding = variantFinding(opened());
  const choices = rotamerChoices(finding);
  assert.ok(choices.some((c) => c.id === "mm-40" && c.clashCount === 0));
  const clashing = choices.find((c) => c.clashCount > 0);
  assert.deepEqual(decisionProblems(finding, draft({ parameters: { rotamer: "mm-40" } })), []);
  assert.deepEqual(decisionProblems(finding, draft({ parameters: { rotamer: clashing.id } })), []);
  assert.match(decisionProblems(finding, draft({ parameters: { rotamer: "zz" } }))[0], /not one of the candidates/);
  const decision = makeDecision(finding, draft({ parameters: { rotamer: "mm-40" } }), "2026-01-01T00:00:00Z");
  assert.deepEqual(decision.parameters, { rotamer: "mm-40" });
});

test("candidate evidence is readable", () => {
  const item = variantFinding(opened()).evidence.find((e) => e.type === "candidate" && e.value === "mm-40");
  assert.match(formatEvidenceValue(item), /^mm-40: chi -65, -65, -40° · library 16 % · 0 clashes/);
});

test("export keeps variants and the variant snapshot; editing variants marks it stale", async () => {
  const state = opened();
  assert.equal(variantsAreStale(state), false);
  const manifest = await buildManifest(state, "2026-01-01T00:00:00Z");
  assert.deepEqual(manifest.variants, variantManifest().variants);
  assert.deepEqual(manifest.variant_snapshot, variantManifest().variant_snapshot);
  assert.equal(manifest.decisions.at(-1).finding_id, "variant_build/R468Q/A:468");
  const edited = setVariants(state, []);
  assert.equal(variantsAreStale(edited), true);
  assert.deepEqual((await buildManifest(edited, "2026-01-01T00:00:00Z")).variants, []);
});

test("a manifest without variants exports without a variants key", async () => {
  const manifest = await buildManifest(createState(miniReport()), "2026-01-01T00:00:00Z");
  assert.equal("variants" in manifest, false);
});
