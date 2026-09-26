// Contract fixture: build manifests with the front end's own modules, for pytest to check
// with simprep's Python validator and CLI (tests/test_frontend_contract.py).
// Usage: node export_manifests.mjs FINDINGS_JSON OUT_DIR "A:45 A:46 ..."
import { readFileSync, writeFileSync } from "node:fs";
import { join } from "node:path";
import { decisionProblems, makeDecision } from "../lib/decisions.js";
import { buildManifest } from "../lib/manifest.js";
import { regionFromDraft } from "../lib/regions.js";
import { createState, recordDecision, setRegions } from "../lib/state.js";

const [findingsPath, outDir, regionResidues] = process.argv.slice(2);
const report = JSON.parse(readFileSync(findingsPath, "utf8"));
const TIME = "2026-01-01T00:00:00Z";

function decide(state, finding, optionId, extra = {}) {
  const draft = { option_id: optionId, rationale: `contract test: ${optionId}`, decided_by: "contract-test", ...extra };
  const problems = decisionProblems(finding, draft);
  if (problems.length) throw new Error(`${finding.id}: ${problems.join(" ")}`);
  return recordDecision(state, makeDecision(finding, draft, TIME));
}

function decideAll(state, findings) {
  return findings.reduce((s, f) => decide(s, f, f.recommended_option), state);
}

async function write(name, state) {
  writeFileSync(join(outDir, `${name}.json`), JSON.stringify(await buildManifest(state, TIME), null, 2));
}

const base = createState(report);
const blocking = report.findings.filter((f) => f.effective_severity === "blocking");
await write("undecided", base);
await write("blocking_decided", decideAll(base, blocking));

const { region, errors } = regionFromDraft({ name: "active_site", description: "contract test region", residues: regionResidues }, []);
if (errors.length) throw new Error(errors.join(" "));
await write("region_all_decided", setRegions(decideAll(base, report.findings), [region]));

const explicit = report.findings.find((f) => f.options.some((o) => o.requires_explicit_choice));
if (explicit) {
  const option = explicit.options.find((o) => o.requires_explicit_choice).id;
  await write("explicit_choice", decide(decideAll(base, blocking), explicit, option, { confirmed_explicit: true }));
}
