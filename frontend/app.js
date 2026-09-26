// simprep review page: wires the pure logic in lib/ to the DOM, Mol* and ajv.
// All state lives in one immutable object (lib/state.js); every change re-renders.
import { decisionProblems, isFinal, isoSeconds, makeDecision, OPTION_PARAMETERS } from "./lib/decisions.js";
import {
  altlocChoices, FAMILY_LABELS, formatAtom, formatEvidenceValue, groupFindings, SEVERITIES, severityCounts,
} from "./lib/findings.js";
import { buildManifest, manifestMismatch } from "./lib/manifest.js";
import { regionFromDraft } from "./lib/regions.js";
import { formatResidueList } from "./lib/residues.js";
import {
  applyDraft, createState, decidedIds, dropOrphans, recordDecision, removeDecision, select, setFilter, setRegions,
  severityIsStale,
} from "./lib/state.js";
import { clearDraft, loadDraft, saveDraft } from "./ui/draft.js";
import { download, fetchBytes, parseJson, readBytes, readStructure } from "./ui/io.js";
import { createValidator } from "./ui/validate.js";
import { createView } from "./ui/view3d.js";

const SCHEMA_BASE = new URL("../schema/", import.meta.url);
const DEMO = {
  structure: new URL("../tests/panel/5FQL.cif.gz", import.meta.url),
  findings: new URL("demo/5FQL.findings.json", import.meta.url),
};

const $ = (id) => document.getElementById(id);
let state = null;
let structure = null;
let view = null;
let validatorPromise = null;
let editingRegion = null;

function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs)) {
    if (value === null || value === undefined || value === false) continue;
    if (key === "class") node.className = value;
    else if (key === "text") node.textContent = value;
    else if (key.startsWith("on")) node.addEventListener(key.slice(2), value);
    else node.setAttribute(key, value === true ? "" : value);
  }
  for (const child of children.flat()) if (child !== null && child !== undefined && child !== false) node.append(child);
  return node;
}

const validator = () => (validatorPromise ??= createValidator(SCHEMA_BASE));

// ---------------------------------------------------------------- opening files

async function openDocuments({ structureName, structureBytes, findingsBytes, manifestBytes }) {
  const validate = await validator();
  const report = parseJson(findingsBytes, "findings.json");
  const reportErrors = validate("findings_report", report);
  if (reportErrors.length) throw new DocumentError("findings.json does not match the findings schema:", reportErrors);
  const opened = await readStructure(structureName, structureBytes);
  if (opened.sha256 !== report.input.sha256) {
    throw new DocumentError("This structure is not the file findings.json was computed for:", [
      `structure ${structureName}: SHA-256 ${opened.sha256}`,
      `findings.json input (${report.input.path}): SHA-256 ${report.input.sha256}`,
    ]);
  }
  const manifest = manifestBytes ? checkedManifest(parseJson(manifestBytes, "manifest.json"), report, validate) : null;
  return { report, manifest, opened };
}

function checkedManifest(manifest, report, validate) {
  const errors = validate("manifest", manifest);
  if (errors.length) throw new DocumentError("manifest.json does not match the manifest schema:", errors);
  const mismatch = manifestMismatch(manifest, report);
  if (mismatch.length) throw new DocumentError("This manifest belongs to another structure:", mismatch);
  return manifest;
}

class DocumentError extends Error {
  constructor(message, details) {
    super(message);
    this.details = details;
  }
}

async function start(documents) {
  structure = documents.opened;
  // A loaded manifest is the record: an older browser draft is offered, never applied
  // over it (and not overwritten until the user chooses). Without a manifest the draft
  // is the only copy of earlier work, so it is restored.
  const opened = { ...createState(documents.report, documents.manifest), restoredDraft: null, pendingDraft: null };
  const draft = loadDraft(documents.report.input.sha256);
  if (!draft) state = opened;
  else if (documents.manifest) state = { ...opened, pendingDraft: draft };
  else state = { ...applyDraft(opened, draft), restoredDraft: draft.savedAt };
  $("open-panel").hidden = true;
  $("review").hidden = false;
  $("top-actions").hidden = false;
  $("subtitle").textContent =
    `${documents.report.structure_summary.name}: ${documents.report.input.path}. Knowledge base ${documents.report.knowledge_base.version}, simprep ${documents.report.simprep.version}.`;
  render();
  await startView();
}

async function startView() {
  try {
    view ??= await createView($("structure-view"));
    view.setStructure(structure.text, structure.format);
    await refreshView();
  } catch (error) {
    $("viewer-status").textContent = error.message;
  }
}

async function refreshView() {
  if (!view) return;
  const finding = currentFinding();
  $("viewer-status").textContent = "Drawing…";
  try {
    await view.show(finding, state.regions);
    $("viewer-status").textContent = finding ? `Showing ${finding.id}` : "Whole structure";
  } catch (error) {
    $("viewer-status").textContent = `3D view failed: ${error.message}`;
  }
}

function showOpenErrors(error) {
  const details = error.details || [];
  $("open-errors").replaceChildren(el("li", { text: error.message }), ...details.map((d) => el("li", { class: "mono", text: d })));
  $("open-status").textContent = "";
}

async function onOpenSubmit(event) {
  event.preventDefault();
  $("open-errors").replaceChildren();
  $("open-status").textContent = "Opening…";
  try {
    const structureFile = $("file-structure").files[0];
    const manifestFile = $("file-manifest").files[0];
    await start(await openDocuments({
      structureName: structureFile.name,
      structureBytes: await readBytes(structureFile),
      findingsBytes: await readBytes($("file-findings").files[0]),
      manifestBytes: manifestFile ? await readBytes(manifestFile) : null,
    }));
  } catch (error) {
    showOpenErrors(error);
  }
}

async function onLoadDemo() {
  $("open-errors").replaceChildren();
  $("open-status").textContent = "Loading the 5FQL example…";
  try {
    await start(await openDocuments({
      structureName: "5FQL.cif.gz",
      structureBytes: await fetchBytes(DEMO.structure),
      findingsBytes: await fetchBytes(DEMO.findings),
      manifestBytes: null,
    }));
  } catch (error) {
    showOpenErrors(error);
  }
}

// ---------------------------------------------------------------- rendering

const currentFinding = () => state.report.findings.find((f) => f.id === state.selectedId) || null;

function update(next, { redraw = false } = {}) {
  state = next;
  if (!state.pendingDraft) saveDraft(state.report.input.sha256, state, isoSeconds(new Date()));
  render();
  if (redraw) refreshView();
}

function render() {
  renderSummary();
  renderNotices();
  renderFilters();
  renderList();
  renderDetail();
  renderRegions();
  renderExport();
}

function renderSummary() {
  const findings = state.report.findings;
  const counts = severityCounts(findings);
  const final = decidedIds(state);
  const decided = final.size;
  const blockingOpen = findings.filter((f) => f.effective_severity === "blocking" && !final.has(f.id)).length;
  $("summary").replaceChildren(
    el("span", { class: "chip neutral", text: `${findings.length} findings` }),
    ...SEVERITIES.map((s) => el("span", { class: `chip ${s}`, text: `${counts[s]} ${s}` })),
    el("span", { class: "chip decided", text: `${decided} of ${findings.length} decided` }),
    el("span", { class: blockingOpen ? "chip blocking" : "chip decided", text: blockingOpen ? `${blockingOpen} blocking undecided` : "all blocking decided" }),
  );
}

function renderNotices() {
  const notices = [];
  if (state.restoredDraft) {
    notices.push(el("div", { class: "notice" },
      `Restored unsaved work from ${state.restoredDraft} (kept in this browser only; the exported manifest is the record). `,
      el("button", { class: "btn", type: "button", text: "Discard it", onclick: onDiscardDraft })));
  }
  if (state.pendingDraft) {
    const draft = state.pendingDraft;
    notices.push(el("div", { class: "notice stale" },
      `This browser also has unsaved work for this structure from ${draft.savedAt} (${Object.keys(draft.decisions).length} decisions). `
      + "The opened manifest is shown. Restoring replaces its decisions and regions with the unsaved work. ",
      el("button", { class: "btn", type: "button", text: "Restore it",
        onclick: () => update({ ...applyDraft(state, draft), pendingDraft: null, restoredDraft: draft.savedAt }, { redraw: true }) }), " ",
      el("button", { class: "btn", type: "button", text: "Discard it", onclick: onDiscardDraft })));
  }
  if (severityIsStale(state)) {
    notices.push(el("div", { class: "notice stale" },
      "Regions changed: the severities shown were computed for the previous regions. Export the manifest and run ",
      el("span", { class: "mono", text: "simprep audit STRUCTURE --manifest manifest.json" }), " to update them."));
  }
  if (state.orphans.length) {
    notices.push(el("div", { class: "notice orphans" },
      `${state.orphans.length} decision(s) in the manifest refer to findings that are not in this audit and will not be exported: `,
      el("span", { class: "mono", text: state.orphans.map((d) => `${d.finding_id} (${d.option_id})`).join(", ") }), " ",
      el("button", { class: "btn", type: "button", text: "Acknowledge and drop them", onclick: () => update(dropOrphans(state)) })));
  }
  $("notices").replaceChildren(...notices);
}

function renderFilters() {
  const severity = el("select", { id: "filter-severity", "aria-label": "Severity filter",
    onchange: (e) => update(setFilter(state, { severities: e.target.value ? [e.target.value] : null })) },
  el("option", { value: "", text: "All severities" }),
  ...SEVERITIES.map((s) => el("option", { value: s, text: s, selected: state.filter.severities?.[0] === s })));
  const family = el("select", { id: "filter-family", "aria-label": "Rule family filter",
    onchange: (e) => update(setFilter(state, { family: e.target.value || null })) },
  el("option", { value: "", text: "All families" }),
  ...Object.entries(FAMILY_LABELS).map(([k, v]) => el("option", { value: k, text: v, selected: state.filter.family === k })));
  const undecided = el("input", { type: "checkbox", id: "filter-undecided", checked: state.filter.undecidedOnly,
    onchange: (e) => update(setFilter(state, { undecidedOnly: e.target.checked })) });
  $("filters").replaceChildren(severity, family, el("label", { for: "filter-undecided" }, undecided, " Undecided only"));
}

function renderList() {
  const groups = groupFindings(state.report.findings, { ...state.filter, decidedIds: decidedIds(state) });
  const nodes = groups.flatMap((group) => [
    el("h3", { text: `${group.label} · ${group.items.length}` }),
    ...group.items.map(findingRow),
  ]);
  $("list").replaceChildren(...(nodes.length ? nodes : [el("p", { class: "msg empty", text: "No findings match these filters." })]));
}

function findingRow(finding) {
  const escalated = finding.base_severity !== finding.effective_severity;
  return el("button", { class: "row", type: "button", "aria-current": String(finding.id === state.selectedId),
    onclick: () => update(select(state, finding.id), { redraw: true }) },
  el("span", { class: `stripe ${finding.effective_severity}` }),
  el("span", {}, el("b", { text: finding.title }), el("span", { class: "mono id", text: finding.id })),
  el("span", { class: "tags" },
    escalated ? el("span", { class: "escalated", text: `${finding.base_severity} →` }) : null,
    el("span", { class: `chip ${finding.effective_severity}`, text: finding.effective_severity }),
    decisionChip(finding)));
}

function decisionChip(finding) {
  const decision = state.decisions[finding.id];
  if (!decision) return null;
  return isFinal(finding, decision.option_id)
    ? el("span", { class: "chip decided", text: "decided" })
    : el("span", { class: "chip warn", text: "not final" });
}

function renderDetail() {
  const finding = currentFinding();
  if (!finding) {
    $("detail").replaceChildren(el("p", { class: "msg", text: "Select a finding." }));
    return;
  }
  $("detail").replaceChildren(
    el("div", {}, el("h2", { text: finding.title }), detailMeta(finding)),
    evidenceTable(finding),
    decisionForm(finding),
    section("Recommendation basis", el("p", { class: "msg", text: finding.recommendation_basis })),
    section("Rule rationale", el("p", { class: "msg", text: finding.rationale })),
    finding.references.length ? section("References", el("ol", { class: "refs" }, ...finding.references.map(reference))) : null,
    finding.verify_flags.length ? section("To verify", el("ul", { class: "verify" }, ...finding.verify_flags.map((v) => el("li", { text: v })))) : null,
  );
}

const section = (title, body) => el("div", {}, el("p", { class: "label", text: title }), body);

function detailMeta(finding) {
  const context = finding.severity_context;
  const where = context.proximity === "no_regions"
    ? "no region set"
    : `${context.proximity}${context.min_distance_angstrom !== null ? `, ${context.min_distance_angstrom} Å from ${context.nearest_region}` : ""}`;
  return el("div", { class: "meta" },
    el("span", { class: `chip ${finding.effective_severity}`, text: finding.effective_severity }),
    el("span", { class: "mono", text: finding.id }),
    el("span", { class: "mono", text: finding.rule_id }),
    el("span", { text: finding.base_severity === finding.effective_severity ? `base ${finding.base_severity}` : `raised from ${finding.base_severity}` }),
    el("span", { text: where }));
}

function evidenceTable(finding) {
  const rows = finding.evidence.map((item) => el("tr", {},
    el("td", { text: item.key }),
    el("td", { class: "v" }, formatEvidenceValue(item),
      item.atoms?.length ? el("div", { class: "mono note", text: item.atoms.map(formatAtom).join(" – ") }) : null,
      item.note ? el("div", { class: "note", text: item.note }) : null)));
  return section(`Evidence · ${finding.evidence.length}`, el("div", { class: "ev-wrap" }, el("table", { class: "ev" }, el("tbody", {}, rows))));
}

function reference(ref) {
  return el("li", {}, ref.citation, ref.doi ? " " : null,
    ref.doi ? el("a", { href: `https://doi.org/${ref.doi}`, target: "_blank", rel: "noopener", text: `doi:${ref.doi}` }) : null);
}

function decisionForm(finding) {
  const saved = state.decisions[finding.id];
  const form = el("form", { class: "form", id: "decision-form", onsubmit: (e) => onSaveDecision(e, finding) },
    el("p", { class: "label", text: saved ? `Decision · recorded ${saved.timestamp}` : "Decision" }),
    el("fieldset", { class: "options" }, el("legend", { class: "label", text: "Treatment" }),
      ...finding.options.map((option) => optionChoice(finding, option, saved))),
    el("div", { id: "decision-extra" }),
    el("label", { class: "field", for: "decision-rationale", text: "Rationale" }),
    el("textarea", { id: "decision-rationale", placeholder: "Why this option, for this study" }),
    el("label", { class: "field", for: "decision-by", text: "Decided by" }),
    el("input", { type: "text", id: "decision-by", autocomplete: "name" }),
    el("div", { class: "form-row" },
      el("button", { class: "btn primary", type: "submit", text: saved ? "Update decision" : "Record decision" }),
      saved ? el("button", { class: "btn", type: "button", text: "Remove decision", onclick: () => update(removeDecision(state, finding.id)) }) : null,
      el("span", { class: "status", id: "decision-status", role: "status" })));
  form.querySelector("#decision-rationale").value = saved?.rationale || "";
  form.querySelector("#decision-by").value = saved?.decided_by || lastDecider();
  form.addEventListener("change", (e) => { if (e.target.name === "option") renderDecisionExtra(form, finding, saved); });
  renderDecisionExtra(form, finding, saved);
  return form;
}

function optionChoice(finding, option, saved) {
  const id = `option-${option.id}`;
  const checked = (saved ? saved.option_id : finding.recommended_option) === option.id;
  return el("label", { class: "opt", for: id },
    el("input", { type: "radio", name: "option", id, value: option.id, checked, "data-final": String(isFinal(finding, option.id)) }),
    el("span", { class: "name" }, option.label, el("span", { class: "tag", text: `cost ${option.cost_hint}` }),
      option.id === finding.recommended_option ? el("span", { class: "tag rec", text: "recommended" }) : null,
      option.requires_explicit_choice ? el("span", { class: "tag explicit", text: "explicit choice only" }) : null),
    el("ul", {}, ...option.trade_offs.map((t) => el("li", { text: t }))));
}

/** Extra inputs for the chosen option: explicit-choice confirmation, altloc parameter. */
function renderDecisionExtra(form, finding, saved) {
  const optionId = form.querySelector('input[name="option"]:checked')?.value;
  const option = finding.options.find((o) => o.id === optionId);
  const extra = form.querySelector("#decision-extra");
  const nodes = [];
  if (option?.requires_explicit_choice) {
    nodes.push(el("label", { class: "confirm", for: "decision-confirm" },
      el("input", { type: "checkbox", id: "decision-confirm", checked: saved?.option_id === optionId }),
      el("span", { text: `I choose "${option.label}" deliberately. It is never applied by default: ${option.trade_offs.join(" ")}` })));
  }
  if (option && !isFinal(finding, optionId)) {
    nodes.push(el("p", { class: "msg", text: `"${option.label}" records that the decision is still open: the finding stays undecided, and prep will not run until it gets a final decision.` }));
  }
  if (OPTION_PARAMETERS[optionId]) nodes.push(...altlocInput(finding, saved));
  extra.replaceChildren(...nodes);
}

function altlocInput(finding, saved) {
  const choices = altlocChoices(finding);
  const current = saved?.parameters?.altloc || "";
  const input = choices.length
    ? el("select", { id: "decision-altloc" }, el("option", { value: "", text: "Choose…" }),
      ...choices.map((c) => el("option", { value: c, text: c, selected: c === current })))
    : el("input", { type: "text", id: "decision-altloc", maxlength: "1", value: current, placeholder: "e.g. A" });
  return [el("label", { class: "field", for: "decision-altloc", text: "Altloc to keep" }), input];
}

const lastDecider = () => Object.values(state.decisions).at(-1)?.decided_by || "";

function onSaveDecision(event, finding) {
  event.preventDefault();
  const form = event.currentTarget;
  const draft = {
    option_id: form.querySelector('input[name="option"]:checked')?.value,
    rationale: form.querySelector("#decision-rationale").value,
    decided_by: form.querySelector("#decision-by").value,
    confirmed_explicit: form.querySelector("#decision-confirm")?.checked || false,
    parameters: { altloc: form.querySelector("#decision-altloc")?.value || "" },
  };
  const problems = decisionProblems(finding, draft);
  if (problems.length) {
    const status = form.querySelector("#decision-status");
    status.className = "status error";
    status.textContent = problems.join(" ");
    return;
  }
  update(recordDecision(state, makeDecision(finding, draft, isoSeconds(new Date()))));
  const status = $("decision-status");
  status.className = "status ok";
  status.textContent = "Decision recorded";
}

function onDiscardDraft() {
  clearDraft(state.report.input.sha256);
  const reopened = createState(state.report, state.baseManifest);
  state = { ...reopened, selectedId: state.selectedId, restoredDraft: null, pendingDraft: null };
  render();
  refreshView();
}

// ---------------------------------------------------------------- regions

function renderRegions() {
  const others = state.regions.filter((r) => r.name !== editingRegion?.name).map((r) => r.name);
  const draft = editingRegion ? { ...editingRegion, residues: formatResidueList(editingRegion.residues) } : {};
  const form = el("form", { class: "form", id: "region-form", onsubmit: (e) => onSaveRegion(e, others) },
    el("label", { class: "field", for: "region-name", text: "Name" }),
    el("input", { type: "text", id: "region-name", value: draft.name || "", placeholder: "active_site" }),
    el("label", { class: "field", for: "region-description", text: "Description" }),
    el("input", { type: "text", id: "region-description", value: draft.description || "", placeholder: "Ca2+ site and catalytic FGly" }),
    el("label", { class: "field", for: "region-residues", text: "Residues (chain:number, ranges allowed)" }),
    el("textarea", { id: "region-residues", placeholder: "A:45, A:46, A:84, A:334-335, A:1551" }),
    el("div", { class: "form-row" },
      el("button", { class: "btn primary", type: "submit", text: editingRegion ? "Save region" : "Add region" }),
      editingRegion ? el("button", { class: "btn", type: "button", text: "Cancel", onclick: () => { editingRegion = null; render(); } }) : null),
    el("ul", { class: "errors", id: "region-errors" }));
  form.querySelector("#region-residues").value = draft.residues || "";
  $("regions").replaceChildren(
    el("h2", { id: "regions-title", text: "Regions of interest" }),
    el("p", { class: "msg", text: "Findings near a region get a higher severity when the audit is re-run with the exported manifest." }),
    ...state.regions.map(regionItem),
    form);
}

function regionItem(region) {
  return el("div", { class: "region" },
    el("div", { class: "form-row" }, el("b", { class: "mono", text: region.name }),
      el("span", { class: "controls" },
        el("button", { class: "btn", type: "button", text: "Edit", onclick: () => { editingRegion = region; render(); } }),
        el("button", { class: "btn danger", type: "button", text: "Delete",
          onclick: () => update(setRegions(state, state.regions.filter((r) => r !== region)), { redraw: true }) }))),
    el("span", { class: "msg", text: region.description }),
    el("span", { class: "mono note", text: formatResidueList(region.residues) }));
}

function onSaveRegion(event, otherNames) {
  event.preventDefault();
  const form = event.currentTarget;
  const { region, errors } = regionFromDraft({
    name: form.querySelector("#region-name").value,
    description: form.querySelector("#region-description").value,
    residues: form.querySelector("#region-residues").value,
  }, otherNames);
  if (errors.length) {
    form.querySelector("#region-errors").replaceChildren(...errors.map((e) => el("li", { text: e })));
    return;
  }
  const kept = state.regions.filter((r) => r !== editingRegion);
  editingRegion = null;
  update(setRegions(state, [...kept, region]), { redraw: true });
}

// ---------------------------------------------------------------- export

function renderExport() {
  $("export").replaceChildren(
    el("h2", { id: "export-title", text: "Export manifest" }),
    el("p", { class: "msg" }, "Builds ", el("span", { class: "mono", text: "manifest.json" }),
      " with the regions, the decisions and a snapshot of these findings, checked against the manifest schema. Then run ",
      el("span", { class: "mono", text: "simprep manifest status manifest.json" }), " to confirm every blocking finding is decided."),
    el("div", { class: "form-row" },
      el("button", { class: "btn primary", type: "button", id: "export-build", text: "Build and download", onclick: onExport }),
      el("button", { class: "btn", type: "button", id: "export-copy", text: "Copy JSON", onclick: onCopy }),
      el("span", { class: "status", id: "export-status", role: "status" })),
    el("ul", { class: "errors", id: "export-errors" }),
    el("pre", { class: "json", id: "export-preview", hidden: true }));
}

async function validManifest() {
  const manifest = await buildManifest(state, isoSeconds(new Date()));
  const errors = (await validator())("manifest", manifest);
  $("export-errors").replaceChildren(...errors.map((e) => el("li", { class: "mono", text: e })));
  if (errors.length) {
    $("export-status").className = "status error";
    $("export-status").textContent = "The manifest does not pass the schema; nothing was exported.";
    return null;
  }
  $("export-preview").hidden = false;
  $("export-preview").textContent = JSON.stringify(manifest, null, 2);
  return manifest;
}

async function onExport() {
  const manifest = await validManifest();
  if (!manifest) return;
  download(manifest, "manifest.json");
  $("export-status").className = "status ok";
  $("export-status").textContent = `Downloaded manifest.json (${manifest.decisions.length} decisions, ${manifest.regions.length} regions).`;
}

async function onCopy() {
  const manifest = await validManifest();
  if (!manifest) return;
  try {
    await navigator.clipboard.writeText(`${JSON.stringify(manifest, null, 2)}\n`);
    $("export-status").className = "status ok";
    $("export-status").textContent = "Copied.";
  } catch {
    $("export-status").className = "status";
    $("export-status").textContent = "Copy was blocked: select the JSON below and copy it.";
  }
}

// ---------------------------------------------------------------- boot

$("open-form").addEventListener("submit", onOpenSubmit);
$("load-demo").addEventListener("click", onLoadDemo);
$("open-other").addEventListener("click", () => {
  $("open-panel").hidden = false;
  $("review").hidden = true;
  $("top-actions").hidden = true;
});
validator().catch((error) => showOpenErrors(error));
