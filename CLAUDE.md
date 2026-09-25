# CLAUDE.md — simprep

Working name: `simprep` (rename freely; keep the package name and this file in sync).

## What this repository is

A general platform that turns an experimental or predicted protein structure into a
curated, fully documented simulation system, and later drives comparative free-energy
landscape studies (WT vs variants) on Google Colab.

The core idea: **target-dependent features** (metal sites, non-standard residues such as
formylglycine, glycans, covalent ligands, buried waters, ambiguous site-bound ions, …)
are detected deterministically, explained with evidence, and paired with treatment
options and a recommended default. The human decides; the decision is recorded.

Instance #1, used to validate every feature: iduronate-2-sulfatase (IDS), PDB 5FQL
(Ca²⁺ site, formylglycine, glycans), WT vs R468Q vs R468W.

## Guiding principle

**Generalize the interface and the plumbing. Never generalize the science.**
- Structure audit, include/exclude decisions, mutation building, provenance, run
  bookkeeping: general.
- Which region to sample, which CVs are meaningful, which metal model is adequate for a
  given question: per-study judgment, captured in the manifest, never automated.
- No feature is built unless IDS or a named regression-panel structure exercises it.

## Architecture

```
schema/        JSON Schemas — the single source of truth
  manifest.schema.json   every decision about a system, with rationale + provenance
  rule.schema.json       format of a knowledge-base rule
  finding.schema.json    what a detector emits
knowledge/     versioned YAML rules, one file per rule family (metals.yaml, …)
src/simprep/
  structure/   parsing (gemmi), biological assembly, residue/atom model
  detectors/   one module per rule family; pure functions: structure -> findings
  severity/    context weighting: distance of a finding to manifest regions
  manifest/    build, validate, hash inputs, serialize
  cli.py       `simprep audit <file>` -> findings.json (+ human-readable report)
frontend/      static site (GitHub Pages): upload, Mol* view, findings panel,
               per-finding decisions, region selection, manifest export
modules/       Colab back-end modules (prep, parameterization, …) — later tasks
notebooks/     thin Colab notebooks that call modules/
tests/
  panel/       regression structures (committed, gzipped) + expected_findings/*.yaml
docs/
  decisions/   ADRs (architecture decision records), numbered
  tasks/       task specs handed to Claude Code
```

Layers communicate only through schema-validated JSON (findings.json, manifest.json).
The front end never transforms coordinates; it records decisions. All coordinate work
happens in Python, driven by the manifest.

## Non-negotiable rules

1. **Detectors are deterministic.** Given the same file and knowledge-base version they
   emit the same findings. No LLM calls inside detection or decision logic.
2. **Flag, don't decide.** Every finding carries evidence (atoms, distances, occupancy,
   B-factors, source records), options with trade-offs, a recommended default, and a
   rationale. Nothing is excluded or modified without an explicit recorded decision.
3. **Unknown chemistry is a hard stop.** Anything not matched by a rule becomes an
   `unrecognized` finding with severity `blocking` and "expert decision required".
   Never fall back silently to a default.
4. **Severity depends on the question.** A finding's base severity comes from its rule;
   its effective severity is escalated by proximity to regions defined in the manifest
   (thresholds are configurable, recorded, and never hardcoded).
5. **Nothing is lost silently.** Every input record (chain, residue, HETATM, water,
   altloc) either survives to the output or is excluded by a recorded decision with a
   reason. Report counts in and counts out.
6. **Provenance on everything.** Manifest records: input-file SHA-256, knowledge-base
   version, simprep version / git commit, schema version, timestamps, every decision.
7. **Science claims need sources.** Rule rationales that cite literature must give a real,
   checkable reference. If a claim or number was not verified in this session, mark it
   `[VERIFY]` in the rule file. Never invent a citation, a PDB ID, or a chemical
   component ID.
8. **Rules are data, not code.** Adding support for a new feature type should mean adding
   a YAML rule (plus a detector function only if no existing matcher fits).

## Tech choices (change only via an ADR in docs/decisions/)

- Python ≥ 3.11. Structure parsing: `gemmi` (mmCIF/PDB, struct_conn, assemblies).
- Validation: `jsonschema` in Python; `ajv` in the browser (pinned CDN build).
- Tests: `pytest`; CI on GitHub Actions for every push/PR.
- Front end: plain HTML/JS (no build step unless an ADR justifies one), Mol* loaded from
  a CDN at a pinned exact version, deployable to GitHub Pages.
- Open question tracked as ADR-0002: run Python detectors in the browser via Pyodide, or
  have the front end consume findings.json produced by the CLI/Colab. Decide by spike,
  not by assumption.

## Code standards

- Pure core, thin I/O shell: detectors and severity logic are pure functions
  (structure + config -> findings); file and network access live at the edges.
- Functions ~20 lines, one job, ≤ 3 arguments (bundle into a dataclass beyond that).
- Units in names: `cutoff_angstrom`, `distance_angstrom`, `energy_kcal_per_mol`.
- Type hints and docstrings on the public surface. No magic numbers — constants come
  from config or the rule file, with a citation comment where they are literature values.
- Raise specific exceptions with actionable messages; never swallow errors. In batch
  loops, count and report failures.
- No dead or commented-out code. Match the structure to the scale — no abstraction that
  does not remove more complexity than it adds.

## Testing

- Every rule family has unit tests on minimal hand-built structures AND regression tests
  on the panel.
- **Expected findings must not be derived from detector output** (that is circular).
  Derive them independently from the file's own annotations (struct_conn / LINK,
  modified-residue records, unobserved-residue records, entity and assembly records)
  or from the primary publication, and mark each expected file `reviewed: false` until
  the maintainer (Bem) reviews it.
- CI must pass before a PR is opened. Tests must not require network access: panel
  structures are committed to the repo.

## Colab back-end conventions (for modules/ and notebooks/)

1. Mount Drive; PROJECT root `/content/drive/MyDrive/googleColab_run/<tool_name>/`
   (demos under `/content/drive/MyDrive/googleColab_test/`).
2. Runs in `./<yymm>/run_<yymmdd>_<hhmmss>/`, each with saved `config.json` (containing
   the manifest and the simprep git commit) and a `run_id`.
3. All configuration in one dict at the top of the notebook; nothing hardcoded deeper.
4. Seeds set explicitly from config.
5. I/O-heavy work on local `/content`; bulk-copy finished artifacts to Drive.
6. Long runs checkpoint to Drive and resume after a disconnect without manual edits.
7. Notebooks clone this repo at a pinned commit; reload modules with
   `importlib.reload(module)`. Never use `%load_ext autoreload` / `%autoreload`.
8. `.ipynb` files must keep real cell separation and newlines; validate with
   `nbformat.validate()` before committing.

## Environment limits (Claude Code on the web)

- No GPU and no long-running jobs here: this environment builds and tests software only.
  Do not attempt MD, WSME-L, or FEP runs.
- Network may be restricted. If fetching from RCSB fails, say so and ask; do not work
  around it. Committed panel files are the source for tests.

## Workflow

- One branch per task; small commits with meaningful messages.
- Record every architectural choice as an ADR (context, options, decision, consequences).
- End each task with a PR whose description contains: what was done, what was not done
  and why, open questions for the maintainer, and every `[VERIFY]` item introduced.
- When a spec in docs/tasks/ conflicts with this file, stop and ask rather than guess.
