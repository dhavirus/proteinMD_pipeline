# TASK-001 — Foundations: schemas, rule format, first detectors, regression panel

Read `CLAUDE.md` first. This task builds the foundation everything else depends on.
No front end beyond a feasibility spike, no Colab modules, no simulation.

## Goal

`simprep audit <structure>` produces a schema-valid `findings.json` and a readable
report for any structure in the regression panel, with deterministic detectors for four
rule families, context-weighted severity, and CI that fails on regression.

## Deliverables

### 1. Repository scaffold
Layout as in `CLAUDE.md`. `pyproject.toml`, pinned dependencies, pytest config,
GitHub Actions workflow (lint + tests), `README.md` with a 10-line quickstart.

### 2. Schemas (`schema/`, JSON Schema draft 2020-12, each with `schema_version`)
- `finding.schema.json`: id, rule_id, rule_family, base_severity
  (`info` | `warn` | `blocking`), effective_severity, locus (chain, residue number,
  insertion code, residue name, atom names), evidence (typed key-values: distances in Å,
  occupancies, B-factors, source records), options (id, label, trade-offs, cost hint),
  recommended_option, rationale, references, verify_flags.
- `rule.schema.json`: rule_id, family, version, matcher (declarative criteria the
  detector evaluates), base_severity, options, recommended_option with conditions,
  rationale, references, verify_flags.
- `manifest.schema.json` (v0.1 — enough to support this task, designed to extend):
  input file + SHA-256, knowledge-base version, simprep version/commit, regions of
  interest (named residue sets), severity-escalation thresholds, findings snapshot,
  decisions (finding_id → chosen option, rationale, decided_by, timestamp).
Include valid and invalid example documents for each schema, tested in CI.

### 3. Knowledge base (`knowledge/`)
Rules for these four families only:
- **metals**: detect metal ions; compute coordination shell (ligand atoms within a
  configurable cutoff, count, mean distance, geometry descriptor); classify as likely
  structural / likely catalytic / ambiguous using stated criteria. Options must include
  at least: nonbonded 12-6, nonbonded 12-6-4, bonded (MCPB.py-style), restrained —
  with trade-offs. Recommended option depends on the manifest region context.
- **nonstandard_residues**: residues in polymer chains that are not the 20 standard
  amino acids, detected from both residue names and modification annotations. Options:
  route to topology builder, revert to parent residue (only if chemically defensible and
  explicitly chosen), exclude segment. Include a specific rule for formylglycine that
  flags the aldehyde vs gem-diol hydration state as an explicit decision. Look up the
  chemical component IDs actually used in 5FQL from the file itself; do not assume them.
- **altlocs**: residues/atoms with alternative locations; evidence = occupancies,
  B-factors, spatial spread. Options: highest occupancy, specific altloc, keep for
  ensemble (multiple systems).
- **missing_residues**: unobserved residues from the file's annotations, split into
  termini vs internal gaps, with gap length and flanking-residue distance. Options:
  model loop, cap termini, truncate — internal gaps near regions of interest escalate.
Every rule YAML must validate against `rule.schema.json` in CI.

### 4. Detector framework (`src/simprep/`)
- Pure functions: `(Structure, RuleSet, AuditConfig) -> list[Finding]`, one module per
  family, registered by family name.
- Anything in the structure that no rule claims (non-water HETATM, unusual polymer
  residue, unexplained covalent link) → `unrecognized` finding, `blocking`.
- Severity module: escalate effective severity by minimum heavy-atom distance to any
  manifest region, thresholds from config. With no regions defined, effective = base.
- Record counts in / counts out per record type in the audit report.
- CLI: `simprep audit FILE [--manifest M] [--out DIR]` → `findings.json` + `report.md`.

### 5. Regression panel (`tests/panel/`)
Five structures, committed gzipped (mmCIF preferred):
- **5FQL** (required).
- A zinc metalloenzyme with a catalytic Zn site.
- A glycoprotein with N-linked glycans resolved.
- A structure with a covalently bound ligand.
- A structure with altlocs and at least one internal missing loop.
For each non-fixed entry: propose 2–3 candidates, verify each against RCSB metadata
(if network access allows), choose one, and document the choice and reasoning in
`tests/panel/PANEL.md`. If you cannot verify an ID, say so — do not guess.
Glycans and covalent ligands have no rule family yet: their findings must appear as
`unrecognized` / `blocking`. That is the expected, tested behavior for this task.

For each structure write `expected_findings/<id>.yaml`, derived independently of the
detectors (file annotations or the primary publication), with `reviewed: false`.
Regression tests compare detector output against these files; mismatches fail CI.

### 6. ADRs (`docs/decisions/`)
- ADR-0001: schema-first architecture and the three-layer contract.
- ADR-0002: Pyodide spike. Test whether gemmi (or an alternative parser) and the
  detector package run in the browser via Pyodide: package availability, load time,
  memory with 5FQL. Record measurements and a recommendation (in-browser detection vs
  front end consuming findings.json). Timebox this; a clear "no, because…" is a valid
  result.

## Acceptance criteria
- `pytest` passes locally and in CI; CI runs without network.
- `simprep audit` on each panel structure produces schema-valid findings.json.
- 5FQL audit surfaces: the Ca²⁺ site with coordination evidence, formylglycine with the
  hydration-state decision, glycans as `unrecognized`, plus any altlocs and missing
  residues actually present.
- Severity escalation demonstrated by a test: the same metal finding is lower severity
  with a distant region than with a region containing it.
- Every rule and schema validated in CI; every `[VERIFY]` item listed in the PR.

## Out of scope
Front-end UI (beyond the ADR-0002 spike), coordinate modification, protonation,
parameterization, Colab modules, rule families other than the four above.

## PR description must include
What was done; what was not and why; the panel choices with verification status; the
ADR-0002 recommendation; all `[VERIFY]` items; open questions for the maintainer.
