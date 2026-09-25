# ADR-0001: Schema-first architecture and the three-layer contract

- Status: accepted (TASK-001)
- Date: 2026-09-25

## Context

simprep turns a deposited or predicted structure into a documented simulation system.
Three kinds of component take part: deterministic Python detectors, a human-facing front
end (static site), and Colab back-end modules that transform coordinates. They are
written at different times, in different languages, and must stay auditable: every
decision needs its evidence, its rationale and its provenance.

## Options

1. **Shared Python objects.** The front end and back end import the same package.
   Rejected: the front end is plain JS, and it couples every layer to one runtime.
2. **Ad-hoc JSON, documented in prose.** Cheap at first; drifts silently, because
   nothing checks it.
3. **Schema-first: JSON Schema as the single source of truth.** Every document that
   crosses a layer boundary has a versioned schema, is validated on write and on read,
   and ships with valid and invalid examples tested in CI.

## Decision

Option 3. Documents that cross a boundary are validated JSON (JSON Schema draft 2020-12):

| document | schema | producer | consumers |
|---|---|---|---|
| rule files (`knowledge/*.yaml`) | `rule.schema.json` | maintainers | detectors |
| finding | `finding.schema.json` | detectors | report, front end, manifest |
| `findings.json` | `findings_report.schema.json` | `simprep audit` | front end, manifest |
| `manifest.json` | `manifest.schema.json` | front end / CLI | Colab modules |

The three layers:

1. **Detection (Python, pure core).** `Structure` model -> findings. Gemmi is used only
   at the edge (`structure/parse.py`); detectors, severity and recommendation logic
   import neither gemmi, PyYAML nor jsonschema (so they also run under Pyodide, see
   ADR-0002).
2. **Decision (front end).** Displays findings and records a human decision per finding
   into the manifest. It never transforms coordinates.
3. **Execution (Colab modules, later tasks).** Apply the manifest's decisions to
   coordinates; every input record either survives or is excluded by a recorded
   decision.

Specific choices:

- Schemas use `$id` URNs (`urn:simprep:schema:<name>`) and cross-reference each other;
  validators load all schemas into one registry (Python `referencing`; `ajv` in the
  browser via `addSchema`).
- Each schema declares `x-schema-version`; documents carry `schema_version` (const).
  A test asserts the two agree. Any breaking change bumps the version.
- Rule files may carry `x-` prefixed top-level keys (YAML anchors shared between rules).
- Finding IDs derive from family + locus (`metals/A:1601`), not from rule IDs, so a
  decision keyed by finding ID survives a rule refactor that keeps the locus.
- A rule's recommendation is data: an ordered list of conditions on region proximity
  and evidence values, falling back to a default. Options flagged
  `requires_explicit_choice` can never be recommended (checked at load time).
- `findings.json` records the audit configuration actually used (regions, thresholds,
  manifest hash) and counts in / passed through / flagged / excluded per record type.

## Consequences

- Adding a field is a schema change with a version bump and example updates; this is
  deliberate friction.
- The front end can validate everything it reads without Python.
- A fourth schema (`findings_report`) was added beyond the three named in TASK-001,
  because the audit output wraps findings with provenance and accounting.
- `knowledge/` and `schema/` live beside `src/`, so simprep runs from a checkout
  (editable install or pinned clone), not from a wheel. Revisit if packaging is needed.
