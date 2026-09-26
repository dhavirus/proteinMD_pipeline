# ADR-0002: Detection in the browser (Pyodide) vs. front end consuming findings.json

- Status: accepted for v0.1 (TASK-001 spike); revisit when the front end is built
- Date: 2026-09-25

## Context

The front end (static site, GitHub Pages) must show findings for an uploaded structure.
Either it runs the Python detectors in the browser via Pyodide, or it consumes a
`findings.json` produced by the CLI / Colab. Decide by measurement.

## Spike (timeboxed)

Code: `spikes/pyodide/` (README there reproduces every number).
Host: Node 22 running the Pyodide 314.0.7 npm package (CPython 3.14.2 on wasm32,
emscripten 5.0.3). Same Pyodide build as the browser; timings are indicative only (no
browser was available in this environment, and cdn.jsdelivr.net was blocked, so no
package could be loaded from the Pyodide CDN).

### 1. Package availability (from `pyodide-lock.json` of 314.0.7, 357 packages)

| package | in Pyodide distribution | note |
|---|---|---|
| gemmi | **no** | no wasm wheel on PyPI either (0.7.5: none), no npm build found |
| biopython | yes (1.87) | would be an alternative parser, different from the CLI's |
| numpy | yes (2.4.6) | not needed by simprep |
| pyyaml | yes (6.0.3) | not needed by the detector core |
| jsonschema | yes (4.26.0) | not needed in the browser (ajv validates there) |

### 2. Does the detector package run?

Yes. After refactoring, the detector core (structure model, detectors, severity,
findings, report) imports nothing compiled: gemmi, PyYAML and jsonschema are confined
to `structure/parse.py`, `knowledge.py` and `schemas.py` (guarded by
`tests/unit/test_pure_core.py`). The spike parses 5FQL with gemmi in CPython, exports the
parser-independent model as JSON, and runs `run_audit` under Pyodide. **Findings are
identical to the CPython run (16/16, byte-for-byte equal dicts; `check_parity.py`).**

### 3. Measurements, 5FQL (8,576 atom records, 1.27 MB model JSON), 3 runs

| step | time |
|---|---|
| load Pyodide runtime | 1.78-1.83 s |
| copy simprep sources into the FS | 7 ms |
| import simprep | 44-45 ms |
| rebuild model from JSON | 58-64 ms |
| audit (all detectors + severity) | 25-26 ms |
| end to end | 1.95-1.98 s |

Memory: Node RSS 177-185 MB, wasm heap 31 MB. Download weight of the Pyodide core
(`pyodide.asm.wasm` 9.6 MB + `python_stdlib.zip` 2.5 MB + `pyodide.asm.mjs` 1.25 MB) is
about 13 MB before compression. For comparison, CPython: gemmi parse 62 ms, audit 9 ms.

## Options

1. **Full in-browser detection.** Needs a parser in the browser. gemmi is unavailable,
   so this means a second parser (Biopython under Pyodide, a pure-Python mmCIF reader,
   or Mol*'s JS parser feeding the model). Two parsers can disagree (altlocs, entity
   types, insertion codes, annotation categories), which breaks "same file -> same
   findings" between the browser and the CLI unless both are kept under a parity test
   on the whole panel.
2. **Front end consumes findings.json from the CLI / Colab.** One parser (gemmi), one
   code path, provenance recorded where the file was parsed. The front end needs no
   Python; it validates findings.json with ajv and records decisions.
3. **Hybrid.** Parse with gemmi (CLI/Colab); the browser re-runs only the context step
   (severity + recommendation) when the user edits regions, via Pyodide on the exported
   model. Costs about 2 s once per session plus about 13 MB download.

## Decision

Option 2 for v0.1: the front end consumes `findings.json` (and, for region editing, a
manifest) produced by `simprep audit`. Keep the detector core Pyodide-compatible (the
guard test stays) so option 3 remains cheap if interactive re-weighting by region turns
out to matter. Option 1 is rejected until gemmi (or a parity-tested parser) exists for
wasm.

## Consequences

- The front end must be given findings.json; it cannot audit an arbitrary upload on
  its own. Running `simprep audit` (locally or in Colab) is the entry point.
- Changing regions in the UI does not re-escalate severities until the audit is re-run,
  unless option 3 is adopted later; the UI must say so.
- Not measured, to redo when a browser is available: cold load over the network from
  the CDN, browser memory, Safari/Firefox behaviour.
