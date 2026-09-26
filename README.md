# simprep

Turns an experimental or predicted protein structure into a curated, fully documented
simulation system. Target-dependent features (metal sites, non-standard residues,
altlocs, missing residues, unknown chemistry) are detected deterministically, explained
with evidence, and paired with treatment options and a recommended default. The human
decides; the manifest records the decision. See `CLAUDE.md` and `docs/decisions/`.

## Quickstart

```bash
python -m pip install -e '.[dev]'                        # Python >= 3.11, from a checkout
simprep audit tests/panel/5FQL.cif.gz --out out/5fql     # -> findings.json + report.md
simprep manifest init tests/panel/5FQL.cif.gz --out m.json
# edit m.json: add regions of interest (named residue sets)
simprep audit tests/panel/5FQL.cif.gz --manifest m.json --out out/5fql_roi
simprep manifest status out/5fql_roi/manifest.json        # exit 3 while blocking findings are undecided
pytest                                                   # unit + panel regression + contract tests
(cd frontend && node --test "tests/*.test.mjs")           # front-end logic (Node 22, no dependencies)
ruff check . && ruff format --check .
```

## Review page

`frontend/` is a static page (published to GitHub Pages from `main`) for reviewing an
audit: open the structure and its `findings.json` (and optionally a `manifest.json`),
inspect each finding in Mol*, record decisions, edit regions of interest, and export
`manifest.json`. Files are read in the browser and never uploaded. Locally:

```bash
python -m http.server 8000        # from the repository root, then open http://localhost:8000/frontend/
python frontend/demo/build_demo.py  # regenerate the 5FQL example after rule changes (CI checks it)
node frontend/tests/smoke/smoke.mjs  # browser smoke test (needs Playwright and network)
```

Layout: `schema/` (JSON Schemas, the source of truth), `knowledge/` (versioned YAML
rules), `src/simprep/` (parser, detectors, severity, manifest, CLI), `frontend/` (review
page), `tests/panel/` (regression structures and independently derived expected findings).
