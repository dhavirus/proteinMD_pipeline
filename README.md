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
# decide every finding (review page below), then apply the decisions:
simprep prep tests/panel/5FQL.cif.gz --manifest decided.json --out out/5fql_prep
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

## Prep

`simprep prep` applies a manifest in which every blocking finding has a final decision
(`expert_review` / `add_rule` are not final). It selects altlocs, deletes excluded
residues, writes `add_link` bonds into `struct_conn`, and never invents coordinates:
loop modelling, capping, chemistry edits and parameterization go into the work order in
`prep_record.json` for later stages. Outputs: `system.cif` + `system.pdb` (or one pair
per altloc with `keep_ensemble`), `prep_record.json` (actions, work order, file hashes,
record accounting) and `prep_report.md`. See ADR-0004.

On Colab, `notebooks/01_prep.ipynb` wraps the same function with the run layout from
`CLAUDE.md` (`modules/prep_module.py`).

## Variants

List variants in the manifest (`"variants"`, or the review page's variant editor, e.g.
`A:468 R>Q`), then run `simprep variants STRUCTURE --manifest M --out DIR`. The first run
writes `DIR/manifest.json` with every side-chain rotamer candidate as a `variant_build`
finding (exit 3); decide them on the review page, then run again with that manifest to
get `wt/`, one directory per variant, `variant_record.json` and `variant_report.md`.
Side chains are built from cited rotamer data in `knowledge/side_chains.yaml`; a rotamer
that clashes with the kept atoms is never built. See ADR-0005 and
`notebooks/02_variants.ipynb`.

Layout: `schema/` (JSON Schemas, the source of truth), `knowledge/` (versioned YAML
rules and residue mappings), `src/simprep/` (parser, detectors, severity, manifest, prep,
CLI), `frontend/` (review page), `modules/` + `notebooks/` (Colab back end),
`tests/panel/` (regression structures, independently derived expected findings, prep
fixtures).
