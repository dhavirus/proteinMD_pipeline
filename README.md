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
pytest                                                   # unit + panel regression tests
ruff check . && ruff format --check .
```

Layout: `schema/` (JSON Schemas, the source of truth), `knowledge/` (versioned YAML
rules), `src/simprep/` (parser, detectors, severity, manifest, CLI), `tests/panel/`
(regression structures and independently derived expected findings).
