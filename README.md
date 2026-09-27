# simprep

Turns an experimental or predicted protein structure into a curated, fully documented
simulation system. Target-dependent features (metal sites, non-standard residues,
altlocs, missing residues, unknown chemistry) are detected deterministically, explained
with evidence, and paired with treatment options and a recommended default. The human
decides; the manifest records the decision. See `CLAUDE.md` and `docs/decisions/`.

## Quickstart

```bash
python -m pip install -e '.[dev,relax]'                  # Python >= 3.11, from a checkout
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
residues, writes `add_link` bonds into `struct_conn`, applies residue mappings that only
delete and rename atoms (`revert_to_parent`; `model_gem_diol`: formylglycine sulfate ALS ->
gem-diol DDZ, ADR-0008), and never invents coordinates:
loop modelling, capping and parameterization go into the work order in
`prep_record.json` for later stages. Outputs: `system.cif` + `system.pdb` (or one pair
per altloc with `keep_ensemble`), `prep_record.json` (actions, work order, file hashes,
record accounting) and `prep_report.md`. See ADR-0004.

On Colab, `notebooks/01_prep.ipynb` wraps the same function with the run layout from
`CLAUDE.md` (`modules/prep_module.py`).

## Modelling

Internal gaps decided `model_loop` are built by `simprep model STRUCTURE --manifest M
--out DIR` (wild type alone) and in-process by `simprep variants`, so every variant
carries the same loop: prep (select) -> model (add) -> variants (build) -> relax.
PDBFixer 1.12.0 places the residues (seeded, Reference platform), OpenMM minimizes the
loop and its restrained flanks (trans-peptide and L-chirality restraints, sterics only),
and the peptide geometry and contacts are checked; a failing check stops the stage
(`modelling.junction_geometry`, exit 3). Outputs: `wt/` (the prepared wild type),
`wt_modelled/` (modelled atoms at occupancy 0.00), `model_record.json` (every modelled
residue and atom, its source, geometry, accounting, work order) and `model_report.md`.
Protocol defaults: `knowledge/modelling.yaml` (override in the manifest's `modelling`).
`truncate` means charged termini, listed for the topology stage. A build takes about two
minutes. See ADR-0007.

## Protonation

State the pH in the manifest (`"protonation": {"ph": 7.2, "ph_rationale": "..."}`; simprep
never defaults it), then run `simprep protonate STRUCTURE --manifest M --out DIR` for the
wild type, or `simprep variants` to protonate every system it writes. PROPKA 3.5.1
estimates pKa values per system; disulfide cysteines and metal ligands (from
`struct_conn`) follow fixed rules; residues near chemistry PROPKA ignores, or within one
pH unit of their pKa, become `protonation` findings with a recorded default; OpenMM adds
every hydrogen (deposited ones are replaced). Outputs: `*_protonated/` systems,
`protonation_record.json` (pH, method, per-residue state and why, differences from the
wild type) and `protonation_report.md`. See ADR-0009.

## Parameterization

`simprep parameterize DIR --manifest M` takes a directory written by `simprep protonate`
or `simprep variants` (with a pH) and builds an OpenMM System per protonated system:
Amber ff14SB + TIP3P as shipped with OpenMM, DDZ from `knowledge/forcefield/ddz.xml`, and
the ion models the manifest decided (`nonbonded_12_6_4`: the Li-Merz 12-6-4 Ca2+ model,
`knowledge/forcefield/lj1264.yaml`). It stops on any residue, bond or angle without
parameters and writes `<name>_parameterized/system.xml` and `parameterization_record.json`.
The DDZ and 12-6-4 data are generated once with AmberTools by `tools/parameterize_ddz.py`
and `tools/lj1264_data.py` (not needed at runtime). See ADR-0010.

## Variants

List variants in the manifest (`"variants"`, or the review page's variant editor, e.g.
`A:468 R>Q`), then run `simprep variants STRUCTURE --manifest M --out DIR`. The first run
writes `DIR/manifest.json` with every side-chain rotamer candidate as a `variant_build`
finding (exit 3); decide them on the review page, then run again with that manifest to
get `wt/` (plus `wt_modelled/` when gaps are modelled), one directory per variant,
`variant_record.json` and `variant_report.md`.
Side chains are built from cited rotamer data in `knowledge/side_chains.yaml`; clashes
flag and do not stop a build. Each variant and a matched wild type are then relaxed around
the site with OpenMM (`pip install -e '.[relax]'`; protocol in the manifest's `relaxation`,
defaults in `knowledge/relaxation.yaml`), giving `wt_relaxed_<site>/` and relaxed variant
directories (rigid builds kept under `unrelaxed/`); clashes that remain are reported as
warnings. See ADR-0005, ADR-0006 and `notebooks/02_variants.ipynb`.

Layout: `schema/` (JSON Schemas, the source of truth), `knowledge/` (versioned YAML
rules and residue mappings), `src/simprep/` (parser, detectors, severity, manifest, prep,
CLI), `frontend/` (review page), `modules/` + `notebooks/` (Colab back end),
`tests/panel/` (regression structures, independently derived expected findings, prep
fixtures).
