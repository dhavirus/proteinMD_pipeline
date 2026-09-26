# TASK-004: Prep v0.1: apply manifest decisions to coordinates

Read `CLAUDE.md` (especially "Colab back-end conventions"), ADR-0001 and TASK-002 first.
TASK-001 to 003 produce a manifest in which every finding has a recorded decision, and
`simprep manifest status` exits 0 only when no blocking finding is undecided. Nothing yet
*applies* those decisions. This task builds the first back-end stage: it turns the
deposited structure plus a complete manifest into prepared system file(s) and a work
order for the later stages, with every input record accounted for.

Status: **accepted** (maintainer, 2026-09-26: defaults accepted for every open
question; the resolved decisions are listed at the end of this file).

## Scope in one sentence

Prep v0.1 **selects, removes and records**. It never invents coordinates: loop
modelling, capping, protonation, chemistry edits and parameterization are later stages
(decision 2). It hands them a machine-readable work order, so no decision is dropped on the way.

## Constraints carried over

- CLAUDE.md: all coordinate work happens in Python, driven by the manifest. Nothing is
  excluded or modified without a recorded decision; counts in and counts out are reported.
  Provenance on everything.
- The front end never transforms coordinates (ADR-0001); prep is the first place they
  change.
- Colab conventions 1–8 (Drive layout, `run_<yymmdd_hhmmss>` runs with `config.json`
  and `run_id`, one config dict, explicit seeds, local I/O then bulk copy, resumable,
  pinned repo commit, `importlib.reload`, `nbformat.validate`).
- No GPU and no MD here: prep must run on CPU in seconds for the panel.

## What each decision does in prep

Options come from the knowledge base; this table is the contract. Rule files gain an
option field `prep_action` so the table is data, not code (schema change, examples).

| option | prep_action | effect in v0.1 |
|---|---|---|
| altlocs `highest_occupancy` | `apply` | keep that altloc (ties by altloc id), drop the others, clear altloc ids |
| altlocs `specific_altloc` | `apply` | keep `parameters.altloc` |
| altlocs `keep_ensemble` | `apply` | one output system per altloc id (decision 4) |
| covalent_contacts `add_link` | `apply` | write the bond into the output's `struct_conn`; geometry unchanged |
| covalent_contacts `treat_noncovalent` | `record` | no change; noted in the work order |
| covalent_contacts `exclude_group`, unrecognized `exclude`, nonstandard `exclude_segment` | `apply` | delete the residues (the whole group for glycans) |
| metals: all four models | `defer` | work order item for parameterization |
| missing `model_loop`, `cap_termini` | `defer` | work order item for modelling |
| missing `truncate`, `exclude_chain` | `record` | nothing to delete; the chain break or terminus is noted for topology |
| nonstandard `route_to_topology_builder`, `keep_sulfate_ester`, unrecognized `parameterize_manually` | `defer` | kept as deposited; parameterization item |
| formylglycine `model_gem_diol`, `model_aldehyde` | `defer` | chemistry edit + parameters item (FGly default is gem-diol) |
| `revert_to_parent` | `apply` only with a mapping (decision 3) | rename and delete atoms per `knowledge/residue_mappings.yaml`; refuse without one |
| `expert_review`, `add_rule` | `unresolved` | prep refuses to run (see §2) |

## Deliverables

### 1. `src/simprep/prep/`: the library (pure core, thin I/O)
- `plan.py`: pure: `(Structure, manifest, RuleSet) -> PrepPlan`, the list of actions per
  finding, the systems to build (one, or one per altloc for `keep_ensemble`) and the work
  order. Raises a specific error listing every unresolved decision.
- `apply.py`: pure: `(Structure, PrepPlan) -> list[Structure]` (selection, deletion,
  link records, mapped renames). Operates on the immutable model; no gemmi.
- `write.py` (I/O edge): prepared structure as mmCIF (primary) and PDB (decision 5), with the
  `add_link` bonds in `struct_conn` / `LINK`.
- `record.py`: the prep record (below).

### 2. Gate
- Prep runs only if `manifest status` would exit 0 **and** no decision uses an option with
  `prep_action: unresolved`. `manifest status` learns the same rule: a decision for
  `expert_review` / `add_rule` counts as *undecided* for blocking findings (today it counts
  as decided; that is a TASK-002 gap this task closes).
- The structure's SHA-256 must match the manifest's input, as in `audit`.

### 3. Prep record: `schema/prep_record.schema.json` (new, with examples)
Per output system: input and manifest SHA-256, knowledge-base version and hash, simprep
version and commit, the applied / recorded / deferred action per finding, the work order,
output file SHA-256s, and **record accounting with the `excluded` column in use**:
`in = passed_through + modified + excluded`, per record type, as in the audit counts. A
test asserts the equation for every panel structure.

### 4. CLI: `simprep prep STRUCTURE --manifest M --out DIR`
Writes `system_<n>.cif` (+ `.pdb`), `prep_record.json` and a short `prep_report.md`.
Deterministic: same inputs give byte-identical coordinate files (no timestamps inside
them).

### 5. Colab: `modules/prep_module.py` + `notebooks/01_prep.ipynb` (decisions 1, 6)
- The module wraps the CLI function with the run bookkeeping: `PROJECT` root
  `/content/drive/MyDrive/googleColab_run/simprep/` (demo under `googleColab_test/`),
  `./<yymm>/run_<yymmdd_hhmmss>/`, `config.json` with the manifest and the simprep commit,
  a `run_id`; work in `/content`, then bulk copy to Drive; re-running a run directory
  skips systems whose outputs and record already exist and match (resume).
- The notebook has one config dict at the top (repo commit, input paths, output root,
  seed), clones the repo at that commit, reloads modules with `importlib.reload`, and
  calls one function. The notebook is validated with `nbformat.validate` in CI.
- Local runs use the same function with a temp directory in place of Drive, which is how
  the tests exercise it (no Colab, no network).

### 6. Tests (offline)
- Unit tests per prep action on hand-built structures.
- Panel: for each structure, a manifest fixture that decides every finding (a
  deliberate mix of options, one fixture per structure, committed), and expected outputs
  derived by hand: residue and atom counts after exclusion, no altlocs left unless
  `keep_ensemble`, `add_link` present in the written mmCIF, work order items listed.
- The accounting equation for every panel structure; re-reading each written file with
  gemmi gives the counts in the record.
- Gate: `expert_review` decision → prep refuses and lists it; `manifest status` exits 3
  for a blocking finding decided with `expert_review`.
- Determinism: two runs give identical coordinate files and records (apart from
  timestamps in the record).

## Acceptance criteria
- `simprep prep` on each panel structure with its fixture manifest writes schema-valid
  records and coordinate files whose re-read counts match the record.
- 5FQL with glycans excluded: the seven glycan groups and their 13 branched + 2 NAG
  residues are gone, counted as excluded, and nothing else is missing.
- 1HZH with `add_link` for both candidates: both ASN–NAG bonds appear in the output
  `struct_conn`.
- 1FO8 with `keep_ensemble` on one residue produces two systems differing only there
  (decision 4).
- An `expert_review` decision blocks prep with a message naming the finding.
- CI green without network; notebook passes `nbformat.validate`.

## Out of scope
Loop modelling, capping, protonation / pKa, the FGly chemistry edit, parameterization,
solvation, mutation building (WT vs R468Q/R468W is the next task), MD, and any GPU
work.

## Resolved decisions (maintainer, 2026-09-26)
1. **Code location**: the library lives in `src/simprep/prep/`; `modules/prep_module.py`
   is a thin Colab wrapper around it.
2. **Scope**: prep v0.1 only selects, removes and records. Every modelling step (loops,
   caps, protonation, chemistry edits) and parameterization is deferred through the work
   order; adding a modelling tool later needs an ADR choosing it.
3. **`revert_to_parent`**: applied only for components with an explicit atom mapping in
   `knowledge/residue_mappings.yaml`; without one, prep refuses and names the component.
   No mapping ships in v0.1 (the ALS parent mismatch, ALA in the file vs CYS in the rule,
   stays a human decision).
4. **`keep_ensemble`**: one system per altloc id, applied to every ensemble residue
   together; no combinatorial expansion.
5. **Output formats**: mmCIF (primary) plus PDB.
6. **Notebook**: `notebooks/01_prep.ipynb` ships in this task, validated with
   `nbformat.validate`.
