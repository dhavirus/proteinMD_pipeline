# TASK-006: Local relaxation of variant sites (and the matching wild type)

Read `CLAUDE.md`, ADR-0004 and ADR-0005 first. TASK-005 builds side chains rigidly in
the unrelaxed wild-type environment; for IDS R468W every rotamer clashes, so no R468W
system exists yet. This task adds a local, restrained energy minimization around each
mutation site, applied identically to the wild type, so that matched WT and variant
systems can be built even where a side chain does not fit as placed.

Status: **draft** for maintainer review. Items marked (Q) are open questions (end of
file). Do not start implementation until they are answered or accepted as defaults.

## Maintainer decisions carried in from the TASK-005 review (2026-09-26)

1. R468W needs a relaxation step.
2. The clash rule only flags. The common 0.4 A overlap setting stays (with the polar
   exception of ADR-0005); a clash no longer stops a build. This replaces TASK-005
   decision 4: `choose_rotamer` may name a clashing rotamer, and relaxation follows.
3. Waters are kept, and they relax with the site.
4. Record accounting by atom name (ADR-0005 decision 7) stays.

## Spike: what a restrained local minimization does at 5FQL A:468 (this session)

Setup: OpenMM 8.6.1 (PyPI), Amber14 protein force field (`amber14-all.xml`, ff14SB) +
TIP3P waters, no solvent model, non-bonded cutoff 1.0 nm. Input: the TASK-005 5FQL WT
system and variants R468Q (rotamer mm-40) and R468W (m95, the fewest-clash rotamer).
Hydrogens are added by OpenMM for the calculation only. Mobile: side chains (not
N, CA, C, O) of residues with a heavy atom within 6 A of residue 468, plus waters in that
shell; everything else is held fixed. Residues without an Amber template (ALS A:84 at
12-14 A, Ca2+ at 14-16 A, Cl- ions) are left out of the calculation.

Clash counts use the TASK-005 criterion on the mutated side chain.

| restraint on mobile heavy atoms (except the mutated side chain) | system | site clashes before -> after | heavy atoms moved | largest move |
|---|---|---|---|---|
| none | WT | 0 -> 0 | 117 | 3.83 A (Lys486 NZ) |
| none | R468W | 10 -> 2 | 118 | 3.99 A (Lys486 NZ) |
| 1,000 kJ/mol/nm^2 | WT | 0 -> 0 | 117 | 2.31 A |
| 1,000 kJ/mol/nm^2 | R468W | 10 -> 2 | 119 | 1.58 A |
| 10,000 kJ/mol/nm^2 | WT | 0 -> 0 | 115 | 0.40 A |
| 10,000 kJ/mol/nm^2 | R468Q | 0 -> 0 | 106 | 0.57 A |
| 10,000 kJ/mol/nm^2 | R468W | 10 -> 4 | 118 | 1.25 A |

(The last three rows are from the deterministic setup below; residual R468W contacts:
NE1 - Trp502 CZ3 2.73 A, CD1 - Trp502 CZ3 2.96 A, CZ3 - Tyr264 O 2.82 A, CD1 - Pro467 C
3.00 A.)

What the spike taught:

- **It works, and it is cheap.** 16-26 s per system on one CPU core; R468W goes from 10
  clashes to 2-4 marginal ones.
- **Without restraints the protocol itself distorts the structure.** In vacuum, Lys479
  and Lys486 swing up to 4 A toward Asp484 in the WT as much as in the variant (charged
  side chains without solvent). Restraints, or a solvent model (Q2), are needed, and
  the WT must go through the same protocol, or protocol artefacts become "variant
  effects".
- **Traps found and fixed:**
  - Rigid TIP3P constraints on fixed (massless) waters silently stop the minimizer:
    nothing moved until `rigidWater=False`.
  - Chains must be split at real gaps (missing residues, excluded residues), or the
    force field bonds across them.
  - Disulfide cysteines need explicit CYX templates when external bonds are ignored.
- **Determinism is not free.** OpenMM places hydrogens at random positions (Python's
  `random`) and minimizes them on the default (multithreaded) platform. Two runs
  differed by up to 0.32 A in heavy atoms. Seeding `random`, and running both the
  hydrogen placement and the minimization on the Reference platform, gave
  byte-identical output twice; the CPU platform did not.

## Scope in one sentence

For every variant site, relax a shell of side chains and waters around the mutated
residue with a restrained force-field minimization, do exactly the same around the same
site in the wild type, record every setting and every moved atom, and flag whatever
clashes remain.

## Deliverables

### 1. Knowledge base and manifest
- `knowledge/relaxation.yaml` (new schema), holding the protocol **defaults** with
  sources:
  - force field files;
  - shell radius;
  - which atoms move;
  - restraint constant;
  - minimizer tolerance and iteration limit;
  - platform;
  - hydrogen pH;
  - random seed.
- Manifest `relaxation` (optional): the protocol actually used for this study. It is
  copied from the defaults when absent, edited per study (CLAUDE.md: per-study
  judgment), and recorded in every output.
- `variant_build` rule: the recommendation no longer switches to `expert_review` when
  every candidate clashes. It recommends the fewest-clash rotamer, and the finding
  states that relaxation will follow (maintainer decision 2). `expert_review` remains
  an option.
- New rule family `relaxation`, finding `relaxation.residual_clashes`:
  - base severity **warn** (clashes only flag);
  - one finding per relaxed system that still has clashes, with the atoms and
    distances;
  - options: `accept`, or `expert_review`.

### 2. `src/simprep/relax/` (I/O edge; OpenMM confined here)
- `shell.py` (pure): mobile and fixed atom sets from the `Structure` and the protocol.
  It also lists residues excluded from the calculation (no force-field template) with
  their distance to the site.
- `openmm_run.py`:
  - Structure -> OpenMM topology: chains split at gaps, CYX for disulfides, waters
    included.
  - Seeded hydrogens on Reference; restraints; minimization on Reference.
  - Back to heavy-atom coordinates.
  - The hydrogens OpenMM added are discarded.
- Refusal rule (unknown chemistry is a hard stop, CLAUDE.md): a residue without a
  template inside the non-bonded cutoff of any mobile atom stops relaxation, naming the
  residue (Q4). Beyond it, it is excluded and recorded.
- Output of a relaxed system:
  - moved heavy atoms take their relaxed coordinates;
  - every residue with a moved atom loses its deposited hydrogens (they no longer match
    its heavy atoms) and gets a `protonation` work-order item, as the mutated residue
    already does (ADR-0005).

### 3. Matched wild type
- For each variant site, `wt_relaxed_<site>/` is the prepared wild type relaxed with
  the identical protocol and shell. The shell is the union of the WT and variant
  shells, so both move the same residues (Q3).
- A test asserts that the relaxed WT and the relaxed variant differ only inside that
  shell.

### 4. CLI, record, report
- `simprep variants` relaxes when the manifest's `relaxation` says so (default: on,
  Q5). Outputs:
  - `wt/` (unrelaxed, as today);
  - `wt_relaxed_<site>/`;
  - `<variant>/` (relaxed; the rigid build is kept as `<variant>/unrelaxed/`).
- `variant_record.json` (schema change) gains, per system:
  - the protocol and the OpenMM version;
  - moved residues and atoms, largest and RMS displacement;
  - clashes before and after;
  - residues left out of the calculation.
- Record accounting counts moved atoms as modified.
- Energies are recorded with a caveat. Fixed-fixed terms dominate the total, so it is
  not a comparable quantity; only the change and the final force on mobile atoms are
  meaningful.

### 5. Dependency and CI
- `openmm==8.6.1` as an optional extra `simprep[relax]`. CI installs it and runs the
  relaxation tests offline. On Colab the notebook installs the extra.
- ADR-0006 records the tool choice (TASK-004 decision 2 requires one).

### 6. Front end
- The `relaxation` family is shown like the others. The review page shows the
  manifest's `relaxation` settings read-only; they are edited in JSON in v0.1 (Q6).

### 7. Tests (offline, CI with the extra installed)
- Unit tests:
  - shell selection on hand-built structures;
  - chain splitting at gaps, and the CYX assignment;
  - the refusal of a template-less residue inside the cutoff;
  - hydrogens removed from moved residues;
  - the accounting.
- 5FQL (expected outcomes written by hand from the spike above, `reviewed: false`):
  - R468W is built and relaxed;
  - site clashes drop from 10 to at most 4, and any residual ones are a warn finding;
  - relaxed WT heavy atoms move at most 0.5 A;
  - relaxed WT and relaxed variant differ only inside the shell;
  - two runs give byte-identical files.
- Runtime budget: three relaxations of about 20 s each in CI.

## Acceptance criteria
- `simprep variants` on 5FQL (TASK-004 fixture decisions, variants R468Q and R468W with
  their fewest-clash rotamers) writes the following, all schema-valid and deterministic:
  - `wt/` and `wt_relaxed_A468/`;
  - `R468Q/` and `R468W/`;
  - the records.
- R468W's residual clashes, if any, appear as a `relaxation.residual_clashes` warn
  finding with atoms and distances. Nothing is silently accepted.
- A template-less residue within the cutoff stops relaxation with a message naming it.
- CI green without network (OpenMM installed from PyPI during setup).

## Out of scope
MD, equilibration, solvation boxes, protonation / pKa (the hydrogens used during
minimization are discarded), parameters for non-standard residues (ALS stays a
chemistry/parameterization work item), backbone relaxation beyond the restraints above,
alchemical/FEP setup, GPU.

## Open questions for the maintainer
1. **Engine and force field (default: OpenMM 8.6.1, Amber ff14SB via `amber14-all.xml`
   + TIP3P, in vacuum with restraints).** ff14SB: Maier et al. 2015,
   doi:10.1021/acs.jctc.5b00255 (PubMed 26574453, read in this session); OpenMM 8:
   Eastman et al. 2024, doi:10.1021/acs.jpcb.3c06662 (PubMed 38154096). The TIP3P
   reference and the file contents of the OpenMM force-field XMLs are `[VERIFY]`.
2. **Electrostatics (default: vacuum + strong restraints, as in the spike).** The
   alternative is an implicit solvent (OpenMM's GBn2 file with amber14), which should
   remove the Lys-Asp collapse and allow weaker restraints. It was not tested in this
   session, so it would need another spike.
3. **Shell and restraints (default):**
   - mobile: side chains (not N, CA, C, O) and waters with a heavy atom within 6 A of
     the mutated residue;
   - restraint: 10,000 kJ/mol/nm^2 on every mobile heavy atom except the mutated side
     chain;
   - minimizer: tolerance 1 kJ/mol/nm, at most 5,000 iterations;
   - the WT uses the union shell.

   The alternatives are a larger shell, a mobile backbone, or weaker restraints; the
   spike shows WT drift grows as restraints weaken.
4. **Template-less chemistry near a site (default):**
   - stop if such a residue is within 1.0 nm (the non-bonded cutoff) of any mobile atom;
   - otherwise exclude it and record it.

   For 5FQL A:468 the nearest (ALS A:84) is 12-14 A from the site. Alternative:
   generic Lennard-Jones-only parameters for such residues, which is more permissive
   and needs sources.
5. **When to relax (default: always, for every variant and its matched WT).** The
   alternatives are only when a clash remains, or a per-variant choice in the manifest.
   "Always" keeps WT/variant comparisons uniform.
6. **Front-end editing of the protocol (default: read-only display in v0.1; edit the
   manifest JSON).**
7. **Platform (default: Reference, for byte-identical output).** It costs about 3x CPU
   time (18 s vs 6.5 s per system here). The alternative is the CPU platform, with
   reproducibility tested to a tolerance (for example 0.05 A) instead of bytes.
