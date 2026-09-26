# ADR-0007: Loop modelling v0.1: PDBFixer placement, restrained sterics-only minimization

- Status: accepted (TASK-007; maintainer, 2026-09-26: the restraint + sterics-only
  protocol is accepted, implicit solvent not pursued; test time of about 8 minutes kept)
- Date: 2026-09-26

## Context

TASK-007 builds unobserved internal loops decided `model_loop`, starting with the 5FQL
SD1-SD2 linker A:444-453 (FRDLEEDPYL). The maintainer decided the system form (intact
single chain), the source (PDBFixer 1.12.0 placement, then an OpenMM minimization on the
Reference platform with the loop mobile and the flanks restrained, every random source
seeded, stop if not byte-identical), the marking (occupancy 0.00, a recorded B, every
modelled atom listed), no caps, and the placement of the stage (prep -> model -> variants
-> relax). TASK-004 decision 2 requires an ADR for any modelling tool. This ADR records the
implementation and what the session's measurements forced.

## Measurements (TASK-007 session, 5FQL, prep fixture decisions)

- **Determinism.** PDBFixer 1.12.0 runs Langevin dynamics on the new atoms when they come
  closer than 0.13 nm. Its integrator is unseeded unless `addMissingAtoms(seed=)` is given,
  and the default platform is multithreaded. With the seed and the Reference platform, two
  runs (also with a different `PYTHONHASHSEED`) gave byte-identical files; with the same
  seed on the default platform they differed. A different seed gives a different loop.
- **Chirality.** PDBFixer's placement had four D residues (Phe444, Leu447, Tyr452, Leu453:
  CA-N-C-CB improper about -30 degrees, against +34 to +38 for the observed flanks);
  PDBFixer minimizes new atoms with its own soft force field (`soft.xml`), whose 22
  impropers all sit on sp2 centres (planarity); nothing holds an sp3 CA's chirality. A plain ff14SB minimization keeps them D and produced cis peptides
  (omega about 4 degrees at three bonds).
- **Vacuum electrostatics.** With full non-bonded terms in vacuum the loop (four acidic side
  chains) stretched: CA-CA up to 3.96 A, C-N-CA up to 132 degrees, a 85 degree CA-C-N at
  the Arg443 flank. Zeroing only the long-range charges (keeping 1-4 terms) was worse;
  zeroing all charges left CA-CA up to 3.92 A.
- **Cost.** Reference platform: about 80-90 ms per energy evaluation of the whole system;
  a converged minimization (tolerance 1 kJ/mol/nm) needs about a thousand. A tolerance of
  10 was fast but not converged (max force 158 kJ/mol/nm, two CA-CA failures, a proline
  with phi +48).

## Decisions

1. **Where the stage lives.** `src/simprep/model/`: `loops.py` (which gaps, insertion,
   shell, torsion restraints; pure), `checks.py` (geometry and contacts; pure),
   `placement.py` (PDBFixer; edge), `run.py` (stage, record), `report.py`. `simprep model`
   writes `wt/` (plain prep), `wt_modelled/`, `model_record.json`, `model_report.md`;
   `simprep variants` builds the same model in-process (both passes, since rotamer
   candidates are evaluated on the system they are built on) and builds every variant on
   it. `simprep prep` never invents coordinates (ADR-0004 unchanged).
2. **`model_loop` is `apply` with `prep_stage: modelling`.** The rule schema now allows an
   apply option to name the later stage that carries it out; prep lists it (action
   `apply`, a `modelling` work item) and applies nothing. The model stage consumes the
   item. Terminal-tail and whole-chain `model_loop` stay `defer` (not built in v0.1). A
   test checks that every apply option either has a prep operation or belongs to the
   modelling stage.
3. **`truncate` means charged termini.** Its work item now ends with the terminus the
   topology stage builds, for example "charged N-terminus (NH3+) at THR A:34"; an internal
   truncate names both ends.
4. **Placement.** PDBFixer 1.12.0 (optional extra `simprep[relax]`, with OpenMM) reads the
   prepared system as PDB with SEQRES, builds only the requested gaps (checked against
   the flanks and the sequence), with `seed` and the Reference platform from the
   protocol; the installed version must equal the protocol's. Only the new residues'
   heavy atoms are taken back; every existing atom keeps its coordinates.
5. **Minimization** (reuses `relax()` of ADR-0006 with a loop shell: every heavy atom of
   the loop free, the two flanks restrained at 10,000 kJ/mol/nm^2, everything else fixed):
   - a bonded-only pre-stage (non-bonded forces removed) so D centres and cis peptides can
     be corrected without clashes in the way;
   - torsion restraints (1,000 kJ/mol/rad^2) in both stages: omega 180 degrees for every
     peptide bond from flank to flank, CA-N-C-CB improper +34 degrees (the observed
     flanks' median) for every non-glycine loop residue. Without them in the final stage,
     clash forces flipped four centres back to D;
   - final stage **sterics only**: Lennard-Jones (ff14SB sigma/epsilon, cutoff 0.6 nm with a
     switch from 0.5 nm) between mobile atoms and the atoms within the cutoff plus a
     1.5 nm margin of their start; no charges, no 1-4 terms, no fixed-fixed pairs. A mobile
     atom moving farther than the margin stops the stage, so the partner list is exact.
     Electrostatics are left to the solvated equilibration (a later task).
   - A non-finite coordinate after any minimization now raises instead of silently
     keeping the old positions (this also guards TASK-006 relaxation).
6. **Checks and findings** (family `modelling`, `knowledge/modelling_checks.yaml`,
   thresholds in the protocol): `modelling.junction_geometry` (blocking) for any peptide
   bond from flank to flank with C-N outside 1.33 +/- 0.05 A or CA-CA outside
   3.8 +/- 0.1 A, or a non-glycine residue with a non-positive improper (D); a blocking
   finding means no `wt_modelled/` is written, the record has `status: rejected`, and the
   command exits 3. `modelling.loop_contacts` (warn) lists clashes of loop atoms with
   anything except the same residue and sequence neighbours (ADR-0005 criterion). phi,
   psi and the improper of every modelled residue are always in the record and report.
7. **Marking and record.** Modelled atoms: occupancy 0.00, B 0.00 (protocol `marking`).
   The loop residues leave the unobserved-residue annotation and count as `added`. The
   record lists every modelled residue with its atoms, the source string (PDBFixer and
   OpenMM versions), the note "no experimental support", geometry, the minimization
   (ADR-0006 format), record accounting against the input, findings, and the work order
   (prep's minus the built gaps, plus protonation items: modelled residues have no
   hydrogens; flanks that moved lost theirs).

## Consequences

- 5FQL A:444-453: all 11 peptide bonds pass (C-N 1.33-1.35 A, CA-CA 3.79-3.88 A, omega
  within 3 degrees of 180), every residue L, 42 loop clashes after placement -> 0, flanks
  and loop moved at most 9.3 A. Some backbone dihedrals are unusual (Asp446 phi about -5,
  several positive phi): the loop is a clash-free, chemically sane starting model, not a
  prediction of its conformation.
- A build takes about two minutes on this machine (PDBFixer about 33 s, minimization about
  75 s); `simprep variants` builds it in both passes. The modelling panel tests build it
  three times (about six minutes).
- The rigid-variant and relaxation panel tests decide 444-453 `truncate` (fixture
  `prep_overrides`) so they do not rebuild the loop.
- Another seed gives another loop. Loop ensembles are out of scope (TASK-007).
