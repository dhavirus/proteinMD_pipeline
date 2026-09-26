# ADR-0006: Relaxation of variant sites with OpenMM

- Status: accepted (TASK-006)
- Date: 2026-09-26

## Context

TASK-005 builds variant side chains rigidly; for IDS R468W every rotamer clashes. The
maintainer decided that clashes only flag, that R468W needs a relaxation step, that
waters are kept and relax, and accepted the TASK-006 defaults (OpenMM, Amber ff14SB +
TIP3P, vacuum with strong restraints, 6 A shell, Reference platform, relax every variant
and its matched wild type). TASK-004 decision 2 requires an ADR for any modelling tool.

## Decisions

1. **Engine: OpenMM 8.6.1** (PyPI, MIT/LGPL), as the optional extra `simprep[relax]`.
   It is imported only in `simprep.relax.openmm_run`, lazily, with an install hint; the
   rest of simprep runs without it. Force field `amber14-all.xml` (ff14SB, Maier et al.
   2015, doi:10.1021/acs.jctc.5b00255) + `amber14/tip3p.xml`; the protocol (files,
   shell, restraint, tolerance, iterations, cutoff, platform, pH, seed) lives in
   `knowledge/relaxation.yaml` and is copied into the manifest's `relaxation`, which
   studies edit and every output records.
2. **Topology from our model, not from a file.** `Structure` -> OpenMM topology directly:
   heavy atoms of standard polymer residues and waters; chains split wherever the
   peptide bond is missing (C-N > 2 A, gaps, left-out residues), so no false bonds;
   `ignoreExternalBonds` for the split ends; explicit `CYX` for disulfide cysteines
   (otherwise CYM also matches); `rigidWater=False` (rigid waters among fixed, massless
   atoms silently stopped the minimizer in the spike).
3. **What moves.** The shell is every residue or water with a heavy atom within 6 A of
   a mutated residue, taken as the union over the wild type and all variants at that
   site, so the relaxed wild type and each relaxed variant move the same residues.
   Side chains (not N, CA, C, O, OXT) and whole waters move; mobile heavy atoms are held
   by a 10,000 kJ/mol/nm^2 restraint to their start, except the mutated residues' side
   chains. Everything else is fixed (zero mass).
4. **Template-less chemistry (refines TASK-006 decision 4).** Residues without an Amber
   template (ALS, ions, ligands, glycans) are left out of the calculation and recorded.
   Decision 4 says such a residue within the 1.0 nm cutoff of a mobile atom stops
   relaxation. At 5FQL A:468 that rule would stop every run: ALS A:84 is 5.1 A from the
   mobile side chain of Tyr348 (12-14 A from residue 468 itself, the distance the spec
   quoted). Instead, a shell residue with a movable atom within the cutoff of
   template-less chemistry is **held fixed** (recorded as `frozen`), so no mobile atom
   ever feels the missing chemistry, which is what decision 4 protects against; the hard
   stop remains when a mutated residue's own side chain is within the cutoff. At 5FQL
   A:468 this freezes Tyr348, Tyr466, Asp478 and Lys479 (all near ALS) in every system.
5. **Hydrogens.** OpenMM adds hydrogens for the calculation (pH 7.0) and they are
   discarded. A residue whose heavy atoms moved (>= 0.001 A, the file precision) comes
   back without hydrogens and gets a `protonation` work-order item; unmoved residues keep
   their deposited hydrogens.
6. **Determinism.** OpenMM places new hydrogens at random positions (Python `random`) and
   minimizes them on the default, multithreaded platform. `random` is seeded from the
   protocol, and both the hydrogen placement and the minimization run on the Reference
   platform: two runs give byte-identical files (tested). Cost on this machine: about
   20 s per system at 5FQL (3 systems about 60 s).
7. **Clashes only flag.** `variant_build` no longer switches to `expert_review` when every
   candidate clashes, and a clashing rotamer is built. After relaxation, remaining
   clashes at a mutated residue (same criterion as ADR-0005) become a
   `relaxation.residual_clashes` finding, base severity warn, recorded with the system.
8. **Outputs.** `wt/` (prepared, unrelaxed), `wt_relaxed_<site>/`, `<variant>/`
   (relaxed) with `<variant>/unrelaxed/` (the rigid build); `variant_record.json` gains
   the protocol, per system what moved, what was frozen or left out, clashes before and
   after, the energy (with the caveat that fixed-fixed terms dominate it) and the final
   force on mobile atoms.

## Consequences

- R468W is built: 10 site clashes before, 2 after (5FQL, default protocol), reported as a
  warning. The relaxed wild type moves at most 0.44 A.
- Vacuum electrostatics pull charged side chains together when restraints are weak (the
  spike: Lys486 up to 4 A without restraints). Weakening the restraint in a manifest
  should come with an implicit-solvent model, which would need its own spike and ADR.
- CI installs OpenMM (a large wheel); the relaxation tests take about 2 minutes.
- A mutation next to template-less chemistry (for example a metal ligand) cannot be
  relaxed until that chemistry is parameterized.
