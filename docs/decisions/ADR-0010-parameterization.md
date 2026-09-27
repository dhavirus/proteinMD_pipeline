# ADR-0010: Parameterization v0.1: generated DDZ parameters, 12-6-4 Ca2+ in OpenMM

- Status: accepted (TASK-010; maintainer, 2026-09-27: all spec defaults)
- Date: 2026-09-27

## Context

After protonation (ADR-0009) three pieces of 5FQL have no ready parameters: DDZ A:84
(gem-diol FGly), the Ca2+ A:1551 the manifest decided `nonbonded_12_6_4`, and the Cl-
ions. OpenMM ships ff14SB, TIP3P, Joung-Cheatham Cl- and Li-Merz 12-6 Ca2+, but no 12-6-4
and nothing for DDZ; OpenMM's ff14SB has no gem-diol carbon terms (no OH-2C-OH angle,
no OH-2C-OH-HO torsion, no H2 type). AmberTools 26.0 (conda) has the 12-6-4 data and
antechamber, but is not pip-installable.

## Decisions

1. **OpenMM force-field XML at runtime, AmberTools only offline.** Two scripts in `tools/`
   run once with AmberTools (and RDKit) and write committed data:
   - `tools/parameterize_ddz.py` -> `knowledge/forcefield/ddz.xml` (+ provenance);
   - `tools/lj1264_data.py` -> `knowledge/forcefield/lj1264.yaml` and ParmEd reference C4
     values in `tests/panel/parameter_fixtures/lj1264_reference.yaml`.
   CI and Colab stay pip-only. Both files are in the knowledge-base hash, so a
   regeneration changes provenance. Rerunning the DDZ script reproduced the file byte
   for byte.
2. **DDZ.**
   - Charges: AM1-BCC on ACE-DDZ-NME (RDKit embedding with a fixed seed; L checked
     with simprep's CA-N-C-CB convention).
   - Backbone (N, H, CA, HA, C, O): ff14SB serine types and charges, read from OpenMM's
     own file.
   - Side chain: charges shifted evenly (+0.0102 e per atom) to a neutral residue;
     GAFF2 types (c3, h2, oh, ho), under unique `ddz-*` OpenMM types.
   - Bonded terms: every term that touches a side-chain atom takes GAFF2 parameters
     (parmchk2 for missing ones: 4 bonds, 8 angles, 13 torsions), so those terms apply
     to DDZ only.
3. **12-6-4 as ParmEd computes it.**
   - The ion's 12-6 terms become the 12-6-4 TIP3P set (Ca2+ Rmin/2 1.642 A, epsilon
     0.10185975 kcal/mol).
   - A CustomNonbondedForce adds -C4/r^4 between the ion and every other particle, with
     C4 = C4(ion, water O) / 1.444 x polarizability x 1.0 (Ca2+: 87.3 kcal/mol/A^4).
   - Each particle's polarizability is found by its LJ terms. ParmEd assumes the same
     (types with equal LJ terms share a polarizability), and simprep checks it.
   - Unit test: C4 per Amber type matches ParmEd's add12_6_4 on a tleap system.
   - Ion-ion pairs are left out; this is enough for one ion per model, as in 5FQL.
4. **Which ion gets which model is the manifest's decision.** `metals` findings decided
   with an option listed in `metal_models` (`nonbonded_12_6_4 -> lj1264`). Other ions
   keep OpenMM's TIP3P parameters (Cl- Joung-Cheatham 12-6).
5. **`simprep parameterize DIR --manifest M`** reads the protonation record in DIR (from
   `protonate` or `variants`) and for each system:
   - builds the OpenMM System in vacuum without a cutoff;
   - stops if a residue has no template, or if any bond or angle has no parameters
     (OpenMM skips unmatched angles silently);
   - writes `<name>_parameterized/system.xml` and `parameterization_record.json` (net
     charge, residues per parameter source, ion models, checks).
   Solvation, cutoffs and the MD-ready system are TASK-011.

## Consequences

- 5FQL wild type (modelled, pH 7.2): 8,770 particles, net charge -15 (equal to the formal
  charges of the recorded states, Arg, termini and ions); 516 ff14SB residues, DDZ, 213
  TIP3P waters, 3 Cl-, Ca2+ 12-6-4 with 8,769 C4 pairs; all 8,686 bonds and 15,151 angles
  parameterized; about 2 s.
- A chain end left by an internal `truncate` has no OXT and no template: parameterize
  stops and says so; the topology stage (TASK-011) must add terminal atoms.
- `[VERIFY]` Li & Merz 2014, Li et al. 2013, Joung & Cheatham 2008 and Miller 1990 are
  cited from AmberTools / OpenMM file headers; PubMed was unavailable in this session.
- `system.xml` is about 12 MB per system.
