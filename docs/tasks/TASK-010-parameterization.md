# TASK-010: Parameterization v0.1: DDZ, the 12-6-4 Ca2+ site, chloride

Read `CLAUDE.md`, ADR-0006 (OpenMM, amber14-all + TIP3P), ADR-0008 (DDZ), ADR-0009
(protonation) first. After TASK-009 every system has its hydrogens and recorded states;
the work order holds the parameterization items of 5FQL:

| item | manifest decision |
|---|---|
| `metals/A:1551` Ca2+ | `nonbonded_12_6_4` (now 5 protein ligands: D45, D46, D334, H335, DDZ OG2) |
| `unrecognized/group/A:1567-1569` Cl- | `parameterize_manually` |
| DDZ A:84 (gem-diol FGly) | no template in any shipped force field |

Status: **draft** (open questions below, each with a default, for the maintainer).

## What the tools ship (read in this session)

- **OpenMM 8.6.1 `amber14/tip3p.xml`**: Cl- from Joung & Cheatham 2008 (`frcmod.ionsjc_tip3p`)
  and Ca2+ from Li et al. 2013 (`frcmod.ions234lm_126_tip3p`, 12-6 only). No 12-6-4.
- **OpenMM's prmtop reader** supports 12-6-4: it reads `LENNARD_JONES_CCOEF` and adds a
  `-c/r^4` CustomNonbondedForce.
- **AmberTools 26.0** (conda-forge; installs in about 9 min):
  - `frcmod.ions234lm_1264_tip3p`: Ca2+ Rmin/2 1.642 A, epsilon 0.10185975 kcal/mol,
    "12-6-4 set for TIP3P water from Li and Merz, JCTC, 2014, 10, 289" (file comment);
  - ParmEd `add1264`: C4(ion, water O) for TIP3P, Ca2+ 87.3 kcal/mol/A^4; C4 with any other
    atom = C4(water) / 1.444 x polarizability of its Amber atom type x a tuning factor
    (default 1.0);
  - `lj_1264_pol.dat`: polarizabilities per Amber atom type ("Miller JACS 112, 8533 (1990)";
    hydroxyl HO/ho 0.000; Ca2+ 0.477 from B3LYP);
  - antechamber AM1-BCC charges + GAFF2 types for DDZ in 3 s.
- `[VERIFY]` Li & Merz 2014 (JCTC 10:289), Li et al. 2013 (JCTC 9:2733) and Miller 1990 are
  cited from file headers only; PubMed was unavailable in this session.

## Spike (this session)

- **ff14SB has no gem-diol carbon.** OpenMM's `protein.ff14SB.xml` has no angle OH-2C-OH, no
  torsion OH-2C-OH-HO, and no H2 hydrogen type (H on a carbon with two electronegative
  neighbours). Serine-like typing covers the backbone and CX-2C-OH, H1-2C-OH only. DDZ
  needs parameters beyond ff14SB whichever route is taken.
- AM1-BCC on the free DDZ (CCD ideal coordinates, net charge 0) works: OG1/OG2 -0.629,
  HG1/HG2 +0.429, CB +0.295, HB +0.081 (GAFF2 types c3, oh, ho, h2). A residue inside a
  chain needs a capped fragment (ACE-DDZ-NME) instead, with the backbone matched to ff14SB.

## Scope in one sentence

Give every atom of the protonated 5FQL systems (wild type and variants) recorded,
sourced parameters, with the Ca2+ site as the manifest decided (12-6-4), and check that
nothing is left unparameterized.

## Deliverables (depend on decision 1)
- DDZ template (atom types, charges, bonded terms) as data with its derivation recorded.
- The 12-6-4 Ca2+ model (12-6 terms and C4 per atom type) as data with sources.
- A `simprep parameterize` check: builds the OpenMM System for each protonated system,
  hard-stops on any atom or bonded term without parameters, records net charge and the
  parameter sources per residue. Solvation and the MD-ready system are TASK-011.
- Tests: the DDZ template matches its derivation; C4 terms reproduce ParmEd's `add1264`
  for a small system (reference numbers committed, computed once with AmberTools);
  5FQL builds with no missing parameters; net charge as hand-derived from the states.

## Open questions (defaults in bold)
1. **Build route.** **OpenMM force-field XML (as now) with committed data: the DDZ template
   and the 12-6-4 C4 table are generated once with AmberTools by a script in `tools/`
   (provenance recorded), and simprep adds the C4 term as a CustomNonbondedForce,
   validated against ParmEd/prmtop energies.** Runtime stays pip-only; CI and Colab need
   no conda. Alternative: the Amber route at runtime (tleap + ParmEd `add_12_6_4` ->
   prmtop -> OpenMM), which is the reference implementation but makes AmberTools (conda)
   a runtime dependency of CI and Colab.
2. **DDZ charges.** **AM1-BCC on a capped ACE-DDZ-NME fragment; backbone atoms (N, H, CA,
   HA, C, O) take ff14SB serine charges and types; the side-chain charges are shifted
   evenly so the residue is neutral.** Alternatives: RESP (HF/6-31G*, needs a QM code such
   as psi4), or charges by analogy with serine (no derivation, weakest).
3. **DDZ side-chain types.** **GAFF2 types for CB, HB, OG1, OG2, HG1, HG2 with parmchk2
   supplying the cross terms to the ff14SB backbone.** Alternative: ff14SB types plus only
   the missing terms borrowed from GAFF2 (smaller change, mixed provenance per term).
4. **12-6-4 details.** **Li-Merz 12-6-4 TIP3P set for Ca2+; C4 for every atom by
   polarizability (tuning factor 1.0), as ParmEd does; Cl- stays Joung-Cheatham 12-6 (its
   C4 with Ca2+ from its polarizability).** Alternative: C4 only for the ligating atoms.
5. **Water model.** **TIP3P** (consistent with the 12-6-4 set chosen and ADR-0006).

## [VERIFY] introduced
- The three ion/polarizability references above (read from AmberTools files, not papers).
