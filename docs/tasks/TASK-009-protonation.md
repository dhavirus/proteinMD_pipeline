# TASK-009: Protonation v0.1: pKa estimates, protonation states and hydrogens

Read `CLAUDE.md`, ADR-0006 (relaxation: hydrogens added for the calculation and
discarded), ADR-0007 (modelled loop), ADR-0008 (DDZ) first. After TASK-008 the work order
holds `protonation` items: the modelled loop A:444-453 and its moved flanks (no
hydrogens), mutated residues, residues moved in relaxation, and DDZ A:84's two hydroxyls.
This task gives every system one consistent, recorded set of protonation states and all
its hydrogens.

Status: **accepted** (maintainer, 2026-09-27: every default accepted; see the end of this file).

## What the sources say (read in this session)

- Demydchuk et al. 2017 (doi:10.1038/ncomms15786, PMC5472762): sulfate elimination is
  expected to be efficient "only at the optimum lysosomal pH of ~ 4.8"; storage buffer pH
  7.5, crystallization pH 6.1; the authors assigned partial charges with PDB2PQR using
  PROPKA for side-chain pKa values. Misfolded variants are retained by ER quality control.
- Olsson et al. 2011, PROPKA3 (J Chem Theory Comput 7:525-37, doi:10.1021/ct100578z,
  PubMed 26596171): pKa RMSD against experiment 0.79 (Asp, Glu), 0.75 (Tyr), 0.65 (Lys),
  1.00 (His).
- `[VERIFY]` ER lumen pH (about 7.2) and lysosomal pH (4.5-5): textbook values, no source
  read in this session.

## Spike (this session): PROPKA 3.5.1 on the modelled 5FQL wild type (with DDZ, no waters)

Measured on a loop built before the ADR-0007 amendment (numpy 2.5.3, one pass); values
for the linker residues A:444-453 will shift with the amended build, the rest should not.

| fact | value |
|---|---|
| deposited hydrogens | 3,960 riding hydrogens in the prepared system (refined with them); none on the loop, DDZ, moved flanks |
| predicted net charge (folded) | +6.5 at pH 5, -13.5 at pH 7, -18.8 at pH 8 |
| residues within 1 pH unit of their pKa (about the method's RMSD) | 43 at pH 4.8 (mostly Asp/Glu), 12 at pH 7.0, 9 at pH 7.4 |
| at pH 7.4 | D252, D308 (7.46), H179, H226, H342, H518, H533, K135 (7.43), K347 (6.72) |
| active site | D45 -8.1, D46 -0.9, D334 4.95, H138 1.1, H229 0.6, H335 -0.7, K135 7.4, K347 6.7 |
| termini | N+ Thr34 8.2, C- Pro550 3.4 (OXT present) |

PROPKA ignores the Ca2+ (and does not parameterize DDZ), so the active-site values are
not credible as they stand: the Ca2+ ligands D45, D46, D334 (carboxylates) and H335 (via
NE2) are bound to a metal, and the lysines around the former sulfate come out 3-4 units
below normal. The pH choice moves the net charge by about 20 e: it is the study's
decision, not a default.

## Scope in one sentence

For each system (wild type, modelled wild type, every variant), estimate pKa values,
turn them into per-residue protonation states at the study's pH with recorded rules and
decisions, and add every hydrogen once, deterministically.

## Deliverables

### 1. Manifest and knowledge base
- Manifest `protonation` (required for this stage): `ph` with a rationale (decision 1),
  and the method (copied from `knowledge/protonation.yaml` like `relaxation`).
- `knowledge/protonation.yaml`: method (PROPKA 3.5.1), ambiguity window (1.0 pH unit,
  from the RMSD above), distance within which PROPKA-blind chemistry makes an estimate
  unreliable (decision 4), and the metal-ligand rules.
- `knowledge/protonation_checks.yaml`, family `protonation`: findings for ambiguous and
  unreliable residues (decision 4); hydrogen definitions for DDZ (HG1, HG2 on OG1, OG2).

### 2. `src/simprep/protonate/`
- pKa estimate (edge: PROPKA as a library, deterministic).
- States (pure): pKa vs pH, then the metal rule (a carboxylate or His nitrogen linked to a
  metal in `struct_conn` is not protonated; His with NE2 bound is HID, with ND1 bound HIE),
  then recorded decisions.
- Hydrogens (edge): OpenMM `Modeller.addHydrogens` with explicit per-residue variants
  (ASH/GLH/HIP/HID/HIE/LYN/CYM), seeded, Reference platform; waters get hydrogens.
- Record: per residue the pKa, the state, and why (pKa, metal rule, decision); counts.

### 3. Findings (family `protonation`)
- `protonation.ambiguous_state` (warn): |pKa - pH| < window; options: the predicted state
  (default), the other state.
- `protonation.unreliable_estimate` (warn): a titratable residue within the distance of
  chemistry PROPKA ignores (metal ions, non-standard residues), not covered by the metal
  rule; options: the standard state at that pH (default), the PROPKA state, expert review.

### 4. Tests (offline)
- Unit: state rules on hand-built residues (pKa vs pH, metal rule, decision override).
- 5FQL panel: H335 is HID, D45/D46/D334 deprotonated; every heavy atom unchanged; no
  heavy atom without its hydrogens; determinism; a variant and the wild type share every
  state except where the record says they differ.
- Expected values written by hand from this spec's measurements (`reviewed: false`).

## Acceptance criteria
- 5FQL (modelled wild type, R468Q, R468W) at the manifest's pH: every residue has a
  recorded state and source; the Ca2+ ligands follow the metal rule; DDZ has HG1, HG2;
  no protonation work items remain; runs are byte-identical; CI green without network.

## Out of scope
Constant-pH MD, tautomers beyond His, Ca2+/DDZ/Cl- parameters (TASK-010), solvation.

## Open questions (defaults in bold)
1. **pH for the IDS study.** No default in code: the manifest must state it. For IDS, the
   variants' question is folding (misfolded variants are retained in the ER), which
   suggests **pH 7.2 (ER lumen, `[VERIFY]`)**; the alternative is 4.8 (lysosome, where the
   enzyme acts; about +20 e more charge and 43 ambiguous residues).
2. **Method.** **PROPKA 3.5.1 for pKa + OpenMM `Modeller` for hydrogens** (both
   pip-installable; the paper used PROPKA). Alternative: PDB2PQR (adds H-bond
   optimization but renames atoms and adds a dependency).
3. **Deposited hydrogens.** **Strip all and add every hydrogen with one method**, so every
   residue follows the same recorded states. Alternative: keep the 3,960 deposited ones
   and add only where missing (mixed provenance; their states came from refinement).
4. **Unreliable estimates.** **Flag titratable residues within 6.0 A (heavy atoms) of a
   metal ion or non-standard residue, apart from metal ligands, and default to the
   standard state**. At 5FQL (measured): R88 2.8, K135 3.5, H138 2.7, H229 4.3, K347 4.4,
   K479 6.0 A from DDZ A:84 / Ca2+ A:1551 (the ligands D45, D46, D334, H335 follow the
   metal rule). Alternative: trust PROPKA.
5. **Variants.** **PROPKA per system**, and the record lists every residue whose state
   differs from the wild type's (the mutation can shift nearby pKa; that is part of the
   comparison). Alternative: apply the wild type's states to every variant.

## [VERIFY] introduced
- ER lumen and lysosomal pH values (decision 1).

## Resolved decisions (maintainer, 2026-09-27)
Every default is accepted: (1) the manifest must state the pH; the IDS study uses pH 7.2
(ER lumen, `[VERIFY]`); (2) PROPKA 3.5.1 for pKa, OpenMM `Modeller` for hydrogens with
explicit variants; (3) all deposited hydrogens are stripped and every hydrogen is added
by one method; (4) metal ligands follow the metal rule, other titratable residues within
6.0 A of a metal ion or non-standard residue become findings defaulting to the standard
state; (5) PROPKA runs per system and the record lists states that differ from the wild
type's.
