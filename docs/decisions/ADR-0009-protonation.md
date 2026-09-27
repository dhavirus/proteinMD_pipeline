# ADR-0009: Protonation v0.1: PROPKA estimates, recorded states, OpenMM hydrogens

- Status: accepted (TASK-009; maintainer, 2026-09-27: all spec defaults)
- Date: 2026-09-27

## Context

After TASK-008 the work order holds protonation items (modelled loop, moved flanks,
mutated residues, DDZ hydroxyls), and the deposited structure carries 3,960 riding
hydrogens from refinement. The spike (PROPKA 3.5.1 on the modelled 5FQL wild type) showed
that the pH moves the net charge by about 20 e between pH 5 and 7, that 9-43 residues sit
within one pH unit of their pKa depending on the pH, and that PROPKA, which ignores the
Ca2+ and DDZ, gives implausible values in the active site (D45 -8.1, H335 -0.7, K347 6.7).

## Decisions

1. **The pH is the study's.** The manifest's `protonation` states `ph` and `ph_rationale`;
   `simprep protonate` refuses without them, and `simprep variants` protonates only when
   they are present. Nothing defaults the pH. The method (PROPKA version, ambiguity
   window 1.0 pH unit, 6.0 A distance, seed, Reference platform) comes from
   `knowledge/protonation.yaml` unless the manifest overrides it.
2. **Tools.** PROPKA 3.5.1 (Olsson et al. 2011, doi:10.1021/ct100578z) as a library, on
   heavy atoms without waters; OpenMM `Modeller.addHydrogens` with an explicit variant per
   titratable residue (seeded, Reference platform). Both are in the `relax` extra. Neutral
   His is left to Modeller's hydrogen-bond tautomer choice; its pH argument is fixed above
   6.5 for that reason and affects no other residue. The record names the tautomer chosen.
3. **All hydrogens are replaced** (decision 3): deposited hydrogens are stripped and every
   hydrogen is added by one method. Heavy atoms keep their coordinates; each hydrogen
   takes its parent's occupancy and B-factor (so modelled-loop hydrogens are 0.00).
4. **Rules before estimates.** A disulfide cysteine is CYX; a residue whose side chain is
   linked to a metal in `struct_conn` is not protonated there (Asp/Glu/Cys deprotonated,
   His HID when NE2 is bound, HIE when ND1 is). Then, per residue: a decision on its
   finding, else the finding's recommended option, else the estimate at the pH.
5. **Findings** (family `protonation`, warn, per system, id
   `protonation/<system>/<residue>`): `unreliable_estimate` for titratable residues within
   6.0 A of a metal ion (elements from `metals.yaml`) or a non-standard polymer residue,
   default the model pKa's state; otherwise `ambiguous_state` within the window, default
   the estimate's state. Options: `predicted_state`, `standard_state`, `specific_state`
   (parameters.state). Undecided findings take the default and the record says so. Arg
   and Tyr have no alternative variant in the force field: their pKa is recorded, they
   get no finding (so R88, 2.8 A from DDZ, is not flagged).
6. **Per system** (decision 5): each system gets its own PROPKA run. `simprep variants`
   protonates the wild type it built the variants on, the relaxed wild types and every
   variant into `<name>_protonated/`; the record lists every state that differs from the
   wild type's. `simprep protonate` does the wild type alone (prepared, or modelled).
7. **Hard stops.** A polymer residue that is neither standard nor defined in
   `residue_definitions`, or a multi-atom non-polymer residue, stops the stage (no
   hydrogen definitions); a titratable residue without an estimate stops it too.

## Consequences

- 5FQL at pH 7.2 (wild type): 102 titratable residues; D45, D46, D334 deprotonated and H335
  HID by the metal rule; four CYX, one free CYS; 13 findings (K135, H138, H229, K347, K479
  unreliable; eight ambiguous); DDZ gets H, HA, HB, HG1, HG2 from the CCD names.
- Protonation takes seconds per system (PROPKA about 7 s, hydrogens about 6 s).
- Parameters for DDZ and Ca2+ (TASK-010) are still needed before a full topology.
