# ADR-0008: Chemistry edits as residue mappings (FGly A:84 -> gem-diol DDZ)

- Status: accepted (TASK-008; maintainer, 2026-09-26: all spec defaults)
- Date: 2026-09-26

## Context

TASK-008 turns the deposited formylglycine sulfate ester (ALS A:84, 5FQL) into the
catalytic gem-diol. The spike showed the edit needs no new coordinates: ALS's CB already
carries both oxygens (CB-OG 1.41 A, CB-OS1 1.40 A, 112 degrees), so the gem-diol is ALS
without S, OS2, OS3, OS4. The CCD (read from RCSB) names it DDZ, 3,3-dihydroxy-L-alanine
(atoms N CA C O CB OG1 OG2, parent ALA); "FGH" in Demydchuk et al. 2017
(doi:10.1038/ncomms15786) is an abbreviation, and CCD `FGH` is an unrelated molecule.

## Decisions

1. **The edit is a prep operation driven by data** (TASK-008 decision 1). A residue
   mapping (`knowledge/residue_mappings.yaml`) deletes and renames atoms, never adds any,
   so ADR-0004 ("prep never invents coordinates") holds. `model_gem_diol` is an `apply`
   option; `revert_to_parent` and `model_gem_diol` share one operation.
2. **Mappings name the option they implement.** The schema gains `option` (one mapping
   per (component, option)) and `to_parent` (the pdbx_struct_mod_residue parent when the
   target is itself non-standard). `ALS -> DDZ` ships: OG -> OG1 (the paper's free
   geminal hydroxyl), OS1 -> OG2 (the Ca2+-bound one); DDZ's CB is not a stereocentre, so
   this is naming only (decision 3).
3. **A mapping keeps the annotations consistent** (this also fixes `revert_to_parent`,
   which renamed atoms but left them stale): links to deleted atoms are dropped, links to
   renamed atoms take the new residue and atom names, the modified-residue record takes
   the new component (or is removed when the target is standard), and the entity
   sequence takes the new component at the residue's position.
4. **Metal sites are noted, not re-decided** (decisions 2, 4). The action note lists each
   removed link; a removed metal link adds "metal coordination changed; the metal
   decision stands". At 5FQL the Ca2+ goes from six ligands to five (OS4 removed, OG2
   still at 2.32 A); nothing is placed in the vacated position.

## Consequences

- Every system after prep (wt/, wt_modelled/, variants) carries DDZ; the chemistry work
  item disappears from the 5FQL work order. DDZ, like ALS, has no Amber template: it is
  left out of relaxation and needs parameters (TASK-010) and hydroxyl hydrogens (TASK-009).
- The free aldehyde (`model_aldehyde`) stays `defer`: no CCD component verified.
