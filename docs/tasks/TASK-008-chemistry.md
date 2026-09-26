# TASK-008: Chemistry v0.1: formylglycine A:84 as the gem-diol (DDZ)

Read `CLAUDE.md`, ADR-0004 (prep never invents coordinates; revert mappings are data),
ADR-0007 (modelling stage) first. After TASK-007 the 5FQL work order still holds:

| stage | item (5FQL prep fixture) |
|---|---|
| chemistry | `nonstandard_residues/A:84`: FGly as the gem-diol (`model_gem_diol`) |
| parameterization | Ca2+ A:1551; Cl- A:1567-1569 |
| protonation | modelled loop A:444-453, moved flanks, mutated residues |

Status: **draft** (open questions below, each with a default, for the maintainer).

## What the primary publication says (read in this session, PubMed PMC5472762)

Demydchuk et al. 2017, Nat Commun 8:15786, doi:10.1038/ncomms15786:

- "FGly functions catalytically as an aldehyde hydrate (FGly aldehyde-hydrate, FGH) with
  two geminal hydroxyl groups; here we observe what we have interpreted as a covalent
  sulfate-ester intermediate (FGly sulfate, FGS) in the IDS active site." An alternative
  phosphate-ester interpretation "cannot be excluded", but a sulfate ester is "more likely".
- The Ca2+ is "coordinated by one nitrogen and five oxygen atoms with approximate
  octahedral geometry: one oxygen each from the side chains of D45, D46 and D334; two
  sulfate oxygens from FGS84 and one nitrogen from the side chain of H335".
- H138 and R88 contact "the free geminal hydroxyl of FGS"; H229 contacts sulfate oxygen OS3.
- The mechanism: FGH84 attacks the substrate sulfur, giving the covalent FGS adduct, then
  sulfate is eliminated and the aldehyde rehydrated; elimination is expected to be
  efficient only at lysosomal pH (about 4.8).
- For their own docking the authors did exactly this edit: "the sulfated FGS84 residue was
  replaced with FGH (Protein Data Bank code DDZ: 3,3-dihydroxy-alanine)", with restraints
  to "preserve Ca2+ coordination geometry".

"FGH" is the paper's abbreviation, **not** a chemical component ID: CCD `FGH` is an
unrelated non-polymer (C27 H35 N3 O2, read in this session). The current option label
"(gem-diol, FGH)" in `knowledge/nonstandard_residues.yaml` and the roadmap in
`docs/HANDOFF.md` must say DDZ.

## What the CCD says (read from RCSB in this session)

| component | name | formula | parent | heavy atoms |
|---|---|---|---|---|
| ALS | (3S)-3-(sulfooxy)-L-serine | C3 H7 N O7 S | ALA | N CA C O CB OG OS1 S OS2 OS3 OS4 (+OXT) |
| DDZ | 3,3-DIHYDROXY L-ALANINE | C3 H7 N O4 | ALA | N CA C O CB OG1 OG2 (+OXT) |

Both are released, L-peptide linking. This clears the `[VERIFY]` on DDZ (and confirms ALS
beyond the 5FQL chem_comp record). The CB of DDZ carries two OH groups and is not a
stereocentre. No CCD component for the free FGly aldehyde was looked up; it stays
`[VERIFY]` and out of scope.

## What the file says (measured in this session)

| fact | value |
|---|---|
| ALS A:84 bonds | CB-OG 1.41 A, CB-OS1 1.40 A, OS1-S 1.50 A (the ester oxygen), OG-CB-OS1 112 degrees |
| Ca2+ A:1551 ligands | D46 OD1 2.29, H335 NE2 2.29, ALS OS1 2.32, D334 OD2 2.32, D45 OD1 2.37, ALS OS4 2.47 A (and S 2.90, D334 OD1 3.32) |
| struct_conn at A:84 | covale1 (V83 C - N), covale2 (C - A85 N), metalc3 (OS1 - Ca), metalc4 (S - Ca), metalc5 (OS4 - Ca) |
| hydrogens on A:84 | none deposited |

## Spike: the edit needs no new coordinates

The gem-diol is ALS without its sulfate: delete S, OS2, OS3, OS4; OG becomes one hydroxyl
and the ester oxygen OS1 the other. CB already has gem-diol geometry (1.40-1.41 A, 112
degrees). The Ca2+ keeps five ligands (OS1, now a hydroxyl, stays at 2.32 A) and loses
OS4: coordination 6 -> 5. That is a delete-and-rename, the same kind of edit as prep's
`revert_to_parent` with a data mapping (ADR-0004 decision 3), not a modelling step.

Found on the way: prep's rename path (`_reverted` in `prep/apply.py`) renames atoms but
does not update `struct_conn` partners, the modified-residue record, or the entity
sequence. No mapping ships yet, so it was never exercised; this task needs all three.

## Scope in one sentence

Turn `model_gem_diol` into a recorded, data-driven edit of the deposited residue (ALS ->
DDZ), with links, annotations and the metal site kept consistent, so the parameterization
stage receives a gem-diol FGly.

## Deliverables

### 1. Knowledge base
- `knowledge/residue_mappings.yaml`: a mapping `ALS -> DDZ` (delete S, OS2, OS3, OS4;
  rename per decision 3), with the paper and the CCD as references. The mapping schema
  gains what it needs (target parent, which links to drop and which to rename follow
  from the atoms).
- `nonstandard_residues.formylglycine`: `model_gem_diol` becomes `apply` (decision 1); its
  label says DDZ, not FGH; the DDZ `[VERIFY]` is removed (CCD read). `model_aldehyde`
  stays `defer`.

### 2. Prep (or a stage, decision 1)
- Apply the mapping: residue renamed ALS -> DDZ, atoms deleted and renamed.
- `struct_conn`: links to deleted atoms are dropped (metalc4 S-Ca, metalc5 OS4-Ca); links
  to renamed atoms follow the new names (metalc3 becomes OG2-Ca; covale1/2 name DDZ).
- `pdbx_struct_mod_residue`: the record becomes DDZ (parent ALA, per the CCD).
- Entity sequence: ALS -> DDZ at A:84.
- Accounting (atoms matched by name, ADR-0004): 1 polymer residue modified; atoms 5
  modified (N CA C O CB), 6 excluded (OG OS1 S OS2 OS3 OS4), 2 added (OG1 OG2); links 2
  excluded (metalc4, metalc5), 3 modified (covale1, covale2, metalc3).
- The same machinery fixes `revert_to_parent` renames (links, mod record, sequence).

### 3. Metal site
- The record notes that the Ca2+ coordination changed (6 -> 5 ligands within the metals
  rule's shell, the sulfate OS4 removed), for the parameterization stage (decision 2).

### 4. Tests (offline)
- Unit tests on hand-built residues: delete/rename, link drop/rename, mod record, sequence.
- 5FQL panel: DDZ at A:84 with N CA C O CB OG1 OG2; no S; Ca2+ links as above; accounting;
  written mmCIF re-read; `simprep variants` and `simprep model` carry DDZ.
- Expected values written by hand from this spec's measurements (`reviewed: false`).

## Acceptance criteria
- 5FQL with `model_gem_diol` for A:84: DDZ at A:84 with exactly the heavy atoms N, CA, C,
  O, CB, OG1, OG2 at the deposited coordinates of N, CA, C, O, CB, OG, OS1 (decision 3);
  no sulfate atoms; one Ca2+ link to A:84 (OG2, 2.32 A), none to S or OS4; the
  modified-residue record and the entity sequence say DDZ; counts balance.
- Nothing else in the file changes; runs are byte-identical.
- CI green without network.

## Out of scope
The free aldehyde, the phosphate-ester reading, hydrogens on the hydroxyls (TASK-009
protonation), DDZ / Ca2+ parameters (TASK-010), active-site water placement.

## Open questions (defaults in bold)
1. **Where the edit runs.** **Prep, as an `apply` operation driven by the mapping** (it
   deletes and renames, never invents coordinates, like `revert_to_parent`; so `wt/`,
   `wt_modelled/` and every variant carry DDZ). Alternative: a separate chemistry stage
   like `simprep model` (more machinery, no coordinates to build).
2. **The vacated Ca2+ site.** **Record it and add nothing** (5 ligands; water from
   solvation fills it, and the parameterization stage sees the note). Alternative: place
   a water at the OS4 position now (an invented coordinate; would need its own decision
   and marking, like the loop).
3. **Which oxygen gets which name.** **OG -> OG1 (the "free geminal hydroxyl" of the paper,
   contacting H138 / R88), OS1 -> OG2 (the Ca2+-bound one).** DDZ's CB is not chiral, so
   this is naming only; recorded in the mapping.
4. **Metal finding.** The manifest's `metals/A:1551` decision (`nonbonded_12_6_4`) was taken
   on a 6-coordinate site. **Keep the decision and note the change in the prep record**
   (no re-audit). Alternative: re-run the metals detector on the edited system and ask
   for a new decision.

## [VERIFY] introduced
- None new. Removed: DDZ identity (CCD read). Kept: FGly aldehyde component ID; ALS as
  FGly sulfate in other entries.
