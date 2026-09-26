# ADR-0005: Variants v0.1: side-chain building, clash criterion, two-pass decisions

- Status: accepted (TASK-005)
- Date: 2026-09-26

## Context

TASK-005 builds matched variant systems (IDS 5FQL: WT vs R468Q vs R468W) from the
prepared wild type. A spike showed that PDBFixer 1.12.0 places a template side chain
without choosing a rotamer (its R468W ring passed 1.0 A from the Ser470 backbone). The
maintainer chose an in-house builder with cited data, hydrogens removed from the mutated
residue, a van der Waals overlap criterion, and a hard stop when every candidate clashes.
This ADR records how those decisions were implemented and where the data forced a
refinement.

## Decisions

1. **Builder.** `simprep.variants.build` places side-chain atoms from the kept N, CA,
   CB with ideal internal coordinates (NeRF, Parsons et al. 2005, doi:10.1002/jcc.20237)
   from `knowledge/side_chains.yaml` (GLN, TRP; schema `side_chains`). Candidates are
   the library rotamers plus the 180-degree flip of a planar amide (Lovell et al. 2000;
   Word et al. 1999: about 20 % of Asn/Gln amides need a flip). The geometry numbers are
   `[VERIFY]` (recalled, not re-read); a test rebuilds every deposited GLN and TRP of the
   panel from its own chi angles and compares with the deposited coordinates (GLN median
   RMSD 0.04 A, TRP 0.16 A), so a wrong number fails CI whatever its source.
2. **Clash criterion (refines TASK-005 decision 3).** A non-polar heavy-atom pair clashes
   when radii sum minus distance reaches 0.4 A, as decided. Applied literally, the
   criterion flags hydrogen bonds: 24 of 26 "clashes" found on deposited side chains
   were N/O pairs. Polar pairs (N/O) therefore clash only below 2.3 A, a limit taken from
   the panel itself (0.1 % of 1592 deposited polar contacts are shorter than 2.27 A);
   the test suite recomputes it. With heavy atoms only, 0.4 A still flags 3.8 % of
   deposited GLN/TRP/ARG side chains (C-O contacts of 2.7-2.8 A); recorded in the rule
   and left for the maintainer (see open questions in the PR). Acceptor-acceptor pairs
   are treated like any polar pair (no donor/acceptor typing in v0.1).
3. **Findings, not choices.** Each mutation site becomes a `variant_build` finding
   (rule `variant_build.rotamer_choice`, blocking) whose evidence lists every candidate
   (new evidence type `candidate`: chi, library frequency, clash count, closest contact).
   The recommended option is `choose_rotamer` with the fewest-clash, most frequent
   rotamer named in evidence, or `expert_review` when every candidate clashes (a rule
   condition, so the policy is data). `simprep variants` refuses a clashing rotamer even
   when a human chooses it (decision 4).
4. **Manifest.** `variants` (name, mutations in author numbering with from/to component
   ids, rationale, references) and `variant_snapshot` (findings plus the RFC 8785 hash of
   the variants they were computed for) are optional manifest fields. Variant decisions
   live in the same `decisions` list. `manifest status` counts variant findings; prep's
   gate does not (the wild type does not depend on them).
5. **Two passes.** Like `audit`, `simprep variants` first writes `manifest.json` with a
   fresh `variant_snapshot` and exits 3 when the snapshot is missing or stale; with a
   current snapshot and final, clash-free decisions it writes `wt/` (the plain prep
   output, byte-identical to `simprep prep`), one directory per variant,
   `variant_record.json` and `variant_report.md`.
6. **What a variant changes.** The mutated residue keeps N, CA, C, O, CB (and OXT) and
   loses every hydrogen and other side-chain atom; the entity sequence is updated; links
   to removed atoms are dropped. A test checks that WT and variant differ only there.
   Each variant's work order is the wild type's plus a `protonation` item.
7. **Accounting by atom name.** Atom records of a changed residue are matched by name:
   kept names are modified, disappearing names excluded, new names added. (TASK-004
   counted a changed residue's kept atoms as modified and could not express additions.)
8. **UniProt position** comes from the file's `_struct_ref_seq` / `_struct_ref` (new
   `SequenceReference` in the model), never assumed; 5FQL maps author 468 to P22304 468.
9. **Purity.** The decision-status logic moved from `simprep.manifest.status` to
   `simprep.decisions` (and `ManifestError` to `simprep.errors`) so the prep and variant
   cores import without jsonschema; the pure-core guard test now covers them.

## Consequences

- R468W cannot be built in v0.1: every rotamer (and a 10-degree chi grid) overlaps the
  kept atoms by at least 0.74 A. Building it needs a decision this task did not make
  (relaxing the surroundings, or accepting a clash with a recorded rationale and
  relaxing at the MD stage), and a tool for it would need its own ADR.
- A new mutation target needs side-chain data and nothing else; a new element near a
  site needs a radius in `variant_build.yaml` (unknown elements stop the build).
- Multi-chain entities: the output entity sequence is taken per chain; a mutation in one
  copy of a homo-oligomer is not yet handled specially (not exercised by the panel).
