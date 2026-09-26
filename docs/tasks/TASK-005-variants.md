# TASK-005: Variant building v0.1: WT vs R468Q vs R468W

Read `CLAUDE.md`, ADR-0001 and ADR-0004 first. TASK-004 turns a decided manifest into a
prepared wild-type (WT) system and a work order. The comparative study CLAUDE.md names
(IDS, 5FQL: WT vs R468Q vs R468W) needs matched variant systems: identical to the
prepared WT in every record except the mutated residue, with the new side chain built,
checked and documented.

Status: **accepted** (maintainer, 2026-09-26: defaults accepted for every open
question; the resolved decisions are listed at the end of this file).

## What the committed 5FQL file says (measured in this session)

| fact | value | source |
|---|---|---|
| author numbering vs UniProt P22304 | identical: author 26-550 = UniProt 26-550 | `_struct_ref_seq` (db_align_beg 26 = pdbx_auth_seq_align_beg 26) |
| residue A:468 | ARG, all 11 heavy atoms modelled, occupancy 1.0, no altloc | `_atom_site` |
| A:468 B-factors | 59.5-64.1 A^2 (median CA B of the chain: 61.3) | `_atom_site` |
| A:468 hydrogens | 13 (the file carries 4094 riding H in total) | `_atom_site` |
| A:468 guanidinium contacts (< 3.5 A) | backbone O of V262 (3.00), A263 (3.18), Y348 (2.92), S477 (2.84); water 2127 (3.24) | geometry |
| nearest other findings | ALS A:84 (FGly) 12.0 A; Ca2+ A:1551 14.2 A; glycans E, F 14.7-14.8 A | geometry |

So the site is fully modelled, unambiguous, and away from every altloc and gap finding.
It is at the escalation distances of the default thresholds (6 / 12 A) only relative to
ALS A:84 (12.0 A); whether a region should include it is a study decision (CLAUDE.md).

## Spike: an off-the-shelf builder is not enough (measured in this session)

PDBFixer 1.12.0 / OpenMM 8.6.1 (both on PyPI) were run on the TASK-004 5FQL fixture
output (`applyMutations`, missing atoms restricted to residue 468):

| variant | time (CPU) | other atoms moved / lost | worst new heavy-atom contact |
|---|---|---|---|
| R468Q | 0.6 s | 0 / 0 | Gln OE1 - Val262 O 2.49 A |
| R468W | 3.9 s | 0 / 0 | Trp CH2 - Ser470 C **1.02 A** (ring through the backbone) |

PDBFixer places a template side chain; it does not search rotamers, and here it
produced an impossible Trp. It also leaves the new residue without hydrogens while the
rest of 5FQL keeps its deposited ones. So v0.1 needs (a) rotamer candidates, (b) a clash
check that turns bad placements into findings, and (c) an explicit hydrogen policy (decision 2).

## Scope in one sentence

Build each variant from the prepared WT system by replacing one residue's side chain
with rotamer candidates from a cited library, report every candidate's contacts as
evidence, let the human choose (recommended default: fewest clashes), and write matched
WT and variant systems with full accounting. No minimization, no MD, no repacking of
neighbours (that would make WT and variant differ in more than the mutation).

## Deliverables

### 1. Manifest: `variants` (schema change, examples)
- `variants: [{name, mutations: [{chain, seq_num, ins_code, from, to}], rationale,
  references}]`; `name` defaults to the conventional form (e.g. `R468Q`).
- Validated against the structure before anything is built, each failure listed:
  `from` matches the modelled residue **and** the entity sequence; the residue is
  observed with N, CA, C and CB; no altlocs left after prep; not excluded, not in a gap,
  not covered by a finding whose decision defers it (e.g. a non-standard residue);
  `to` is a residue the builder has data for.
- Mapping to UniProt numbering is **reported**, from the file's `_struct_ref_seq`
  (5FQL: identity), never assumed. A file without that record says so.

### 2. Knowledge base: side-chain data as rules-are-data
- `knowledge/side_chains.yaml` (new schema): per residue type, the atoms, ideal internal
  coordinates (bond lengths, angles, dihedral references) and the rotamers (chi means +
  labels) with the source cited (decision 1). v0.1 ships **GLN and TRP** only (CLAUDE.md: no
  feature a panel structure or IDS does not exercise); other targets are refused with
  "add side-chain data". Numbers not verified in this session carry `[VERIFY]`.
- `knowledge/variant_build.yaml`: rule family `variant_build`, e.g.
  `variant_build.steric_clash` (base severity blocking) with the clash criterion as data
  (decision 3), options `choose_rotamer` (apply; parameter `rotamer`), `expert_review`
  (unresolved), and a recommendation condition "fewest clashes, then library
  probability" (decision 4).

### 3. `src/simprep/variants/` (pure core, thin I/O)
- `build.py`: pure; place side-chain atoms from the WT backbone (N, CA, C, CB kept) with
  the ideal geometry and a chosen chi set (NeRF). No gemmi, no external builder.
- `check.py`: pure; heavy-atom contacts of each candidate with every other atom kept in
  the system (altloc-aware, as in TASK-003), emitted as `variant_build` findings with the
  atoms and distances as evidence.
- `apply.py`: pure; the prepared WT `Structure` + a chosen candidate -> variant
  `Structure`: side chain replaced, hydrogens per Q2, the entity full sequence updated at
  that position, `struct_conn` records touching the replaced atoms dropped and counted.
- Outputs reuse `prep/write.py` and `prep/record.py` accounting: the first use of the
  `added` column for atom records (new side-chain atoms); removed WT side-chain atoms
  are `excluded`.

### 4. CLI and record
- `simprep variants STRUCTURE --manifest M --out DIR` runs prep (WT) and then builds
  every variant into `DIR/<name>/`, plus `DIR/wt/`. Two passes, as in audit:
  1. without rotamer decisions it writes `variant_findings.json` (candidates + clashes)
     and stops with exit 3 if any blocking `variant_build` finding is undecided;
  2. with decisions recorded in the manifest it writes the systems.
- `variant_record.json` (new schema) per variant: the mutation(s), UniProt mapping,
  chosen rotamer and its contacts, file hashes, accounting against WT, provenance.
- A test asserts WT and each variant differ **only** at the mutated residue (and its
  sequence entry and dropped connections).

### 5. Front end (small)
- The review page shows `variant_build` findings like any other family (FAMILY_LABELS,
  3D view of the candidate's clashing atoms), and edits the `variants` list (decision 5).

### 6. Colab
- `notebooks/02_variants.ipynb` + `modules/variants_module.py`, same conventions as
  TASK-004 (one run directory holds `wt/` and every variant). (decision 6)

### 7. Tests (offline)
- Builder geometry checked against **deposited** data, not against itself: rebuild every
  GLN and TRP in the panel from its own observed chi angles and backbone; heavy-atom RMSD
  to the deposited side chain below a stated tolerance (the tolerance is a test constant
  with a rationale, Q1).
- Unit tests: NeRF placement on hand-built backbones; clash detection (altloc-aware,
  bonded neighbours excluded); every validation refusal; recommendation ordering.
- 5FQL R468Q and R468W: candidates listed, clash findings for the spike's placements
  reproduced as a negative control (a Trp through Ser470 must be flagged), the chosen
  rotamer written, accounting equations, re-read counts, WT/variant diff limited to
  A:468, determinism.
- Expected findings for R468Q/R468W are written by hand from geometry measured
  independently of the builder, `reviewed: false`.

## Acceptance criteria
- For 5FQL with the TASK-004 fixture decisions and variants R468Q and R468W:
  `simprep variants` first stops with the rotamer choices as blocking findings, then,
  with decisions, writes `wt/`, `R468Q/`, `R468W/` whose re-read counts match the
  records and which differ from WT only at A:468.
- Every candidate's closest contacts are recorded; no candidate is silently chosen.
- A mismatched `from` (e.g. `K468X` on 5FQL) and a target without side-chain data are
  refused with actionable messages.
- CI green without network.

## Out of scope
Minimization or relaxation of the variant, repacking neighbours, protonation / pKa,
alchemical (hybrid) topologies for FEP, glycine/proline or non-standard targets,
insertions/deletions, multi-chain assembly, MD. Any claim about what R468Q or R468W do
to IDS (stability, activity, disease severity) belongs to the study, with its own
citations, not to this task.

## Resolved decisions (maintainer, 2026-09-26)
1. **Builder and data**: an in-house NeRF builder (pure Python) places side chains from
   ideal internal coordinates and rotamer chi values stored in
   `knowledge/side_chains.yaml`, GLN and TRP only in v0.1. Candidate sources are the
   penultimate rotamer library (Lovell et al. 2000) or the Dunbrack
   backbone-independent library for chi values, and Engh & Huber for ideal geometry;
   every number stays `[VERIFY]` until its source is read in a session and cited with a
   DOI. PDBFixer and external packers (FASPR, SCWRL4) are not used. The choice is
   recorded in an ADR (TASK-004 decision 2).
2. **Hydrogens**: the mutated residue of each variant carries no hydrogens (the WT side
   chain's hydrogens are excluded with it), and each variant's work order gains
   "protonation: A:468 has no hydrogens". WT keeps its deposited hydrogens; the
   protonation stage later treats every system alike.
3. **Clash criterion**: a heavy-atom pair is a clash when its distance is below the sum
   of the two van der Waals radii minus 0.4 A; radii from a cited table, `[VERIFY]`
   until read. Bonded neighbours (1-2, 1-3) and the residue's own backbone are never
   clashes.
4. **Recommendation**: fewest clashes, ties broken by library probability. If every
   candidate clashes, the recommendation is `expert_review` and the finding stays
   blocking: a hard stop in v0.1 (no acceptance with a rationale, no relaxation here).
5. **Front end**: the review page gets a small variant-list editor; mutations are typed
   as `A:468 R>Q` and validated against the loaded structure.
6. **Colab**: a separate `notebooks/02_variants.ipynb` runs WT and all variants in one
   run directory.
7. **Naming**: variant names use the author numbering (`R468Q`); the UniProt position
   from `_struct_ref_seq` is recorded alongside (identical for 5FQL).
