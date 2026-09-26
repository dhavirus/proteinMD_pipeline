# TASK-003: Covalent-contact candidates the annotations miss

Read `CLAUDE.md` and ADR-0001 first. This task adds the rule agreed in the TASK-001 review
(1HZH): flag chemically plausible covalent attachments that the file does not annotate.

Status: **accepted** (maintainer, 2026-09-26: defaults accepted for every open
question; the resolved decisions are listed at the end of this file).

## Why

Detection today is annotation-driven: an attachment exists only if `struct_conn` says so.
In 1HZH (IgG1 b12) neither Fc glycan has an ASN–NAG record. Both glycans are therefore
reported as unrecognized groups **without** an attachment, and a preparation step would
simulate them as free sugars. The deposited geometry says otherwise:

| glycan | NAG 1 C1 to | distance |
|---|---|---|
| chain A glycan | ASN H314 ND2 | 2.64 Å |
| chain B glycan | ASN K314 ND2 | 2.45 Å |

(Measured on the committed `tests/panel/1HZH.cif.gz` during TASK-001.)

## What a naive distance check gets wrong (measured on the panel)

A survey of every unannotated non-polymer/branched ↔ polymer heavy-atom pair under 3.0 Å:

| entry | pairs < 3.0 Å | < 2.6 Å | what they are |
|---|---:|---:|---|
| 5FQL | 38 | 21 | almost all atoms *next to* an annotated bond (NAG O5 ↔ ASN ND2 at 2.3 Å across the annotated C1–ND2 link) |
| 3KS3 | 2 | 0 | glycerol hydrogen-bonding backbone O |
| 1HZH | 14 | 3 | the two real attachments above, plus a clash (NAG O7 ↔ VAL H277 CG1, 2.53 Å) |
| 6OIM | 11 | 0 | neighbours of the annotated MOV–Cys12 bond; GDP hydrogen bonds |
| 1FO8 | 9 | 5 | methylmercury altloc B (0.10) overlapping Phe316 altloc A (0.95): alternative states, never simultaneous |

So the rule must (1) skip residue pairs that already have an annotated covalent link,
(2) compare only altloc-compatible atoms (same altloc, or either has none), and (3) match
chemically plausible attachment patterns, not just short distances.

## Goal

A new rule family, `covalent_contacts`, whose detector reports each unannotated,
altloc-compatible atom pair that matches a declared attachment pattern within that
pattern's distance window, as a finding that says "this looks bonded but is not
annotated", with the geometry as evidence and a decision to make.

## Deliverables

### 1. Knowledge base: `knowledge/covalent_contacts.yaml`
- One rule, `covalent_contacts.attachment_candidate`, base severity **blocking**: the
  topology cannot be built until someone decides whether the bond exists.
- Matcher (new kind `attachment_pattern`, schema change + example documents):
  a list of patterns, each `{name, ligand: {comp_ids | "*", atom_names}, polymer:
  {res_names, atom_names}, max_distance_angstrom}`. v0.1 patterns, each exercised or
  explicitly negative-tested on the panel:
  - `n_glycosylation`: ligand NAG C1 ↔ ASN ND2 (exercised: 1HZH, and must stay silent on
    5FQL where the links are annotated).
  - `cysteine_adduct`: any non-polymer carbon ↔ CYS SG (negative control: 6OIM, where the
    MOV–Cys12 bond is annotated and its neighbours must not fire).
  - Further patterns (O-glycosylation SER OG / THR OG1, lysine adducts) only when a panel
    structure exercises them (CLAUDE.md), so not in v0.1 (decision 1).
- Distance windows (decision 2): one upper limit per pattern, 3.0 Å for both v0.1
  patterns (above the 2.64 Å measured in 1HZH, below typical non-bonded C–N / C–S
  contacts). It is a working threshold recorded in the rule file with that rationale,
  not a literature value; no van der Waals radius table in v0.1.
- Options: `add_link` (record the bond; preparation builds it), `treat_noncovalent`
  (keep the group as a separate molecule; the short contact stays a modelling issue),
  `exclude_group`, `expert_review`. Recommended: `add_link` when the pattern is
  `n_glycosylation` **and** the ASN is in an N-X-S/T sequon (checked from the entity
  sequence, X ≠ P); otherwise `expert_review`. Base severity blocking (decision 5),
  escalated by region like every other finding (it cannot rise further).

### 2. Detector: `src/simprep/detectors/covalent_contacts.py`
- Pure function, registered as `covalent_contacts`, over the existing `Structure` model.
- Candidates: for each (non-polymer or branched residue, polymer residue) pair **without**
  any annotated covalent `struct_conn` between them, each altloc-compatible atom pair
  matching a pattern within its window.
- Evidence: the distance (with both atoms and altlocs), the pattern name, the expected
  bond, the sequon check result (`flag`) for `n_glycosylation`, the partner group's
  unrecognized finding id, and a `source_record` noting that `struct_conn` has no such
  link.
- Needs the polymer sequence per chain for the sequon check. Add `entity_poly_seq` (or the
  gemmi entity full sequence) to the `Structure` model; record its presence in
  `annotation_categories`.
- Candidates are not claims: the glycan stays an `unrecognized` finding. The candidate
  finding's anchors are both residues.
- Neighbour search: brute force is fine at panel size (tens of ligand atoms × thousands
  of polymer atoms); keep it pure Python (ADR-0002 guard test).

### 3. Expected findings (derived independently of the detector)
- 1HZH: two candidates, NAG A1 ↔ ASN H314 and NAG B1 ↔ ASN K314, both `n_glycosylation`,
  both in a sequon. Independent derivation: the sequon from the file's `entity_poly_seq`
  (derive script) plus the measured C1–ND2 geometry are accepted as the independent
  evidence (decision 3). The statement that ASN 314 is the IgG1 Fc Asn297 (EU numbering)
  site carries `[VERIFY]` in the expected-findings file until a source is read; the
  1HZH paper's abstract (Saphire et al. 2001, doi:10.1126/science.1061692) does not
  state it.
- 5FQL, 3KS3, 6OIM, 1FO8: **no** candidates (negative controls, each one a regression
  test for one of the false-positive classes in the survey above).
- `derive_expected.py` gains the sequon derivation; the candidate list itself is written
  by hand into `curated_checks` with its source, not computed by a copy of the detector.

### 4. Front end
- `FAMILY_LABELS` gains `covalent_contacts`; the 3D view shows both atoms of the pair.
  No other UI change; regenerate the demo only if 5FQL output changes (it must not).

### 5. Tests
- Unit tests on hand-built structures: annotated pair skipped; altloc-incompatible pair
  skipped; outside the window skipped; pattern mismatch skipped; sequon on/off changes the
  recommendation; proline in the sequon's X position disables it.
- Panel regression as in §3. Schema examples for the new matcher kind.

## Acceptance criteria
- 1HZH yields exactly the two `n_glycosylation` candidates, blocking, recommended
  `add_link`, with the sequon flag true.
- 5FQL, 3KS3, 6OIM and 1FO8 yield none.
- The new rule file validates; the detector passes the pure-core guard test; CI is green
  without network.

## Out of scope
Building the bond (a preparation-stage task), generic steric-clash reporting (decision 4),
patterns no panel structure exercises, modelling missing glycan residues.

## Resolved decisions (maintainer, 2026-09-26)
1. **Patterns in v0.1**: `n_glycosylation` (exercised by 1HZH) and `cysteine_adduct`
   (negative-tested on 6OIM) only. Other chemistries wait for a panel structure.
2. **Distance windows**: a single upper limit per pattern (3.0 Å for both), a documented
   working threshold, not a literature value.
3. **1HZH evidence**: the sequon from the file's own sequence plus the measured geometry
   count as independent evidence; the Asn297 literature link stays `[VERIFY]` and the
   expected file stays `reviewed: false` until a source is read.
4. **Steric clashes** (e.g. 1HZH NAG O7 ↔ VAL H277 CG1, 2.53 Å): out of scope; a
   separate `clashes` family can come later if needed.
5. **Severity**: blocking.
