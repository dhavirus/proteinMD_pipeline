# TASK-007: Modelling v0.1: the 444-453 linker and the chain ends of IDS

Read `CLAUDE.md`, ADR-0004 (prep never invents coordinates), ADR-0006 (OpenMM relaxation)
first. After TASK-006 the 5FQL work order still holds, in stage order:

| stage | item (5FQL prep fixture) |
|---|---|
| modelling | `missing_residues/A:26-33`: unobserved N-terminal residues (8) |
| modelling | `missing_residues/A:444-453`: internal gap (10 residues) |
| chemistry | `nonstandard_residues/A:84`: FGly as the gem-diol (FGH) |
| parameterization | Ca2+ A:1551; Cl- A:1567-1569 |

This task is the modelling stage. It is the first task that adds atoms with no experimental
support, so it has to say where they come from, mark them, and use the same model in the
wild type and every variant.

Status: **accepted** (maintainer, 2026-09-26: defaults accepted for every open
question; the resolved decisions are listed at the end of this file).

## What the primary publication says (read in this session)

Demydchuk et al. 2017, Nat Commun 8:15786, doi:10.1038/ncomms15786 (PubMed 28593992,
full text PMC5472762), the 5FQL paper:

- "The first 33 residues, encoding the N-terminal signal peptide and propeptide, are
  cleaved during enzyme secretion and therefore are not present in the structure."
  So A:26-33 are **not** missing from the mature enzyme: Thr34 is its real N-terminus.
- "Lack of interpretable electron density prevented the modelling of one loop region
  (residues 444-453)." It links SD1 (residues 34-443, the "heavy" 42 kDa chain) to
  SD2 (455-550, the "light" 14 kDa chain), the fragments of lysosomal proteolytic
  processing. It is "consistent with a disordered region that is susceptible to
  proteolytic cleavage", but "in our structure this loop, although disordered, is
  likely to remain intact as the protein is purified from conditioned medium and has
  not trafficked to lysosomes".
- The paper also notes that misfolded proteins are retained by ER quality control, where
  the precursor is still one chain. That favours the intact linker for a folding study,
  but it is the study's call (decision 1).

## What the file says (measured in this session)

| fact | value |
|---|---|
| gap sequence | A:444-453 FRDLEEDPYL (Phe-Arg-Asp-Leu-Glu-Glu-Asp-Pro-Tyr-Leu) |
| span to bridge | CA443 - CA454 = 17.4 A, over 11 peptide steps (up to about 42 A extended) |
| flank disorder | CA B-factors 108-137 A^2 (440-443) and 73-150 A^2 (454-457); median chain CA B 61 A^2 |
| distance to the study site | 32-36 A from R468 (heavy atoms); about 30 A from Ca2+ A:1551 and ALS A:84 |
| C-terminus | A:550 observed; the entity sequence ends there (nothing missing) |

## Spike: PDBFixer 1.12.0 loop building (this session)

Run on the prepared 5FQL WT (TASK-004 fixture), only the internal gap requested:

- 10 residues added in 18 s; no pre-existing atom moved; peptide bonds 1.29-1.43 A.
- **Geometry is not usable as is.** CA443 - CA444 = 2.61 A (a trans peptide gives about
  3.8 A), CA444 - CA445 = 4.46 A, and several phi/psi values are odd (for example Phe444
  phi +140, Asp450 phi +178). Two contacts with Pro454 are flagged by the TASK-005
  criterion.
- **Not reproducible.** Two runs with Python's `random` seeded differ by up to 6.7 A.
  The randomness comes from elsewhere (numpy or the platform); not investigated yet.

## Network (this session)

AlphaFold DB (`alphafold.ebi.ac.uk`) and UniProt (`rest.uniprot.org`) could not be
reached from this environment. A predicted model or UniProt feature annotations can only
be used if the maintainer supplies the file (as with the panel structures).

## Scope in one sentence

Turn `model_loop` decisions into modelled residues. They come from a recorded, deterministic
source and are marked as modelled everywhere. They are built once on the prepared wild
type and shared by every variant. The chain-end decisions (`truncate`, `cap_termini`)
get an exact, recorded meaning for the topology stage.

## Deliverables

### 1. Knowledge base and manifest
- `knowledge/modelling.yaml` (new schema) holds the defaults for loop modelling:
  - method;
  - junction geometry checks (C-N, CA-CA tolerances);
  - relaxation of the loop (mobile backbone and side chains of the loop and its two flanking residues);
  - seeds.
- Manifest `modelling` (optional): the protocol used, copied from the defaults and
  recorded, like `relaxation`.
- `missing_residues` options:
  - `model_loop` becomes `apply` (it was `defer`). It needs this ADR (TASK-004
    decision 2).
  - `truncate` stays `record`, but its meaning is fixed: the chain ends there, with
    charged termini (NH3+ / COO-), listed for the topology stage.
  - `cap_termini` stays `defer` until a study needs caps (decision 4).

### 2. `src/simprep/model/` (the loop builder at the edge, geometry checks pure)
- The loop builder takes the prepared structure, the gap (from the unobserved-residue
  annotations and the entity sequence) and the protocol, and returns the structure with
  the loop residues added.
- **Junction checks** (pure), run before the model is accepted:
  - C-N peptide bonds and CA-CA distances within tolerance;
  - no heavy-atom clash with the rest (TASK-005 criterion, flag only);
  - the loop residues' phi/psi are reported.
- **Marking**:
  - modelled atoms are listed in the record;
  - they are written with a recorded occupancy convention (decision 3);
  - the residues are removed from the `pdbx_unobs_or_zero_occ_residues` annotation, and
    the accounting counts them as `added`.
- **One model for every system**: the loop is built once on the prepared wild type.
  Variants and relaxation build on the modelled wild type, so the linker is identical in
  WT, R468Q and R468W.

### 3. Where it runs (decision 5)
- The default is a library stage used in-process by `simprep variants`, plus a
  `simprep model` command for the wild type alone:
  prep (select) -> model (add) -> variants (build) -> relax.
- Outputs: `wt_modelled/` (or `wt/` with a modelled flag); the record lists every
  modelled residue and its source.

### 4. Findings
- A new `modelling` family:
  - `modelling.junction_geometry` (blocking if a junction fails the checks);
  - `modelling.loop_contacts` (warn), with atoms and distances.
- The loop is marked "no experimental support" in the report, with a note that
  analyses (contact maps, native-structure models) should treat it as modelled.

### 5. Tests (offline)
- The builder puts the right sequence in the gap.
- The junctions pass the checks.
- Pre-existing atoms outside the loop's relaxation shell are unmoved.
- Two runs give identical files.
- Accounting: 10 polymer residues added, and the unobserved annotation shrinks by 10.
- WT and R468W share the loop coordinates exactly.
- Expected values are written by hand from this spec's measurements (`reviewed: false`).

## Acceptance criteria
- 5FQL with `model_loop` for A:444-453 and `truncate` for A:26-33:
  - 10 residues FRDLEEDPYL appear between 443 and 454;
  - every junction C-N is 1.33 +/- 0.05 A and every CA-CA 3.8 +/- 0.15 A (trans; widened
    from 0.1 by the maintainer on 2026-09-27, see the ADR-0007 amendment);
  - the output is deterministic;
  - the residues are marked as modelled;
  - the same loop appears in WT, R468Q and R468W.
- A:26-33 decided `truncate` produces no atoms, and the topology note "charged N-terminus
  at Thr34" is recorded.
- CI green without network.

## Out of scope
Loop conformational sampling or ensembles, the processed two-chain form as a separate
system (unless decision 1 picks it), caps, the FGly chemistry edit, parameterization,
protonation, MD.

## Resolved decisions (maintainer, 2026-09-26)
1. **System form**: the intact single chain, with the 444-453 linker modelled. This is the
   crystallized, secreted precursor, and the form that folds in the ER (Demydchuk et al.
   2017). The lysosomal two-chain form is not built in v0.1; building it later needs a
   source for the processing sites (`[VERIFY]`).
2. **Loop source**: PDBFixer 1.12.0 places the residues, then an OpenMM minimization
   (Reference platform) runs with the loop backbone and side chains mobile and the
   flanking residues restrained. Every random source is seeded. The first
   implementation step is to find and seed PDBFixer's randomness (the spike differed by
   6.7 A between runs); if byte-identical output cannot be achieved, the task stops and
   reports instead of shipping a non-deterministic builder. MODELLER, an AlphaFold DB
   graft and in-house loop closure are not used. Recorded in an ADR.
3. **Marking modelled atoms**: occupancy 0.00 in the written files. B-factors are set to
   a recorded value. The record lists every modelled residue and atom with its source
   and protocol.
4. **Caps**: none in v0.1. `cap_termini` stays `defer`; IDS needs no caps once A:26-33
   is decided `truncate` (the propeptide is cleaved; Demydchuk et al. 2017).
5. **Pipeline placement**: modelling is its own stage between prep and variants:
   prep (select) -> model (add) -> variants (build) -> relax. It is used in-process by
   `simprep variants` and exposed as `simprep model` for the wild type alone. `simprep
   prep` still never invents coordinates (ADR-0004).
6. **IDS study manifest**: not changed by this task. Recording the study's decisions
   (`truncate` for A:26-33 citing the paper, `model_loop` for 444-453) is for the
   maintainer; the test fixtures stay test inputs, and the acceptance tests use their
   own manifests with those decisions.
