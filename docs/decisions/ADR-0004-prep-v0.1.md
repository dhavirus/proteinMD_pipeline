# ADR-0004: Prep v0.1: how decisions become coordinates

- Status: accepted (TASK-004)
- Date: 2026-09-26

## Context

TASK-004 adds the first stage that changes coordinates: `simprep prep` applies a
fully decided manifest to its structure. The maintainer settled the scope (select,
remove and record only), the code location (`src/simprep/prep/`, a thin Colab wrapper in
`modules/`), `revert_to_parent` only with an explicit mapping, one system per altloc for
`keep_ensemble`, and mmCIF + PDB output. This ADR records the implementation choices
those decisions left open.

## Decisions

1. **The option table is data.** Every knowledge-base option carries `prep_action`
   (`apply` / `record` / `defer` / `unresolved`), and `record`/`defer` options carry
   `prep_stage` (rule schema enforces both). Prep reads them from the manifest's
   findings snapshot, so a manifest is applied with the options it was decided on. A
   snapshot without `prep_action` (made before knowledge base 0.3.0) is refused with a
   request to re-audit rather than guessed.
2. **Apply operations are keyed by option id** (`OPERATIONS` in `prep/plan.py`). A test
   checks that every `apply` option in `knowledge/` has one, so adding such an option
   without code fails CI instead of failing at prep time.
3. **Plan, then apply, then write.** `plan.py` (manifest -> actions, work order,
   exclusions, altloc choices, new links, reverts) and `apply.py` (plan -> systems) are
   pure and run on the immutable `Structure` model; gemmi is used only in `write.py`.
   The plan collects *every* problem (non-final decisions, missing altlocs, missing
   mappings, conflicts) and raises one `PrepError` listing them all.
4. **Conflicts are errors, not precedence rules.** A residue that one decision excludes
   and another links (`add_link`) or keeps for a later stage (`defer`) stops prep. No
   decision silently wins.
5. **Gate.** Prep runs only when no blocking finding is undecided and no decision uses a
   non-final option (`prep_action: unresolved`, e.g. `expert_review`), on any finding.
   Non-blocking undecided findings pass through unchanged and appear in the work order
   under stage `decision`.
6. **Accounting by comparison.** Counts are computed by comparing each output system
   with the input, per record type (residues per entity class, atom records,
   struct_conn, pdbx_struct_mod_residue, unobserved polymer residues):
   `in = passed_through + modified + excluded`, `out = passed_through + modified + added`.
   Atoms kept in a modified residue count as modified; atoms a modification dropped
   (other altlocs, mapped deletions) count as excluded. Tests re-read every written
   mmCIF and compare with the `out` column, and derive expected atom counts from the raw
   `atom_site` table independently of simprep.
7. **Writer.** gemmi writes the model as built from our `Structure`:
   - mmCIF is the primary output and round-trips the model exactly (entity classes,
     full polymer sequences, connections, modified-residue and unobserved-residue
     annotations). PDB has no branched entities or REMARK 465, so it carries residues,
     atoms and LINK/SSBOND records only.
   - No unit cell is written (the model has none; a placeholder CRYST1 would be read as
     a periodic box by MD tools).
   - `pdbx_dist_value` is recomputed by gemmi from the coordinates, not copied.
   - No timestamps, so identical inputs give byte-identical files (tested).
   - Connection types outside covale/disulf/hydrog/metalc raise `WriteError` instead of
     being written as something else.
8. **Colab module.** `modules/prep_module.py` adds the run bookkeeping only (Drive
   layout, `config.json` with config + manifest + input SHA-256 + simprep commit,
   `run_id`, local work dir then bulk copy). Resume is per run: a run whose recorded
   outputs all exist with matching SHA-256 is skipped; otherwise prep (seconds on the
   panel) runs again and replaces the outputs. Resuming with a changed config, manifest
   or input is refused.

## Consequences

- A later stage (modelling, chemistry, parameterization) consumes `prep_record.json`'s
  work order; adding a modelling tool needs its own ADR (TASK-004 decision 2).
- New `apply` options need an operation in `OPERATIONS` (enforced by the test above).
- The record's atom `modified` count is per residue, not per atom; if a later task needs
  per-atom change tracking, it changes that definition here.
