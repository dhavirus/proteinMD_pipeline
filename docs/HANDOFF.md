# Handoff: state of simprep after TASK-006 (2026-09-26)

Written at the end of the cloud session so work can continue from a local checkout
(Claude Code on a PC). `CLAUDE.md` still governs everything; this file is the map.

## Where things are

| | |
|---|---|
| `main` | `0fc1f0e`: TASK-001 to TASK-006 merged (PRs #1-#6), CI green |
| branch `claude/quirky-mayer-t7x0sl` | `main` + the accepted TASK-007 spec + this file + `spikes/task007_loop/` (not merged yet) |
| next task | **TASK-007 modelling**, spec accepted: `docs/tasks/TASK-007-modelling.md` |

Done so far:

| task | what it gave | ADR |
|---|---|---|
| 001 | audit: deterministic detectors, findings.json, panel regression (5FQL, 3KS3, 1HZH, 6OIM, 1FO8) | 0001, 0002 |
| 002 | review page (static, Mol*, ajv): decisions, regions, manifest export | 0003 |
| 003 | covalent-contact candidates the annotations miss (1HZH N-glycans) | - |
| 004 | `simprep prep`: apply decisions (select, remove, record), work order, record accounting | 0004 |
| 005 | `simprep variants`: rotamer candidates as findings, WT vs R468Q/R468W builds | 0005 |
| 006 | relaxation (OpenMM, Reference platform, byte-identical) of variants + matched WT | 0006 |

## Local setup

```bash
git clone https://github.com/dhavirus/proteinMD_pipeline.git
cd proteinMD_pipeline
git checkout claude/quirky-mayer-t7x0sl       # or main, once this branch is merged
python3.11 -m venv .venv && . .venv/bin/activate
python -m pip install -e '.[dev,relax]'       # relax = openmm==8.6.1 (about 50 MB)
pytest                                        # about 340 tests, 2-3 min (OpenMM relaxation tests included)
(cd frontend && node --test "tests/*.test.mjs")   # Node 22
ruff check . && ruff format --check .
```

Optional: the browser smoke test (`node frontend/tests/smoke/smoke.mjs`) needs Playwright
and a Chromium; with network it loads Mol* 5.11.0 and ajv-dist 8.17.1 from the pinned CDN
(in the cloud it ran against a local copy via `SIMPREP_CDN_DIR`). It is not part of CI.

The panel structures are committed (`tests/panel/*.cif.gz`); no test needs the network.

## How the work has been run (keep doing this)

1. Draft a spec in `docs/tasks/TASK-00N-*.md` with a short **spike** first (measure on the
   real data before deciding), open questions with defaults, `[VERIFY]` on anything not
   read in the session.
2. Maintainer accepts or changes the defaults; the spec is finalized with a "Resolved
   decisions" section.
3. Implement on a branch, small commits, ADR for any architectural or tool choice,
   expected test values written by hand (never from the implementation), `reviewed: false`.
4. PR with: what was done, what was not and why, open questions, `[VERIFY]` items. CI green
   before the PR.

## Next: TASK-007 (modelling), first steps

- Step 1 (decision 2): make PDBFixer loop building deterministic. Start from
  `spikes/task007_loop/` (script + findings). If byte-identical output cannot be
  reached, stop and report.
- Then: `src/simprep/model/` stage between prep and variants; `model_loop` becomes
  `apply`; modelled atoms written with occupancy 0.00 and listed in the record; one loop
  shared by WT, R468Q, R468W; `truncate` = charged termini recorded for topology.
- Acceptance: 5FQL with `truncate` for A:26-33 (the paper: signal peptide + propeptide
  are cleaved, Thr34 is the mature N-terminus) and `model_loop` for A:444-453.
- `pdbfixer==1.12.0` will need to join the `relax` extra (or a new `model` extra), via the ADR.

## Rough roadmap after TASK-007 (estimate, revise as tasks land)

| # | stage | notes |
|---|---|---|
| 008 | chemistry: FGly (ALS A:84) -> gem-diol (FGH, CCD `DDZ` per the paper, `[VERIFY]` in CCD) | work-order item `nonstandard_residues/A:84` |
| 009 | protonation / pKa | the paper: lysosomal pH ~4.8; modelled and moved residues have no H |
| 010 | parameterization | Ca2+ (12-6-4 per the manifest decision), Cl-, FGH |
| 011 | system assembly | solvation, ions, full topology: MD-ready WT and variants |
| 012 | equilibration / MD on Colab | GPU, checkpoint/resume, Colab conventions in CLAUDE.md |
| 013+ | comparative free-energy study | CVs, sampling method (WSME-L / enhanced MD / FEP): study decisions |

About six more tasks after TASK-007, give or take one or two.

## Open items for the maintainer (not blocking the code)

- **Merge** branch `claude/quirky-mayer-t7x0sl` (spec + handoff + spike) into `main`.
- **Review fixtures** marked `reviewed: false`: `tests/panel/expected_findings/*.yaml`,
  `tests/panel/prep_fixtures/*.yaml`, `tests/panel/variant_fixtures/5FQL.yaml`.
- **Study manifest**: record the IDS study's own decisions (for example `truncate` for
  A:26-33 citing Demydchuk et al. 2017, doi:10.1038/ncomms15786; `model_loop` for 444-453).
- **Colab**: notebooks `01_prep` and `02_variants` have never been run on Colab; set
  `repo_commit` to a real commit first.
- **GitHub Pages** for the review page (workflow exists; Pages must be enabled in the
  repository settings), and one smoke run against the live CDN.
- **`[VERIFY]` items** in the knowledge base (not read in a session): `metals.yaml`
  (Harding 2006 shell cutoffs, open-site heuristics, 12-6-4 for Na+/K+),
  `nonstandard_residues.yaml` (DDZ / ALS identity in the CCD, gem-diol as default),
  `missing_residues.yaml` (3.8 A CA-CA), `side_chains.yaml` (Engh & Huber geometry,
  Lovell rotamer table), `variant_build.yaml` (Bondi radii, 0.4 A on heavy atoms),
  `relaxation.yaml` (TIP3P reference, amber14-all.xml = ff14SB), and the 1HZH Asn297
  attachment note in `tests/panel/expected_findings/1HZH.yaml`.
  Verified in PubMed during the sessions: Lovell 2000, Word 1999, Parsons 2005,
  Shapovalov & Dunbrack 2011, Maier 2015 (ff14SB), Eastman 2024 (OpenMM 8),
  Demydchuk 2017 (5FQL, full text).

## Things that did not survive the cloud container

Scratch files in `/tmp` (spike outputs, the CDN copy used by the smoke test, the
PDBFixer venv) are gone; everything needed to reproduce them is in the repository
(`spikes/`, the task specs' "Spike" sections, and the commands above).
