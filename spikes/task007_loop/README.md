# TASK-007 spike: PDBFixer loop building (not part of the package)

`pdbfixer_loop.py` builds the unobserved internal loop(s) of a prepared structure with
PDBFixer (terminal gaps are skipped) and writes a PDB. It is the script behind the
"Spike" section of `docs/tasks/TASK-007-modelling.md`.

```bash
python -m venv .spike && .spike/bin/pip install openmm==8.6.1 pdbfixer==1.12.0
# input: a prepared wild type written by `simprep prep` (system.pdb, which carries SEQRES)
.spike/bin/python spikes/task007_loop/pdbfixer_loop.py OUT/wt/system.pdb wt_loop.pdb
```

Findings (TASK-006/007 session, 5FQL prep fixture, gap A:444-453):
- 10 residues added in about 18 s; pre-existing atoms unmoved; C-N bonds 1.29-1.43 A.
- Junction geometry distorted: CA443-CA444 2.61 A, CA444-CA445 4.46 A; odd phi/psi.
- Two runs differ by up to 6.7 A even with Python `random` seeded: the first
  implementation step of TASK-007 is to find and seed the remaining randomness
  (numpy? the OpenMM platform used by `addMissingAtoms`?). If byte-identical output
  cannot be reached, TASK-007 stops and reports (decision 2).

Resolved in the TASK-007 implementation session (ADR-0007):
- The remaining randomness was PDBFixer's own Langevin integrator (unseeded unless
  `addMissingAtoms(seed=...)`) on the default multithreaded platform. With the seed and
  `PDBFixer(..., platform=Reference)` two runs are byte-identical.
- The placement also contains D residues (Phe444, Leu447, Tyr452, Leu453), so simprep
  follows it with a restrained, sterics-only minimization (trans omega and L chirality
  restraints); see `src/simprep/model/` and ADR-0007.
