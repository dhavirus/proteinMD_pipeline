"""TASK-007 spike: build the internal unobserved loop(s) of a prepared structure with
PDBFixer and write a PDB. Terminal gaps are skipped. Not part of the package.

Usage: python pdbfixer_loop.py INPUT.pdb OUTPUT.pdb
"""

import random
import sys
import time

from openmm.app import PDBFile
from pdbfixer import PDBFixer

SEED = 1  # Python's random only; the spike showed this is not enough for determinism


def main(source: str, target: str) -> None:
    random.seed(SEED)
    start = time.time()
    fixer = PDBFixer(filename=source)
    fixer.findMissingResidues()
    print("missing (chain index, position):", {k: len(v) for k, v in fixer.missingResidues.items()})
    chain_length = len(list(next(fixer.topology.chains()).residues()))
    fixer.missingResidues = {
        key: names for key, names in fixer.missingResidues.items() if 0 < key[1] < chain_length
    }
    print("kept internal only:", fixer.missingResidues)
    fixer.findMissingAtoms()
    fixer.missingAtoms = {}
    fixer.missingTerminals = {}
    fixer.addMissingAtoms()
    with open(target, "w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    print(f"built in {time.time() - start:.1f}s")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
