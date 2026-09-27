"""pKa estimates with PROPKA (I/O edge; TASK-009, ADR-0009).

PROPKA reads the system as PDB, heavy atoms only and without waters (all deposited
hydrogens are stripped first, TASK-009 decision 3). It ignores metal ions and residues it
has no parameters for; states.py flags residues near them.
"""

from __future__ import annotations

import tempfile
from dataclasses import replace
from importlib.metadata import version as installed_version
from pathlib import Path

from simprep.prep.write import write_system
from simprep.protonate.states import Estimate, ProtonationError
from simprep.structure.model import ResidueClass, ResidueId, Structure

INPUT_NAME = "input"
AVERAGE_CONFORMATION = "AVR"


def without_hydrogens(structure: Structure) -> Structure:
    residues = tuple(replace(r, atoms=r.heavy_atoms) for r in structure.residues)
    return replace(structure, residues=residues)


def estimate(structure: Structure, method: dict) -> list[Estimate]:
    """PROPKA's pKa per titratable group of ``structure`` (termini excluded)."""
    try:
        import propka.run
    except ImportError as error:
        raise ProtonationError(
            "protonation needs PROPKA: pip install 'simprep[relax]' (propka==3.5.1)"
        ) from error
    if installed_version("propka") != method["pka_tool_version"]:
        raise ProtonationError(
            f"the protonation method names PROPKA {method['pka_tool_version']}, installed is "
            f"{installed_version('propka')}; install that version or change the method"
        )
    dry = [r for r in structure.residues if r.residue_class is not ResidueClass.WATER]
    heavy = without_hydrogens(replace(structure, residues=tuple(dry)))
    with tempfile.TemporaryDirectory() as directory:
        write_system(heavy, Path(directory), INPUT_NAME)
        path = Path(directory) / f"{INPUT_NAME}.pdb"
        molecule = propka.run.single(str(path), optargs=["--quiet"], write_pka=False)
    return [
        Estimate(
            ResidueId(g.atom.chain_id, g.atom.res_num, g.atom.icode.strip()),
            g.residue_type,
            g.pka_value,
            g.model_pka,
        )
        for g in molecule.conformations[AVERAGE_CONFORMATION].groups
        if g.residue_type not in ("N+", "C-")
    ]
