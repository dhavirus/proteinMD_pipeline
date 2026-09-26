"""Initial placement of loop residues with PDBFixer (I/O edge; TASK-007, ADR-0007).

PDBFixer reads the prepared system as PDB (with SEQRES) and builds only the gaps decided
``model_loop``; nothing else is added. Its Langevin integrator is seeded and it runs on
the Reference platform, so the placement is byte-identical between runs. Only the new
residues' heavy atoms are taken back; every existing atom keeps its own coordinates.
"""

from __future__ import annotations

import tempfile
from importlib.metadata import version as installed_version
from pathlib import Path

from simprep.model.loops import Gap, ModelError
from simprep.prep.write import write_system
from simprep.structure.model import Atom, ResidueId, Structure

COORDINATE_DECIMALS = 3  # the precision of the written files
INPUT_NAME = "input"


def _pdbfixer(version: str):
    try:
        import openmm
        import pdbfixer
        from openmm import app, unit
    except ImportError as error:
        raise ModelError(
            "loop modelling needs PDBFixer and OpenMM: pip install 'simprep[relax]' "
            "(pdbfixer==1.12.0, openmm==8.6.1)"
        ) from error
    if installed_version("pdbfixer") != version:
        raise ModelError(
            f"the modelling protocol names PDBFixer {version}, installed is "
            f"{installed_version('pdbfixer')}; install that version or change the protocol"
        )
    return openmm, app, unit, pdbfixer.PDBFixer


def place_loops(structure: Structure, gaps: tuple[Gap, ...], protocol: dict) -> dict:
    """Heavy atoms of every gap residue, by residue id, placed by PDBFixer under the
    protocol's ``placement`` settings."""
    openmm, app, unit, fixer_class = _pdbfixer(protocol["tool_version"])
    platform = openmm.Platform.getPlatformByName(protocol["platform"])
    with tempfile.TemporaryDirectory() as directory:
        write_system(structure, Path(directory), INPUT_NAME)
        fixer = fixer_class(filename=str(Path(directory) / f"{INPUT_NAME}.pdb"), platform=platform)
    fixer.findMissingResidues()
    fixer.missingResidues = _requested(fixer, gaps)
    fixer.findMissingAtoms()
    fixer.missingAtoms, fixer.missingTerminals = {}, {}
    fixer.addMissingAtoms(seed=protocol["seed"])
    positions = fixer.positions.value_in_unit(unit.angstrom)
    return {rid: atoms for gap in gaps for rid, atoms in _placed(fixer, gap, positions).items()}


def _key(residue) -> tuple[str, int, str]:
    return residue.chain.id, int(residue.id), (residue.insertionCode or "").strip()


def _flank_key(rid: ResidueId) -> tuple[str, int, str]:
    return rid.chain, rid.seq_num, rid.ins_code


def _requested(fixer, gaps: tuple[Gap, ...]) -> dict:
    """PDBFixer's missing residues restricted to ``gaps``; ModelError if one is absent or
    PDBFixer reads a different sequence."""
    chains = list(fixer.topology.chains())
    found = {}
    for (chain_index, index), names in fixer.missingResidues.items():
        residues = list(chains[chain_index].residues())
        if 0 < index < len(residues):
            found[(_key(residues[index - 1]), _key(residues[index]))] = (
                (chain_index, index),
                names,
            )
    requested = {}
    for gap in gaps:
        key = tuple(_flank_key(rid) for rid in gap.flanks)
        expected = [member.res_name for member in gap.run.members]
        if key not in found or found[key][1] != expected:
            raise ModelError(
                f"{gap.finding_id}: PDBFixer did not find the gap {gap.label} with sequence "
                f"{'-'.join(expected)} between the flanks (found {found.get(key, (None, None))[1]})"
            )
        requested[found[key][0]] = found[key][1]
    return requested


def _placed(fixer, gap: Gap, positions) -> dict:
    """The new residues between the gap's flanks, as heavy atoms in A."""
    before, after = (_flank_key(rid) for rid in gap.flanks)
    residues = [r for chain in fixer.topology.chains() for r in chain.residues()]
    keys = [_key(r) for r in residues]
    start, end = keys.index(before), keys.index(after)
    new = residues[start + 1 : end]
    if [r.name for r in new] != [m.res_name for m in gap.run.members]:
        raise ModelError(f"{gap.finding_id}: PDBFixer's output does not hold {gap.label}")
    return {
        member.residue: tuple(
            _atom(a, positions) for a in residue.atoms() if a.element.symbol != "H"
        )
        for member, residue in zip(gap.run.members, new, strict=True)
    }


def _atom(atom, positions) -> Atom:
    position = tuple(round(c, COORDINATE_DECIMALS) for c in positions[atom.index])
    return Atom(atom.name, atom.element.symbol.upper(), position, 1.0, 0.0)
