"""Write a prepared :class:`Structure` as mmCIF (primary) and PDB via gemmi (I/O edge).

The output carries what prep keeps: atoms, entity classes, full polymer sequences,
connection records (annotated ones that survived plus ``add_link`` bonds) and
modified-residue annotations. It carries no timestamps, so the same structure always
gives byte-identical files. It carries no unit cell either: the model has none, and a
placeholder CRYST1 would be read as a periodic box by MD tools. ``pdbx_dist_value`` is
recomputed by gemmi from the coordinates rather than copied from the input.
"""

from __future__ import annotations

from pathlib import Path

import gemmi

from simprep.provenance import sha256_file
from simprep.structure.model import Link, LinkPartner, Residue, ResidueClass, Structure

ENTITY_TYPES = {
    ResidueClass.POLYMER: gemmi.EntityType.Polymer,
    ResidueClass.NONPOLYMER: gemmi.EntityType.NonPolymer,
    ResidueClass.BRANCHED: gemmi.EntityType.Branched,
    ResidueClass.WATER: gemmi.EntityType.Water,
}
CONNECTION_TYPES = {
    "covale": gemmi.ConnectionType.Covale,
    "disulf": gemmi.ConnectionType.Disulf,
    "hydrog": gemmi.ConnectionType.Hydrog,
    "metalc": gemmi.ConnectionType.MetalC,
}
FORMAT_SUFFIXES = (("mmcif", ".cif"), ("pdb", ".pdb"))
FORMAT_BY_SUFFIX = {suffix: fmt for fmt, suffix in FORMAT_SUFFIXES}
UNOBSERVED_CATEGORY = "_pdbx_unobs_or_zero_occ_residues."


class WriteError(ValueError):
    """The prepared structure holds something the output formats cannot carry faithfully."""


def write_system(structure: Structure, out_dir: Path, name: str) -> list[dict]:
    """Write ``<name>.cif`` and ``<name>.pdb`` into ``out_dir``; return path, format and
    SHA-256 of each file (paths relative to ``out_dir``). Unobserved residues are written
    to the mmCIF only; the PDB file has no branched entities and no REMARK 465."""
    gemmi_structure = to_gemmi(structure)
    files = []
    for fmt, suffix in FORMAT_SUFFIXES:
        path = out_dir / f"{name}{suffix}"
        _write(gemmi_structure, path, structure.unobserved_residues)
        files.append({"path": path.name, "format": fmt, "sha256": sha256_file(path)})
    return files


def _write(gemmi_structure: gemmi.Structure, path: Path, unobserved: tuple) -> None:
    fmt = FORMAT_BY_SUFFIX[path.suffix]
    if fmt == "mmcif":
        groups = gemmi.MmcifOutputGroups(True)
        groups.cell = groups.symmetry = False
        document = gemmi_structure.make_mmcif_document(groups)
        _add_unobserved(document.sole_block(), unobserved)
        document.write_file(str(path))
    else:
        gemmi_structure.write_pdb(str(path), gemmi.PdbWriteOptions(cryst1_record=False))


def to_gemmi(structure: Structure) -> gemmi.Structure:
    """A single-model gemmi structure with entities, sequences and connections set up."""
    gemmi_structure = gemmi.Structure()
    gemmi_structure.name = structure.name
    model = gemmi.Model("1")
    for chain_name in dict.fromkeys(r.id.chain for r in structure.residues):
        chain = gemmi.Chain(chain_name)
        for residue in (r for r in structure.residues if r.id.chain == chain_name):
            chain.add_residue(_gemmi_residue(residue))
        model.add_chain(chain)
    gemmi_structure.add_model(model)
    gemmi_structure.setup_entities()
    _set_sequences(gemmi_structure, dict(structure.polymer_sequences))
    gemmi_structure.connections = gemmi.ConnectionList([_connection(k) for k in structure.links])
    gemmi_structure.mod_residues = [_mod_residue(m) for m in structure.modified_residues]
    return gemmi_structure


def _gemmi_residue(residue: Residue) -> gemmi.Residue:
    gemmi_residue = gemmi.Residue()
    gemmi_residue.name = residue.name
    gemmi_residue.seqid = gemmi.SeqId(residue.id.seq_num, residue.id.ins_code or " ")
    gemmi_residue.entity_type = ENTITY_TYPES[residue.residue_class]
    gemmi_residue.het_flag = "A" if residue.residue_class is ResidueClass.POLYMER else "H"
    if residue.label_seq is not None:
        gemmi_residue.label_seq = residue.label_seq
    for atom in residue.atoms:
        gemmi_residue.add_atom(_gemmi_atom(atom))
    return gemmi_residue


def _gemmi_atom(atom) -> gemmi.Atom:
    gemmi_atom = gemmi.Atom()
    gemmi_atom.name = atom.name
    gemmi_atom.element = gemmi.Element(atom.element)
    gemmi_atom.pos = gemmi.Position(*atom.position)
    gemmi_atom.occ = atom.occupancy
    gemmi_atom.b_iso = atom.b_iso
    gemmi_atom.altloc = atom.altloc or "\0"
    return gemmi_atom


def _set_sequences(gemmi_structure: gemmi.Structure, sequences: dict) -> None:
    """Full entity sequences (observed + unobserved) from the input, per polymer entity."""
    subchain_chain = {
        residue.subchain: chain.name for chain in gemmi_structure[0] for residue in chain
    }
    for entity in gemmi_structure.entities:
        if entity.entity_type != gemmi.EntityType.Polymer:
            continue
        chains = [subchain_chain[s] for s in entity.subchains if subchain_chain.get(s) in sequences]
        if chains:
            entity.full_sequence = list(sequences[chains[0]])


def _connection(link: Link) -> gemmi.Connection:
    if link.conn_type not in CONNECTION_TYPES:
        raise WriteError(
            f"struct_conn {link.conn_id} has type {link.conn_type!r}, which the writer does "
            f"not map (known: {', '.join(CONNECTION_TYPES)}); add it to CONNECTION_TYPES"
        )
    connection = gemmi.Connection()
    connection.name = link.conn_id
    connection.type = CONNECTION_TYPES[link.conn_type]
    connection.partner1 = _address(link.partner1)
    connection.partner2 = _address(link.partner2)
    if link.distance_angstrom is not None:
        connection.reported_distance = link.distance_angstrom
    return connection


def _address(partner: LinkPartner) -> gemmi.AtomAddress:
    residue = partner.residue
    return gemmi.AtomAddress(
        residue.chain,
        gemmi.SeqId(residue.seq_num, residue.ins_code or " "),
        partner.res_name,
        partner.atom_name,
        partner.altloc or "\0",
    )


def _mod_residue(modified) -> gemmi.ModRes:
    mod = gemmi.ModRes()
    mod.chain_name = modified.residue.chain
    mod.res_id = gemmi.ResidueId()
    mod.res_id.name = modified.res_name
    mod.res_id.seqid = gemmi.SeqId(modified.residue.seq_num, modified.residue.ins_code or " ")
    mod.parent_comp_id = modified.parent_res_name or ""
    mod.details = modified.details or ""
    return mod


def _add_unobserved(block: gemmi.cif.Block, unobserved: tuple) -> None:
    """pdbx_unobs_or_zero_occ_residues rows for the unobserved residues (flag 1)."""
    if not unobserved:
        return
    rows = [
        {
            "id": str(index),
            "PDB_model_num": "1",
            "polymer_flag": "Y" if u.is_polymer else "N",
            "occupancy_flag": "1",
            "auth_asym_id": u.residue.chain,
            "auth_comp_id": u.res_name,
            "auth_seq_id": str(u.residue.seq_num),
            "PDB_ins_code": u.residue.ins_code or None,
            "label_seq_id": None if u.label_seq is None else str(u.label_seq),
        }
        for index, u in enumerate(unobserved, start=1)
    ]
    block.set_mmcif_category(UNOBSERVED_CATEGORY, {key: [r[key] for r in rows] for key in rows[0]})
