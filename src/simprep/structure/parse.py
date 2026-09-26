"""Read mmCIF / PDB files with gemmi into the immutable :class:`Structure` model."""

from __future__ import annotations

from pathlib import Path

import gemmi

from simprep.structure.model import (
    Atom,
    Link,
    LinkPartner,
    ModifiedResidue,
    Residue,
    ResidueClass,
    ResidueId,
    SequenceReference,
    Structure,
    UnobservedResidue,
)

STRUCT_CONN = "_struct_conn."
MOD_RESIDUE = "_pdbx_struct_mod_residue."
UNOBSERVED = "_pdbx_unobs_or_zero_occ_residues."
CRYSTAL_GROW = "_exptl_crystal_grow."
POLY_SEQ = "_entity_poly_seq."
STRUCT_REF = "_struct_ref."
STRUCT_REF_SEQ = "_struct_ref_seq."
ANNOTATION_CATEGORIES = (STRUCT_CONN, MOD_RESIDUE, UNOBSERVED, POLY_SEQ)
UNOBSERVED_OCCUPANCY_FLAG = "1"  # mmCIF: 1 = unobserved, 0 = zero occupancy

ENTITY_CLASSES = {
    gemmi.EntityType.Polymer: ResidueClass.POLYMER,
    gemmi.EntityType.NonPolymer: ResidueClass.NONPOLYMER,
    gemmi.EntityType.Branched: ResidueClass.BRANCHED,
    gemmi.EntityType.Water: ResidueClass.WATER,
}


class StructureReadError(ValueError):
    """The input file cannot be represented faithfully by the structure model."""


def read_structure(path: Path) -> Structure:
    """Parse ``path`` (mmCIF or PDB, optionally gzipped) into a :class:`Structure`.

    Only single-model files are accepted: silently auditing model 1 of an ensemble
    would drop records without a decision.
    """
    gemmi_structure = gemmi.read_structure(str(path))
    if len(gemmi_structure) != 1:
        raise StructureReadError(
            f"{path}: {len(gemmi_structure)} models found; simprep v0.1 audits exactly one. "
            "Extract the model to audit into its own file first."
        )
    gemmi_structure.setup_entities()
    block = _annotation_block(path, gemmi_structure)
    return Structure(
        name=gemmi_structure.name,
        residues=_residues(gemmi_structure[0]),
        links=tuple(_link(row) for row in _rows(block, STRUCT_CONN)),
        modified_residues=tuple(_modified(row) for row in _rows(block, MOD_RESIDUE)),
        unobserved_residues=_unobserved(block),
        annotation_categories=frozenset(
            name.strip("_.") for name in ANNOTATION_CATEGORIES if _has_category(block, name)
        ),
        crystallization_details=_crystallization_details(block),
        polymer_sequences=_polymer_sequences(gemmi_structure),
        sequence_references=_sequence_references(block),
    )


def _sequence_references(block: gemmi.cif.Block) -> tuple[SequenceReference, ...]:
    """struct_ref_seq alignments with their database (struct_ref); rows lacking author or
    database ranges are skipped (they cannot map a position)."""
    databases = {row["id"]: row for row in _rows(block, STRUCT_REF)}
    references = []
    for row in _rows(block, STRUCT_REF_SEQ):
        database = databases.get(row.get("ref_id"), {})
        fields = (
            row.get("pdbx_auth_seq_align_beg"),
            row.get("pdbx_auth_seq_align_end"),
            row.get("db_align_beg"),
            row.get("pdbx_strand_id"),
        )
        if None in fields or not database.get("db_name"):
            continue
        references.append(
            SequenceReference(
                chain=row["pdbx_strand_id"],
                db_name=database["db_name"],
                db_accession=row.get("pdbx_db_accession")
                or database.get("pdbx_db_accession")
                or "",
                auth_begin=int(fields[0]),
                auth_end=int(fields[1]),
                db_begin=int(fields[2]),
            )
        )
    return tuple(references)


def _polymer_sequences(gemmi_structure: gemmi.Structure) -> tuple[tuple[str, tuple[str, ...]], ...]:
    """Full entity sequence per chain that has a polymer part."""
    sequences = []
    for chain in gemmi_structure[0]:
        span = chain.get_polymer()
        entity = gemmi_structure.get_entity_of(span) if span else None
        if entity is not None and entity.full_sequence:
            sequences.append((chain.name, tuple(entity.full_sequence)))
    return tuple(sequences)


def _crystallization_details(block: gemmi.cif.Block) -> str | None:
    """Crystallization conditions (all rows of exptl_crystal_grow), whitespace-normalized."""
    details = [row.get("pdbx_details") for row in _rows(block, CRYSTAL_GROW)]
    text = "; ".join(" ".join(d.split()) for d in details if d)
    return text or None


def _annotation_block(path: Path, gemmi_structure: gemmi.Structure) -> gemmi.cif.Block:
    """Return the mmCIF block holding annotations (converted from PDB if needed)."""
    if gemmi_structure.input_format == gemmi.CoorFormat.Pdb:
        return gemmi_structure.make_mmcif_document().sole_block()
    return gemmi.cif.read(str(path)).sole_block()


def _residues(model: gemmi.Model) -> tuple[Residue, ...]:
    return tuple(_residue(chain.name, gemmi_residue) for chain in model for gemmi_residue in chain)


def _residue(chain_name: str, gemmi_residue: gemmi.Residue) -> Residue:
    seqid = gemmi_residue.seqid
    return Residue(
        id=ResidueId(chain_name, seqid.num, seqid.icode.strip()),
        name=gemmi_residue.name,
        residue_class=_residue_class(gemmi_residue),
        atoms=tuple(_atom(atom) for atom in gemmi_residue),
        label_seq=gemmi_residue.label_seq,
    )


def _residue_class(gemmi_residue: gemmi.Residue) -> ResidueClass:
    if gemmi_residue.is_water():
        return ResidueClass.WATER
    try:
        return ENTITY_CLASSES[gemmi_residue.entity_type]
    except KeyError:
        raise StructureReadError(
            f"residue {gemmi_residue.name} {gemmi_residue.seqid} has no entity type; "
            "the file's entity records are missing or inconsistent."
        ) from None


def _atom(gemmi_atom: gemmi.Atom) -> Atom:
    pos = gemmi_atom.pos
    return Atom(
        name=gemmi_atom.name,
        element=gemmi_atom.element.name.upper(),
        position=(pos.x, pos.y, pos.z),
        occupancy=gemmi_atom.occ,
        b_iso=gemmi_atom.b_iso,
        altloc=gemmi_atom.altloc if gemmi_atom.has_altloc() else "",
    )


def _has_category(block: gemmi.cif.Block, category: str) -> bool:
    return category in block.get_mmcif_category_names()


def _rows(block: gemmi.cif.Block, category: str) -> list[dict[str, str | None]]:
    """Rows of an mmCIF category as dicts; CIF nulls ('?' and '.') become None."""
    table = block.find_mmcif_category(category)
    tags = [tag[len(category) :] for tag in table.tags]
    return [{tag: _value(row[index]) for index, tag in enumerate(tags)} for row in table]


def _value(raw: str) -> str | None:
    return None if raw in ("?", ".") else gemmi.cif.as_string(raw)


def _residue_id(row: dict, chain_key: str, seq_key: str, icode_key: str) -> ResidueId:
    return ResidueId(row[chain_key], int(row[seq_key]), row.get(icode_key) or "")


def _partner(row: dict, index: int) -> LinkPartner:
    prefix = f"ptnr{index}_"
    return LinkPartner(
        residue=_residue_id(
            row, f"{prefix}auth_asym_id", f"{prefix}auth_seq_id", f"pdbx_{prefix}PDB_ins_code"
        ),
        res_name=row.get(f"{prefix}auth_comp_id") or row[f"{prefix}label_comp_id"],
        atom_name=row[f"{prefix}label_atom_id"],
        altloc=row.get(f"pdbx_{prefix}label_alt_id") or "",
    )


def _link(row: dict) -> Link:
    distance = row.get("pdbx_dist_value")
    return Link(
        conn_id=row["id"],
        conn_type=row["conn_type_id"],
        partner1=_partner(row, 1),
        partner2=_partner(row, 2),
        distance_angstrom=float(distance) if distance else None,
    )


def _modified(row: dict) -> ModifiedResidue:
    return ModifiedResidue(
        residue=_residue_id(row, "auth_asym_id", "auth_seq_id", "PDB_ins_code"),
        res_name=row.get("auth_comp_id") or row["label_comp_id"],
        parent_res_name=row.get("parent_comp_id"),
        details=row.get("details"),
    )


def _unobserved(block: gemmi.cif.Block) -> tuple[UnobservedResidue, ...]:
    """Unobserved (not zero-occupancy) residues; zero-occupancy ones are still modelled."""
    return tuple(
        UnobservedResidue(
            residue=_residue_id(row, "auth_asym_id", "auth_seq_id", "PDB_ins_code"),
            res_name=row["auth_comp_id"],
            label_seq=int(row["label_seq_id"]) if row.get("label_seq_id") else None,
            is_polymer=row.get("polymer_flag") == "Y",
        )
        for row in _rows(block, UNOBSERVED)
        if row.get("occupancy_flag") == UNOBSERVED_OCCUPANCY_FLAG
    )
