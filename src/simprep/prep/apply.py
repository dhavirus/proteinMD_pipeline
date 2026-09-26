"""Apply a prep plan to a structure (pure: immutable model in, immutable models out)."""

from __future__ import annotations

from dataclasses import dataclass, replace

from simprep.prep.plan import PrepPlan
from simprep.structure.model import Link, Residue, ResidueId, Structure


@dataclass(frozen=True)
class System:
    """One prepared system: ``altloc`` is the ensemble altloc it keeps (None without one).
    What changed is derived by comparing ``structure`` with the input (see record.py)."""

    name: str
    altloc: str | None
    structure: Structure


def apply_plan(structure: Structure, plan: PrepPlan) -> list[System]:
    """One system, or one per ensemble altloc (TASK-004 decision 4)."""
    altlocs = plan.ensemble_altlocs or (None,)
    return [
        System(_system_name(altloc), altloc, _prepared(structure, plan, altloc))
        for altloc in altlocs
    ]


def _system_name(ensemble_altloc: str | None) -> str:
    return "system" if ensemble_altloc is None else f"system_altloc_{ensemble_altloc}"


def _prepared(structure: Structure, plan: PrepPlan, ensemble_altloc: str | None) -> Structure:
    kept = tuple(
        _prepared_residue(residue, plan, ensemble_altloc)
        for residue in structure.residues
        if residue.id not in plan.excluded
    )
    kept_ids = {residue.id for residue in kept}
    links = tuple(
        _collapse_link(link, plan, ensemble_altloc)
        for link in structure.links + plan.new_links
        if _link_survives(link, kept_ids, plan, ensemble_altloc)
    )
    modified = tuple(
        m for m in structure.modified_residues if m.residue in kept_ids - set(plan.reverts)
    )
    return replace(structure, residues=kept, links=links, modified_residues=modified)


def kept_altloc(residue_id: ResidueId, plan: PrepPlan, ensemble_altloc: str | None) -> str | None:
    """The altloc this residue keeps in this system, or None if its altlocs stay."""
    if residue_id in plan.ensemble_residues:
        return ensemble_altloc
    return plan.altloc_choice.get(residue_id)


def _prepared_residue(residue: Residue, plan: PrepPlan, ensemble_altloc: str | None) -> Residue:
    altloc = kept_altloc(residue.id, plan, ensemble_altloc)
    atoms = residue.atoms
    if altloc is not None:
        atoms = tuple(replace(a, altloc="") for a in atoms if a.altloc in ("", altloc))
    if residue.id in plan.reverts:
        return _reverted(replace(residue, atoms=atoms), plan.reverts[residue.id])
    return replace(residue, atoms=atoms)


def _reverted(residue: Residue, mapping: dict) -> Residue:
    """Rename / delete atoms per an explicit mapping (knowledge/residue_mappings.yaml)."""
    elements = mapping.get("rename_elements", {})
    atoms = []
    for atom in residue.atoms:
        if atom.name in mapping["delete"]:
            continue
        name = mapping["rename"].get(atom.name, atom.name)
        atoms.append(replace(atom, name=name, element=elements.get(name, atom.element)))
    return replace(residue, name=mapping["to"], atoms=tuple(atoms))


def _link_survives(link: Link, kept: set[ResidueId], plan: PrepPlan, ensemble_altloc) -> bool:
    for partner in (link.partner1, link.partner2):
        if partner.residue not in kept:
            return False
        altloc = kept_altloc(partner.residue, plan, ensemble_altloc)
        if altloc is not None and partner.altloc not in ("", altloc):
            return False
    return True


def _collapse_link(link: Link, plan: PrepPlan, ensemble_altloc: str | None) -> Link:
    """Clear partner altlocs on residues whose altlocs were collapsed."""

    def partner(p):
        collapsed = kept_altloc(p.residue, plan, ensemble_altloc) is not None
        return replace(p, altloc="") if collapsed else p

    return replace(link, partner1=partner(link.partner1), partner2=partner(link.partner2))
