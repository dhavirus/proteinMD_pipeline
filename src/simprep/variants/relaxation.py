"""Relaxation of variant sites and the matched wild type (TASK-006).

Variants are grouped by the residues they mutate (their site). Per site, one shell is
used for everyone: the union of the wild type's shell and each variant's shell, so the
relaxed wild type and every relaxed variant move exactly the same residues.
"""

from __future__ import annotations

from dataclasses import dataclass

from simprep.findings import Finding, atom_ref, ev_count, ev_distance, ev_label, residue_locus
from simprep.relax.openmm_run import RelaxOutcome, relax
from simprep.relax.shell import Shell, build_shell, shell_residues
from simprep.rules import RuleSet
from simprep.severity import apply_context
from simprep.structure.model import ResidueId, Structure
from simprep.variants.build import BACKBONE_KEPT
from simprep.variants.check import ClashCriterion, Contact, clashes, contacts, neighbours

FAMILY = "relaxation"
CRITERION_RULE = "variant_build.rotamer_choice"
ENERGY_NOTE = (
    "Total potential energy of the calculation, dominated by fixed-fixed terms; not "
    "comparable between systems. Only the change and the final force on mobile atoms "
    "are meaningful."
)


@dataclass(frozen=True)
class Relaxed:
    """One relaxed system with what the relaxation did and what clashes remain."""

    name: str
    outcome: RelaxOutcome
    shell: Shell
    clashes_before: tuple[Contact, ...]
    clashes_after: tuple[Contact, ...]


def site_label(sites: frozenset[ResidueId]) -> str:
    return "_".join(f"{r.chain}{r.seq_num}{r.ins_code}" for r in sorted(sites))


def union_shell(structures: list[Structure], sites: frozenset, protocol: dict) -> frozenset:
    """Residues in the shell of any of ``structures`` (wild type and its variants)."""
    radius = protocol["mobile_shell_angstrom"]
    return frozenset().union(*(shell_residues(s, sites, radius) for s in structures))


def relax_system(name: str, structure: Structure, where: tuple, context: tuple) -> Relaxed:
    """Relax ``structure`` at ``sites`` over ``residues``; ``context`` is (protocol, ruleset)."""
    sites, residues = where
    protocol, ruleset = context
    shell = build_shell(structure, sites, residues, protocol["nonbonded_cutoff_nm"])
    outcome = relax(structure, shell, protocol)
    criterion = _criterion(ruleset)
    return Relaxed(
        name,
        outcome,
        shell,
        site_clashes(structure, sites, criterion),
        site_clashes(outcome.structure, sites, criterion),
    )


def _criterion(ruleset: RuleSet) -> ClashCriterion:
    rule = next(r for r in ruleset.rules if r.rule_id == CRITERION_RULE)
    return ClashCriterion.from_matcher(rule.matcher)


def site_clashes(structure: Structure, sites: frozenset, criterion: ClashCriterion) -> tuple:
    """Clashes of the site residues' side chains (beyond CB) with everything else."""
    found = []
    for rid in sorted(sites):
        residue = structure.residue(rid)
        side = tuple(a for a in residue.atoms if a.is_heavy and a.name not in BACKBONE_KEPT)
        nearby = [
            (r, a) for r, a in neighbours(structure, rid) if a.element in criterion.radii_angstrom
        ]
        found += clashes(contacts(side, nearby, criterion), criterion)
    return tuple(found)


def residual_findings(relaxed: Relaxed, structure: Structure, context: tuple) -> list[dict]:
    """A warn finding when clashes remain after relaxation (clashes only flag)."""
    ruleset, config = context
    if not relaxed.clashes_after:
        return []
    rule = ruleset.family(FAMILY)[0]
    sites = sorted(relaxed.shell.sites)
    finding = Finding(
        id=f"{FAMILY}/{relaxed.name}/{site_label(frozenset(sites))}",
        rule=rule,
        title=f"{relaxed.name}: {len(relaxed.clashes_after)} clash(es) remain after relaxation",
        locus=residue_locus(structure.residue(sites[0])),
        anchor_residues=tuple(sites),
        evidence=(
            ev_label("system", relaxed.name),
            ev_count("clashes_before", len(relaxed.clashes_before)),
            ev_count("clashes_after", len(relaxed.clashes_after)),
            *(_clash_item(c, relaxed.outcome.structure, sites) for c in relaxed.clashes_after),
        ),
    )
    return [f.to_dict() for f in apply_context([finding], structure, config)]


def _clash_item(contact: Contact, structure: Structure, sites: list) -> dict:
    site = next(
        structure.residue(r) for r in sites if contact.new_atom in structure.residue(r).atoms
    )
    return ev_distance(
        "residual_clash",
        contact.distance_angstrom,
        [atom_ref(site, contact.new_atom), atom_ref(contact.residue, contact.atom)],
        f"{contact.new_atom.name} - {contact.residue.name} {contact.residue.id.label()} "
        f"{contact.atom.name}, overlap {contact.overlap_angstrom:.2f} A",
    )


def relaxation_entry(relaxed: Relaxed) -> dict:
    """The record's account of one relaxation."""
    outcome, shell = relaxed.outcome, relaxed.shell
    initial, final = outcome.energy_kj_per_mol
    return {
        "sites": [rid.to_dict() for rid in sorted(shell.sites)],
        "shell_residues": len(shell.residues),
        "mobile_atoms": len(shell.mobile),
        "frozen": [
            {**rid.to_dict(), "near": near, "distance_angstrom": distance}
            for rid, near, distance in shell.frozen
        ],
        "left_out": list(outcome.left_out),
        "moved_residues": [rid.to_dict() for rid in outcome.moved_residues],
        "moved_atoms": outcome.moved_atoms,
        "max_displacement_angstrom": round(outcome.max_displacement_angstrom, 3),
        "rms_displacement_angstrom": round(outcome.rms_displacement_angstrom, 3),
        "clashes_before": len(relaxed.clashes_before),
        "clashes_after": len(relaxed.clashes_after),
        "energy_kj_per_mol": {
            "initial": round(initial, 1),
            "final": round(final, 1),
            "note": ENERGY_NOTE,
        },
        "max_mobile_force_kj_per_mol_nm": round(outcome.max_mobile_force_kj_per_mol_nm, 2),
        "openmm_version": outcome.openmm_version,
    }


def protonation_items(relaxed: Relaxed, finding_id: str, skip: frozenset) -> list[dict]:
    """Work-order items for residues whose hydrogens were dropped because they moved."""
    return [
        {
            "stage": "protonation",
            "finding_id": finding_id,
            "option_id": None,
            "description": f"{rid.label()}: moved in relaxation, hydrogens removed",
            "residues": [rid.to_dict()],
        }
        for rid in relaxed.outcome.moved_residues
        if rid not in skip
    ]
