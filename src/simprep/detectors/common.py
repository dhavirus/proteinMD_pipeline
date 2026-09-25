"""Helpers shared by several detectors."""

from __future__ import annotations

from simprep.rules import Rule, RuleSet
from simprep.structure.model import Atom, Residue, ResidueClass

GENERIC_NONSTANDARD_KIND = "nonstandard_polymer_residue"


def standard_polymer_residues(ruleset: RuleSet) -> frozenset[str]:
    """The standard residue names, defined once by the generic non-standard-residue rule."""
    for rule in ruleset.family("nonstandard_residues"):
        if rule.matcher["kind"] == GENERIC_NONSTANDARD_KIND:
            return frozenset(rule.matcher["standard_residues"])
    raise LookupError(
        "knowledge base has no nonstandard_polymer_residue rule defining the standard residues"
    )


def is_standard_polymer(residue: Residue, standard: frozenset[str]) -> bool:
    return residue.residue_class is ResidueClass.POLYMER and residue.name in standard


def primary_atom(atoms: tuple[Atom, ...]) -> Atom:
    """Highest-occupancy record; ties broken by altloc identifier (deterministic)."""
    return sorted(atoms, key=lambda atom: (-atom.occupancy, atom.altloc))[0]


def rules_of_kind(rules: tuple[Rule, ...], kind: str) -> tuple[Rule, ...]:
    return tuple(rule for rule in rules if rule.matcher["kind"] == kind)
