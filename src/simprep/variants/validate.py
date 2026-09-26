"""Check requested variants against the prepared wild-type system (pure; TASK-005).

Every problem is collected, so one run lists all of them.
"""

from __future__ import annotations

from dataclasses import dataclass

from simprep.prep.plan import PrepPlan
from simprep.structure.model import Residue, ResidueId, Structure

REQUIRED_BACKBONE = ("N", "CA", "C", "CB")
KEPT_ACTIONS = ("apply",)


class VariantError(ValueError):
    """The requested variants cannot be built as specified; the message lists every reason."""


@dataclass(frozen=True)
class Site:
    """One mutation of one variant, resolved against the prepared wild type."""

    variant: str
    residue: Residue
    to: str
    uniprot: dict | None

    @property
    def label(self) -> str:
        return f"{self.residue.id.label()} {self.residue.name}>{self.to}"


def mutation_id(mutation: dict) -> ResidueId:
    return ResidueId(mutation["chain"], mutation["seq_num"], mutation["ins_code"])


def resolve_sites(
    variants: list[dict], context: tuple[Structure, Structure, PrepPlan], side_chains: dict
) -> list[Site]:
    """Sites for every mutation; raise VariantError listing every problem.

    ``context`` is (input structure, prepared wild type, prep plan)."""
    source, wild_type, plan = context
    problems = _name_problems(variants)
    sites = []
    for variant in variants:
        for mutation in variant["mutations"]:
            found = _mutation_problems(mutation, (source, wild_type, plan), side_chains)
            problems += [f"{variant['name']}: {problem}" for problem in found]
            if not found:
                residue = wild_type.residue(mutation_id(mutation))
                sites.append(
                    Site(
                        variant["name"],
                        residue,
                        mutation["to"],
                        uniprot_position(source, residue.id),
                    )
                )
    if problems:
        raise VariantError("variants cannot be built:\n  " + "\n  ".join(problems))
    return sites


def _name_problems(variants: list[dict]) -> list[str]:
    problems, seen = [], set()
    for variant in variants:
        if variant["name"] in seen:
            problems.append(f"{variant['name']}: variant name used twice")
        seen.add(variant["name"])
        positions = [mutation_id(m) for m in variant["mutations"]]
        if len(set(positions)) != len(positions):
            problems.append(f"{variant['name']}: the same residue is mutated twice")
    return problems


def _mutation_problems(mutation: dict, context: tuple, side_chains: dict) -> list[str]:
    source, wild_type, plan = context
    rid = mutation_id(mutation)
    where = f"{rid.label()} {mutation['from']}>{mutation['to']}"
    if rid not in wild_type.residue_index:
        return [f"{where}: residue not in the prepared system (unobserved or excluded)"]
    residue = wild_type.residue(rid)
    problems = []
    if residue.name != mutation["from"]:
        problems.append(f"{where}: the structure has {residue.name} there")
    sequence_name = _sequence_name(source, residue)
    if sequence_name not in (None, mutation["from"]):
        problems.append(f"{where}: the entity sequence has {sequence_name} there")
    if mutation["to"] == mutation["from"]:
        problems.append(f"{where}: from and to are the same residue")
    if mutation["to"] not in side_chains:
        problems.append(
            f"{where}: no side-chain data for {mutation['to']} in "
            "knowledge/side_chains.yaml (add it with its source)"
        )
    problems += _residue_problems(residue, plan, where)
    return problems


def _residue_problems(residue: Residue, plan: PrepPlan, where: str) -> list[str]:
    problems = []
    names = {atom.name for atom in residue.atoms}
    missing = [name for name in REQUIRED_BACKBONE if name not in names]
    if missing:
        problems.append(f"{where}: missing {', '.join(missing)}")
    if residue.altlocs:
        problems.append(f"{where}: altlocs {', '.join(residue.altlocs)} are still undecided")
    for action in plan.actions:
        if residue.id in action.residues and action.prep_action not in KEPT_ACTIONS:
            problems.append(
                f"{where}: covered by {action.finding_id} "
                f"({action.prep_action}); decide it with an applied option first"
            )
    return problems


def _sequence_name(source: Structure, residue: Residue) -> str | None:
    sequence = source.sequence(residue.id.chain)
    if residue.label_seq is None or not sequence or residue.label_seq > len(sequence):
        return None
    return sequence[residue.label_seq - 1]


def uniprot_position(source: Structure, residue_id: ResidueId) -> dict | None:
    """UniProt accession and position of a residue from the file's own struct_ref_seq."""
    for reference in source.sequence_references:
        position = reference.db_position(residue_id.seq_num)
        if reference.chain == residue_id.chain and reference.db_name == "UNP" and position:
            return {"accession": reference.db_accession, "position": position}
    return None
