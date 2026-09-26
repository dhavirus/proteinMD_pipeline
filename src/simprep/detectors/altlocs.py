"""Alternative locations: occupancies, B-factors and spatial spread per residue."""

from __future__ import annotations

import itertools
from collections import defaultdict

from simprep.config import AuditConfig
from simprep.detectors import register
from simprep.findings import (
    Finding,
    Locus,
    atom_ref,
    ev_b_factor,
    ev_count,
    ev_distance,
    ev_label,
    ev_occupancy,
    residue_locus,
)
from simprep.rules import Rule, RuleSet
from simprep.structure.geometry import distance
from simprep.structure.model import Atom, Residue, Structure


@register("altlocs")
def detect_altlocs(structure: Structure, ruleset: RuleSet, config: AuditConfig) -> list[Finding]:
    """Per-residue findings, plus one finding per aggregating rule (e.g. waters)."""
    rules = ruleset.family("altlocs")
    findings: list[Finding] = []
    aggregated: dict[str, list[Residue]] = defaultdict(list)
    for residue in structure.residues:
        if not residue.altlocs:
            continue
        rule = _matching_rule(residue, rules)
        if rule is None:
            continue
        if rule.matcher["aggregate"]:
            aggregated[rule.rule_id].append(residue)
        else:
            findings.append(_residue_finding(residue, rule))
    by_id = {rule.rule_id: rule for rule in rules}
    findings += [_aggregate_finding(members, by_id[rid]) for rid, members in aggregated.items()]
    return findings


def max_spread(residue: Residue) -> tuple[float, tuple[Atom, Atom] | None]:
    """Largest distance between altloc records of the same atom name."""
    by_name: dict[str, list[Atom]] = defaultdict(list)
    for atom in residue.atoms:
        by_name[atom.name].append(atom)
    pairs = [pair for atoms in by_name.values() for pair in itertools.combinations(atoms, 2)]
    if not pairs:
        return 0.0, None
    best = max(
        pairs,
        key=lambda pair: (
            distance(pair[0].position, pair[1].position),
            pair[0].name,
            pair[0].altloc,
            pair[1].altloc,
        ),
    )
    return distance(best[0].position, best[1].position), best


def _matching_rule(residue: Residue, rules: tuple[Rule, ...]) -> Rule | None:
    spread, _ = max_spread(residue)
    for rule in rules:
        matcher = rule.matcher
        upper = matcher["max_spread_angstrom"]
        if (
            residue.residue_class.value in matcher["residue_classes"]
            and spread >= matcher["min_spread_angstrom"]
            and (upper is None or spread < upper)
        ):
            return rule
    return None


def _residue_finding(residue: Residue, rule: Rule) -> Finding:
    spread, pair = max_spread(residue)
    evidence = [
        ev_label("altloc_ids", ",".join(residue.altlocs)),
        ev_count("atoms_with_altlocs", len({a.name for a in residue.atoms if a.altloc})),
    ]
    for altloc in residue.altlocs:
        atoms = [atom for atom in residue.atoms if atom.altloc == altloc]
        evidence.append(
            ev_occupancy(
                f"occupancy_{_key(altloc)}",
                _mean([a.occupancy for a in atoms]),
                note=f"mean over altloc {altloc}",
            )
        )
        evidence.append(
            ev_b_factor(
                f"b_factor_{_key(altloc)}",
                _mean([a.b_iso for a in atoms]),
                note=f"mean over altloc {altloc}",
            )
        )
    if pair is not None:
        evidence.append(
            ev_distance(
                "max_spread", spread, [atom_ref(residue, pair[0]), atom_ref(residue, pair[1])]
            )
        )
    return Finding(
        id=f"altlocs/{residue.id.label()}",
        rule=rule,
        title=f"{rule.title}: {residue.name} {residue.id.label()} ({spread:.2f} A)",
        locus=residue_locus(residue, tuple(sorted({a.name for a in residue.atoms if a.altloc}))),
        anchor_residues=(residue.id,),
        evidence=tuple(evidence),
    )


def _aggregate_finding(residues: list[Residue], rule: Rule) -> Finding:
    spreads = [max_spread(residue)[0] for residue in residues]
    return Finding(
        id=f"altlocs/{rule.rule_id.split('.', 1)[1]}",
        rule=rule,
        title=f"{rule.title}: {len(residues)} residues",
        locus=Locus(None, None, (), tuple((r.id, r.name) for r in residues)),
        anchor_residues=tuple(r.id for r in residues),
        evidence=(
            ev_count("residues_with_altlocs", len(residues)),
            ev_label("residues", " ".join(f"{r.name}{r.id.label()}" for r in residues)),
            ev_distance("largest_spread", max(spreads), [], note="over all listed residues"),
        ),
    )


def _key(altloc: str) -> str:
    return altloc.lower() if altloc.isalnum() else f"x{ord(altloc)}"


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)
