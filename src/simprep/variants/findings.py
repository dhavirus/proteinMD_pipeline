"""variant_build findings: one per mutation site, every rotamer candidate as evidence
(pure; TASK-005). The human chooses the rotamer; nothing is picked silently."""

from __future__ import annotations

from dataclasses import dataclass

from simprep.findings import Finding, atom_ref, ev_count, ev_flag, ev_label, residue_locus
from simprep.rules import RuleSet
from simprep.structure.model import Atom, Structure
from simprep.variants.build import Candidate, build_side_chain, candidates
from simprep.variants.check import ClashCriterion, Contact, clashes, contacts, neighbours
from simprep.variants.validate import Site

FAMILY = "variant_build"
DECIMALS = 3


@dataclass(frozen=True)
class Evaluation:
    """One candidate placed at a site, with its contacts and clashes."""

    candidate: Candidate
    atoms: tuple[Atom, ...]
    contacts: tuple[Contact, ...]
    clashes: tuple[Contact, ...]


def evaluate_site(site: Site, wild_type: Structure, ruleset: RuleSet) -> list[Evaluation]:
    """Every candidate of the target residue built at ``site`` and checked."""
    rule = ruleset.family(FAMILY)[0]
    criterion = ClashCriterion.from_matcher(rule.matcher)
    data = ruleset.side_chains[site.to]
    nearby = neighbours(wild_type, site.residue.id)
    evaluations = []
    for candidate in candidates(data):
        atoms = build_side_chain(site.residue, data, candidate.chi_degree)
        found = contacts(atoms, nearby, criterion)
        evaluations.append(
            Evaluation(candidate, atoms, tuple(found), tuple(clashes(found, criterion)))
        )
    return evaluations


def best(evaluations: list[Evaluation]) -> Evaluation:
    """Fewest clashes, then highest library frequency, then rotamer id (deterministic)."""
    return min(
        evaluations,
        key=lambda e: (len(e.clashes), -e.candidate.frequency_percent, e.candidate.rotamer_id),
    )


def finding_id(site: Site) -> str:
    return f"{FAMILY}/{site.variant}/{site.residue.id.label()}"


def site_finding(site: Site, evaluations: list[Evaluation], ruleset: RuleSet) -> Finding:
    """The rotamer-choice finding for ``site`` (severity context not yet applied)."""
    rule = ruleset.family(FAMILY)[0]
    chosen = best(evaluations)
    evidence = [
        ev_label("mutation", site.label),
        ev_label("uniprot_position", _uniprot_label(site)),
        ev_count("candidates", len(evaluations)),
        *(_candidate_item(site, e) for e in evaluations),
        ev_label(
            "recommended_rotamer",
            chosen.candidate.rotamer_id,
            "fewest clashes, then library frequency",
        ),
        ev_flag("all_candidates_clash", all(e.clashes for e in evaluations)),
    ]
    return Finding(
        id=finding_id(site),
        rule=rule,
        title=f"{site.variant}: side chain for {site.label}",
        locus=residue_locus(site.residue),
        anchor_residues=(site.residue.id,),
        evidence=tuple(evidence),
    )


def _uniprot_label(site: Site) -> str:
    if site.uniprot is None:
        return "not mapped (no UniProt struct_ref_seq range covers this residue)"
    return f"{site.uniprot['accession']} {site.uniprot['position']}"


def _candidate_item(site: Site, evaluation: Evaluation) -> dict:
    candidate, worst = evaluation.candidate, (evaluation.contacts or (None,))[0]
    item = {
        "key": "candidate",
        "type": "candidate",
        "value": candidate.rotamer_id,
        "chi_degree": [round(chi, 1) for chi in candidate.chi_degree],
        "frequency_percent": candidate.frequency_percent,
        "clash_count": len(evaluation.clashes),
        "max_overlap_angstrom": None if worst is None else round(worst.overlap_angstrom, DECIMALS),
    }
    if worst is not None:
        item["atoms"] = [
            atom_ref(site.residue, worst.new_atom) | {"res_name": site.to},
            atom_ref(worst.residue, worst.atom),
        ]
        item["note"] = (
            f"closest: {worst.new_atom.name} - {worst.residue.name} "
            f"{worst.residue.id.label()} {worst.atom.name} "
            f"{worst.distance_angstrom:.2f} A"
        )
    return item
