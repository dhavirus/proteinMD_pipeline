"""Everything no family rule claimed: unknown chemistry groups and unexplained links.

Runs after the family detectors and receives the set of residues they claimed.
Standard polymer residues and waters need no claim.
"""

from __future__ import annotations

from simprep.detectors.common import is_standard_polymer, rules_of_kind, standard_polymer_residues
from simprep.findings import (
    Finding,
    Locus,
    ev_b_factor,
    ev_count,
    ev_label,
    ev_occupancy,
    ev_source_record,
    partner_ref,
)
from simprep.rules import Rule, RuleSet
from simprep.structure.model import Link, Residue, ResidueClass, ResidueId, Structure

WILDCARD = "*"


def detect_unrecognized(
    structure: Structure, ruleset: RuleSet, claimed: frozenset[ResidueId]
) -> list[Finding]:
    """Blocking findings for unclaimed chemistry groups and unexplained covalent links."""
    rules = ruleset.family("unrecognized")
    (group_rule,) = rules_of_kind(rules, "unclaimed_group")
    (link_rule,) = rules_of_kind(rules, "unexplained_link")
    unclaimed = unclaimed_residues(structure, ruleset, (claimed, group_rule))
    groups = connected_groups(unclaimed, structure.links, group_rule.matcher["group_by_conn_types"])
    grouped = frozenset(residue.id for group in groups for residue in group)
    findings = [_group_finding(structure, group, group_rule) for group in groups]
    findings += [
        _link_finding(link, link_rule)
        for link in structure.links
        if _is_unexplained(structure, link, (link_rule, grouped))
    ]
    return findings


def unclaimed_residues(
    structure: Structure, ruleset: RuleSet, claimed_rule: tuple[frozenset[ResidueId], Rule]
) -> list[Residue]:
    """Residues of the rule's classes, other than standard polymer residues and waters,
    that no family rule claimed."""
    claimed, rule = claimed_rule
    standard = standard_polymer_residues(ruleset)
    classes = set(rule.matcher["residue_classes"])
    return [
        residue
        for residue in structure.residues
        if residue.id not in claimed
        and residue.residue_class is not ResidueClass.WATER
        and not is_standard_polymer(residue, standard)
        and residue.residue_class.value in classes
    ]


def connected_groups(
    residues: list[Residue], links: tuple[Link, ...], conn_types: list[str]
) -> list[tuple[Residue, ...]]:
    """Group residues joined (transitively) by links of ``conn_types``; deterministic order."""
    parent = {residue.id: residue.id for residue in residues}

    def root(rid: ResidueId) -> ResidueId:
        while parent[rid] != rid:
            rid = parent[rid]
        return rid

    for link in links:
        a, b = link.partner1.residue, link.partner2.residue
        if link.conn_type in conn_types and a in parent and b in parent:
            parent[max(root(a), root(b))] = min(root(a), root(b))
    groups: dict[ResidueId, list[Residue]] = {}
    for residue in residues:
        groups.setdefault(root(residue.id), []).append(residue)
    return [tuple(sorted(members, key=lambda r: r.id)) for _, members in sorted(groups.items())]


def _group_finding(structure: Structure, group: tuple[Residue, ...], rule: Rule) -> Finding:
    ids = {residue.id for residue in group}
    attachments = [
        link
        for link in structure.links
        if (link.partner1.residue in ids) != (link.partner2.residue in ids)
    ]
    internal = [
        link
        for link in structure.links
        if link.partner1.residue in ids and link.partner2.residue in ids
    ]
    atoms = [atom for residue in group for atom in residue.atoms]
    names = "-".join(residue.name for residue in group)
    first = group[0]
    evidence = [
        ev_label("components", names),
        ev_label("residue_classes", ",".join(sorted({r.residue_class.value for r in group}))),
        ev_count("residues", len(group)),
        ev_count("atom_records", len(atoms)),
        ev_count("internal_links", len(internal)),
        ev_occupancy("mean_occupancy", sum(a.occupancy for a in atoms) / len(atoms)),
        ev_b_factor("mean_b_factor", sum(a.b_iso for a in atoms) / len(atoms)),
        ev_label("decision_required", "expert decision required"),
    ]
    evidence += [_link_record("attachment", link) for link in attachments]
    evidence += [_link_record("internal_link", link) for link in internal]
    evidence += _site_annotations(structure, attachments, ids)
    return Finding(
        id=f"unrecognized/group/{first.id.label()}",
        rule=rule,
        title=f"{rule.title}: {names} at {first.id.label()}"
        + (f" ({len(attachments)} external link(s))" if attachments else ""),
        locus=Locus(first.id, first.name, (), tuple((r.id, r.name) for r in group)),
        anchor_residues=tuple(r.id for r in group),
        evidence=tuple(evidence),
        claims=frozenset(ids),
    )


def _site_annotations(structure: Structure, attachments: list[Link], ids: set) -> list[dict]:
    """Modification records on the residues this group attaches to (e.g. glycosylation)."""
    sites = {
        p.residue
        for link in attachments
        for p in (link.partner1, link.partner2)
        if p.residue not in ids
    }
    return [
        ev_source_record(
            "attachment_site_annotation",
            "pdbx_struct_mod_residue",
            {
                "residue": mod.residue.label(),
                "res_name": mod.res_name,
                "parent_comp_id": mod.parent_res_name,
                "details": mod.details,
            },
        )
        for mod in structure.modified_residues
        if mod.residue in sites
    ]


def _is_unexplained(structure: Structure, link: Link, rule_grouped: tuple) -> bool:
    rule, grouped = rule_grouped
    if link.conn_type not in rule.matcher["conn_types"]:
        return False
    if link.partner1.residue in grouped or link.partner2.residue in grouped:
        return False
    return not any(
        _recognized(structure, link, pattern) for pattern in rule.matcher["recognized_links"]
    )


def _recognized(structure: Structure, link: Link, pattern: dict) -> bool:
    """True if ``link`` matches ``pattern`` with partners in either order."""
    if link.conn_type != pattern["conn_type"]:
        return False
    partners = (link.partner1, link.partner2)
    return any(
        all(
            _partner_matches(structure, partner, (pattern, index))
            for index, partner in enumerate(order)
        )
        for order in (partners, partners[::-1])
    )


def _partner_matches(structure: Structure, partner, pattern_index: tuple[dict, int]) -> bool:
    pattern, index = pattern_index
    name = pattern["res_names"][index]
    if name != WILDCARD and partner.res_name != name:
        return False
    if partner.atom_name != pattern["atom_names"][index]:
        return False
    classes = pattern.get("residue_classes")
    if classes is None:
        return True
    residue = structure.residue_index.get(partner.residue)
    return residue is not None and residue.residue_class.value == classes[index]


def _link_record(key: str, link: Link) -> dict:
    return ev_source_record(
        key,
        "struct_conn",
        {
            "id": link.conn_id,
            "conn_type": link.conn_type,
            "distance_angstrom": link.distance_angstrom,
        },
    ) | {"atoms": [partner_ref(link.partner1), partner_ref(link.partner2)]}


def _link_finding(link: Link, rule: Rule) -> Finding:
    p1, p2 = link.partner1, link.partner2
    return Finding(
        id=f"unrecognized/link/{link.conn_id}",
        rule=rule,
        title=f"{rule.title}: {p1.res_name}{p1.residue.label()}.{p1.atom_name} - "
        f"{p2.res_name}{p2.residue.label()}.{p2.atom_name} ({link.conn_type})",
        locus=Locus(
            p1.residue,
            p1.res_name,
            (p1.atom_name,),
            ((p1.residue, p1.res_name), (p2.residue, p2.res_name)),
        ),
        anchor_residues=(p1.residue, p2.residue),
        evidence=(
            _link_record("link", link),
            ev_label("decision_required", "expert decision required"),
        ),
    )
