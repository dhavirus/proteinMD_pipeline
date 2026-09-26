"""Turn a decided manifest into a prep plan (pure: no files, no gemmi).

Each finding in the manifest's snapshot becomes one :class:`Action`. Options whose
``prep_action`` is ``apply`` are carried out by an operation registered in
:data:`OPERATIONS`, unless they name a ``prep_stage``: then that later stage carries them
out (``model_loop``: the modelling stage, ADR-0007) and prep only lists them. ``record``
and ``defer`` options, and apply options of a later stage, become work-order items;
findings without a final decision are passed through unchanged and listed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from simprep.decisions import decision_status
from simprep.rules import RuleSet
from simprep.structure.model import Link, LinkPartner, Residue, ResidueId, Structure

ADDED_LINK_TYPE = "covale"
METAL_LINK = "metalc"
TRUNCATE = "truncate"
# The fixed meaning of `truncate` for the topology stage (TASK-007): charged termini.
TERMINUS_NOTES = {
    "n_terminal": "charged N-terminus (NH3+) at {after}",
    "c_terminal": "charged C-terminus (COO-) at {before}",
    "internal": "chain break: charged C-terminus (COO-) at {before}, charged N-terminus "
    "(NH3+) at {after}",
}
ADDED_LINK_PREFIX = "prep_link"


class PrepError(ValueError):
    """The manifest cannot be applied as it stands; the message lists every reason."""


@dataclass(frozen=True)
class Action:
    """What prep does for one finding."""

    finding_id: str
    rule_id: str
    option_id: str | None
    prep_action: str  # apply | record | defer | undecided
    prep_stage: str | None
    residues: tuple[ResidueId, ...]
    note: str


@dataclass(frozen=True)
class WorkItem:
    stage: str
    finding_id: str
    option_id: str | None
    description: str
    residues: tuple[ResidueId, ...]


@dataclass
class PlanDraft:
    """Mutable accumulator the operations fill; frozen into :class:`PrepPlan`."""

    excluded: set[ResidueId] = field(default_factory=set)
    altloc_choice: dict[ResidueId, str] = field(default_factory=dict)
    ensemble: set[ResidueId] = field(default_factory=set)
    new_links: list[Link] = field(default_factory=list)
    reverts: dict[ResidueId, dict] = field(default_factory=dict)
    problems: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class PrepPlan:
    actions: tuple[Action, ...]
    work_order: tuple[WorkItem, ...]
    excluded: frozenset[ResidueId]
    altloc_choice: dict[ResidueId, str]
    ensemble_residues: frozenset[ResidueId]
    ensemble_altlocs: tuple[str, ...]
    new_links: tuple[Link, ...]
    reverts: dict[ResidueId, dict]


def build_plan(structure: Structure, manifest: dict, ruleset: RuleSet) -> PrepPlan:
    """The prep plan for ``structure`` under ``manifest``; raises PrepError listing every
    reason it cannot be applied (undecided blocking findings, conflicts, missing data)."""
    status = decision_status(manifest, include_variants=False)
    if status.blocking_undecided or status.unresolved_decisions:
        raise PrepError(_gate_message(status))
    findings = {f["id"]: f for f in manifest["findings_snapshot"]["findings"]}
    decisions = {d["finding_id"]: d for d in manifest["decisions"]}
    draft = PlanDraft()
    context = OperationContext(structure, findings, ruleset, draft)
    actions = [_action(context, finding, decisions.get(fid)) for fid, finding in findings.items()]
    ensemble_altlocs = _ensemble_altlocs(structure, draft)
    draft.problems += _conflicts(actions, draft)
    if draft.problems:
        raise PrepError("prep cannot apply this manifest:\n  " + "\n  ".join(draft.problems))
    return PrepPlan(
        actions=tuple(actions),
        work_order=tuple(item for a in actions if (item := _work_item(a, findings)) is not None),
        excluded=frozenset(draft.excluded),
        altloc_choice=dict(draft.altloc_choice),
        ensemble_residues=frozenset(draft.ensemble),
        ensemble_altlocs=ensemble_altlocs,
        new_links=tuple(draft.new_links),
        reverts=dict(draft.reverts),
    )


def _gate_message(status) -> str:
    lines = [
        f"[{f['effective_severity']}] {f['id']}"
        for f in status.undecided_findings
        if f["effective_severity"] == "blocking"
    ]
    lines += [
        f"{d['finding_id']}: {d['option_id']} is not a final decision"
        for d in status.unresolved_decisions
    ]
    return (
        "prep needs a final decision for every blocking finding, and no decision may use a "
        "non-final option (see `simprep manifest status`):\n  " + "\n  ".join(lines)
    )


def finding_residues(finding: dict) -> tuple[ResidueId, ...]:
    """Residues the finding covers (its extent), as ids."""
    extent = finding["locus"]["extent"]
    return tuple(ResidueId(r["chain"], r["seq_num"], r["ins_code"]) for r in extent)


def _action(context: OperationContext, finding: dict, decision: dict | None) -> Action:
    """The action for one finding; ``apply`` options run their operation on the draft."""

    def action(option_id, prep_action, prep_stage, note) -> Action:
        return Action(
            finding["id"],
            finding["rule_id"],
            option_id,
            prep_action,
            prep_stage,
            finding_residues(finding),
            note,
        )

    if decision is None:
        return action(None, "undecided", None, "no decision recorded; passed through unchanged")
    option = next(o for o in finding["options"] if o["id"] == decision["option_id"])
    prep_action = option.get("prep_action")
    if prep_action is None:
        context.draft.problems.append(
            f"{finding['id']}: option {option['id']} has no prep_action; the findings snapshot "
            "predates knowledge base 0.3.0, so re-run `simprep audit --manifest` first"
        )
        return action(option["id"], "undecided", None, "no prep_action in the snapshot")
    if prep_action == "apply" and "prep_stage" in option:
        stage = option["prep_stage"]
        return action(option["id"], "apply", stage, f"{option['label']} ({stage} stage)")
    if prep_action == "apply":
        note = OPERATIONS[option["id"]](context, finding, decision)
        return action(option["id"], "apply", None, note)
    return action(option["id"], prep_action, option.get("prep_stage"), option["label"])


def _work_item(action: Action, findings: dict) -> WorkItem | None:
    finding = findings[action.finding_id]
    if action.prep_action == "undecided":
        return WorkItem(
            "decision",
            action.finding_id,
            action.option_id,
            f"Decide: {finding['title']}",
            action.residues,
        )
    if action.prep_stage is not None:
        return WorkItem(
            action.prep_stage,
            action.finding_id,
            action.option_id,
            f"{action.note}: {finding['title']}{_terminus_note(action, finding)}",
            action.residues,
        )
    return None


def _terminus_note(action: Action, finding: dict) -> str:
    """For `truncate` on unobserved residues: which termini the topology stage builds."""
    position = _evidence(finding, "position")
    if action.option_id != TRUNCATE or position not in TERMINUS_NOTES:
        return ""
    before, after = (_flank(finding, side) for side in ("before", "after"))
    return "; " + TERMINUS_NOTES[position].format(before=before, after=after)


def _flank(finding: dict, side: str) -> str | None:
    """``THR A:34`` from the flank evidence (``THRA:34``) and the anchor residues."""
    value = _evidence(finding, f"flank_{side}")
    for ref in finding["anchor_residues"]:
        label = ResidueId(ref["chain"], ref["seq_num"], ref["ins_code"]).label()
        if value and value.endswith(label):
            return f"{value[: -len(label)]} {label}"
    return value


# ---------------------------------------------------------------- apply operations


@dataclass(frozen=True)
class OperationContext:
    structure: Structure
    findings: dict
    ruleset: RuleSet
    draft: PlanDraft


def _present(context: OperationContext, residues) -> list[Residue]:
    index = context.structure.residue_index
    return [index[rid] for rid in residues if rid in index]


def _highest_occupancy(context: OperationContext, finding: dict, decision: dict) -> str:
    chosen = []
    for residue in _present(context, finding_residues(finding)):
        if residue.altlocs:
            altloc = best_altloc(residue)
            context.draft.altloc_choice[residue.id] = altloc
            chosen.append(f"{residue.id.label()}:{altloc}")
    return "kept altloc " + ", ".join(chosen)


def best_altloc(residue: Residue) -> str:
    """Altloc with the highest mean occupancy; ties broken by altloc id."""

    def mean_occupancy(altloc: str) -> float:
        values = [a.occupancy for a in residue.atoms if a.altloc == altloc]
        return sum(values) / len(values)

    return sorted(residue.altlocs, key=lambda alt: (-mean_occupancy(alt), alt))[0]


def _specific_altloc(context: OperationContext, finding: dict, decision: dict) -> str:
    altloc = decision.get("parameters", {}).get("altloc")
    if not altloc:
        context.draft.problems.append(f"{finding['id']}: specific_altloc needs parameters.altloc")
        return "altloc missing"
    for residue in _present(context, finding_residues(finding)):
        if residue.altlocs and altloc not in residue.altlocs:
            context.draft.problems.append(
                f"{finding['id']}: {residue.id.label()} has no altloc {altloc!r} "
                f"(has {', '.join(residue.altlocs)})"
            )
        elif residue.altlocs:
            context.draft.altloc_choice[residue.id] = altloc
    return f"kept altloc {altloc}"


def _keep_ensemble(context: OperationContext, finding: dict, decision: dict) -> str:
    residues = [r for r in _present(context, finding_residues(finding)) if r.altlocs]
    context.draft.ensemble.update(r.id for r in residues)
    return "one system per altloc for " + ", ".join(r.id.label() for r in residues)


def _exclude(context: OperationContext, finding: dict, decision: dict) -> str:
    residues = [r.id for r in _present(context, finding_residues(finding))]
    context.draft.excluded.update(residues)
    return "excluded " + ", ".join(rid.label() for rid in residues)


def _exclude_group(context: OperationContext, finding: dict, decision: dict) -> str:
    """Exclude the ligand of a covalent-contact candidate, with its whole glycan group."""
    group_id = _evidence(finding, "ligand_group_finding")
    ligand = finding_residues(finding)[:1]
    group = finding_residues(context.findings[group_id]) if group_id in context.findings else ligand
    residues = [r.id for r in _present(context, group)]
    context.draft.excluded.update(residues)
    return "excluded " + ", ".join(rid.label() for rid in residues)


def _add_link(context: OperationContext, finding: dict, decision: dict) -> str:
    item = next(e for e in finding["evidence"] if e["key"] == "contact_distance")
    partners = [
        LinkPartner(
            ResidueId(a["chain"], a["seq_num"], a["ins_code"]),
            a["res_name"],
            a["atom_name"],
            a["altloc"],
        )
        for a in item["atoms"]
    ]
    conn_id = f"{ADDED_LINK_PREFIX}{len(context.draft.new_links) + 1}"
    if any(link.conn_id == conn_id for link in context.structure.links):
        context.draft.problems.append(f"{finding['id']}: input already has a link {conn_id}")
    link = Link(conn_id, ADDED_LINK_TYPE, partners[0], partners[1], item["value"])
    context.draft.new_links.append(link)
    ends = " - ".join(f"{p.res_name} {p.residue.label()} {p.atom_name}" for p in partners)
    return f"added {ADDED_LINK_TYPE} link {conn_id}: {ends}"


def _mapped(context: OperationContext, finding: dict, decision: dict) -> str:
    """Delete and rename atoms per the mapping for (component, option) in
    knowledge/residue_mappings.yaml (revert_to_parent, model_gem_diol)."""
    option = decision["option_id"]
    mappings = {(m["from"], m["option"]): m for m in context.ruleset.residue_mappings}
    notes = []
    for residue in _present(context, finding_residues(finding)):
        mapping = mappings.get((residue.name, option))
        if mapping is None:
            context.draft.problems.append(
                f"{finding['id']}: {option} needs an atom mapping for {residue.name} in "
                "knowledge/residue_mappings.yaml; none exists, so prep will not guess"
            )
            continue
        context.draft.reverts[residue.id] = mapping
        notes.append(_mapping_note(residue, mapping, context.structure.links))
    return "; ".join(notes)


def _mapping_note(residue: Residue, mapping: dict, links: tuple[Link, ...]) -> str:
    """What the mapping does to ``residue``, including links it removes; a removed metal
    link means the metal site's coordination changes (TASK-008 decision 4)."""
    renamed = ", ".join(f"{old}->{new}" for old, new in mapping["rename"].items())
    note = f"{residue.name} {residue.id.label()} -> {mapping['to']}"
    deleted = ", ".join(mapping["delete"]) or "nothing"
    note += f" (deleted {deleted}; renamed {renamed or 'nothing'})"
    removed = [
        (link, other)
        for link in links
        for mine, other in ((link.partner1, link.partner2), (link.partner2, link.partner1))
        if mine.residue == residue.id and mine.atom_name in mapping["delete"]
    ]
    for link, other in removed:
        note += f"; link {link.conn_id} to {other.res_name} {other.residue.label()} removed"
        if link.conn_type == METAL_LINK:
            note += " (metal coordination changed; the metal decision stands)"
    return note


def _evidence(finding: dict, key: str):
    return next((e.get("value") for e in finding["evidence"] if e["key"] == key), None)


OPERATIONS = {
    "highest_occupancy": _highest_occupancy,
    "specific_altloc": _specific_altloc,
    "keep_ensemble": _keep_ensemble,
    "exclude": _exclude,
    "exclude_segment": _exclude,
    "exclude_group": _exclude_group,
    "add_link": _add_link,
    "revert_to_parent": _mapped,
    "model_gem_diol": _mapped,
}


def _ensemble_altlocs(structure: Structure, draft: PlanDraft) -> tuple[str, ...]:
    """Altloc ids of the ensemble; every ensemble residue must have all of them."""
    residues = [structure.residue_index[rid] for rid in sorted(draft.ensemble)]
    altlocs = tuple(sorted({alt for r in residues for alt in r.altlocs}))
    for residue in residues:
        missing = sorted(set(altlocs) - set(residue.altlocs))
        if missing:
            draft.problems.append(
                f"keep_ensemble: {residue.id.label()} has no altloc {', '.join(missing)}; "
                "an ensemble system needs every ensemble residue in every altloc"
            )
    return altlocs


def _conflicts(actions: list[Action], draft: PlanDraft) -> list[str]:
    """Decisions that keep, link or defer residues another decision excludes."""
    problems = []
    link_residues = {p.residue for link in draft.new_links for p in (link.partner1, link.partner2)}
    for rid in sorted(draft.excluded & link_residues):
        problems.append(f"{rid.label()} is excluded by one decision and linked by add_link")
    problems += _link_altloc_conflicts(draft)
    for action in actions:
        if action.prep_action == "defer" and set(action.residues) & draft.excluded:
            clash = ", ".join(r.label() for r in sorted(set(action.residues) & draft.excluded))
            problems.append(
                f"{action.finding_id}: {action.option_id} keeps {clash}, "
                "which another decision excludes"
            )
    return problems


def _link_altloc_conflicts(draft: PlanDraft) -> list[str]:
    """An added bond to an altloc atom needs that altloc kept in every system."""
    problems = []
    for link in draft.new_links:
        for partner in (link.partner1, link.partner2):
            label = f"{link.conn_id} ({partner.residue.label()} {partner.atom_name})"
            if not partner.altloc:
                continue
            if partner.residue in draft.ensemble:
                problems.append(
                    f"{label} bonds altloc {partner.altloc}, but keep_ensemble "
                    "puts the residue in a system per altloc"
                )
            elif draft.altloc_choice.get(partner.residue, partner.altloc) != partner.altloc:
                problems.append(
                    f"{label} bonds altloc {partner.altloc}, but the altloc "
                    f"decision keeps {draft.altloc_choice[partner.residue]}"
                )
    return problems
