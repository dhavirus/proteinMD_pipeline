"""The prep record (prep_record.json) and its short human-readable report (pure).

Counts compare each output system with the input, per record type:
``in = passed_through + modified + excluded`` and ``out = passed_through + modified + added``.
An atom record counts as modified when its residue was modified (altlocs collapsed,
renamed by a mapping); atoms a modification dropped count as excluded.
"""

from __future__ import annotations

from dataclasses import dataclass

from simprep import SCHEMA_VERSION
from simprep.prep.apply import System
from simprep.prep.plan import Action, PrepPlan, WorkItem
from simprep.rules import RuleSet
from simprep.structure.model import ResidueClass, ResidueId, Structure


@dataclass(frozen=True)
class RecordContext:
    """Provenance stamped onto the prep record."""

    input: dict
    manifest_sha256: str
    ruleset: RuleSet
    simprep: dict
    generated_at: str


def system_entry(structure: Structure, system: System, files: list[dict]) -> dict:
    """One ``systems`` entry: the written files and the accounting against ``structure``."""
    return {
        "name": system.name,
        "altloc": system.altloc,
        "files": files,
        "counts": system_counts(structure, system.structure),
    }


def build_prep_record(plan: PrepPlan, systems: list[dict], context: RecordContext) -> dict:
    """The prep_record.json document (validate with the prep_record schema)."""
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": context.generated_at,
        "simprep": context.simprep,
        "knowledge_base": {"version": context.ruleset.version, "sha256": context.ruleset.sha256},
        "input": context.input,
        "manifest_sha256": context.manifest_sha256,
        "actions": [_action(a) for a in plan.actions],
        "work_order": [_work_item(w) for w in plan.work_order],
        "systems": systems,
    }


def _residues(ids: tuple[ResidueId, ...]) -> list[dict]:
    return [rid.to_dict() for rid in ids]


def _action(action: Action) -> dict:
    return {
        "finding_id": action.finding_id,
        "rule_id": action.rule_id,
        "option_id": action.option_id,
        "prep_action": action.prep_action,
        "prep_stage": action.prep_stage,
        "residues": _residues(action.residues),
        "note": action.note,
    }


def _work_item(item: WorkItem) -> dict:
    return {
        "stage": item.stage,
        "finding_id": item.finding_id,
        "option_id": item.option_id,
        "description": item.description,
        "residues": _residues(item.residues),
    }


def system_counts(before: Structure, after: Structure) -> list[dict]:
    """Record accounting of ``after`` (a prepared system) against ``before`` (the input)."""
    rows = [_residue_row(before, after, cls) for cls in ResidueClass]
    rows.append(_atom_row(before, after))
    rows.append(_keyed_row("struct_conn records", _links(before), _links(after)))
    rows.append(_keyed_row("pdbx_struct_mod_residue records", _mods(before), _mods(after)))
    rows.append(
        _keyed_row(
            "unobserved polymer residues (annotation)", _unobserved(before), _unobserved(after)
        )
    )
    return rows


def _residue_row(before: Structure, after: Structure, cls: ResidueClass) -> dict:
    def of_class(structure: Structure) -> dict:
        return {r.id: r for r in structure.residues if r.residue_class is cls}

    return _keyed_row(f"{cls.value} residues", of_class(before), of_class(after))


def _atom_row(before: Structure, after: Structure) -> dict:
    out = after.residue_index
    row = dict.fromkeys(("passed_through", "modified", "excluded"), 0)
    for residue in before.residues:
        kept = out.get(residue.id)
        kept_atoms = 0 if kept is None else len(kept.atoms)
        row["excluded"] += len(residue.atoms) - kept_atoms
        row["passed_through" if kept == residue else "modified"] += kept_atoms
    total_in = sum(len(r.atoms) for r in before.residues)
    total_out = sum(len(r.atoms) for r in after.residues)
    return _row("atom records", (total_in, total_out), row)


def _keyed_row(record_type: str, before: dict, after: dict) -> dict:
    kept = before.keys() & after.keys()
    row = {
        "passed_through": sum(1 for key in kept if before[key] == after[key]),
        "modified": sum(1 for key in kept if before[key] != after[key]),
        "excluded": len(before.keys() - after.keys()),
    }
    return _row(record_type, (len(before), len(after)), row)


def _row(record_type: str, totals: tuple[int, int], row: dict) -> dict:
    total_in, total_out = totals
    added = total_out - row["passed_through"] - row["modified"]
    return {"record_type": record_type, "in": total_in, **row, "added": added, "out": total_out}


def _links(structure: Structure) -> dict:
    return {link.conn_id: link for link in structure.links}


def _mods(structure: Structure) -> dict:
    return {(m.residue, m.res_name): m for m in structure.modified_residues}


def _unobserved(structure: Structure) -> dict:
    return {u.residue: u for u in structure.unobserved_residues if u.is_polymer}


def render_prep_report(record: dict) -> str:
    """Short Markdown summary: systems and counts, actions, work order."""
    lines = [
        f"# Prep report: {record['input']['path']}",
        "",
        f"- input SHA-256 `{record['input']['sha256']}`",
        f"- manifest SHA-256 `{record['manifest_sha256']}`",
        f"- knowledge base {record['knowledge_base']['version']}, "
        f"simprep {record['simprep']['version']} ({record['simprep']['git_commit']})",
        "",
    ]
    for system in record["systems"]:
        lines += _system_lines(system)
    lines += ["## Actions", "", "| finding | option | prep action | note |", "|---|---|---|---|"]
    lines += [
        f"| {a['finding_id']} | {a['option_id'] or '-'} | {a['prep_action']} | {a['note']} |"
        for a in record["actions"]
    ]
    lines += ["", "## Work order for later stages", ""]
    lines += [
        f"- **{w['stage']}**: {w['description']} ({w['finding_id']})" for w in record["work_order"]
    ] or ["- nothing"]
    return "\n".join(lines) + "\n"


def _system_lines(system: dict) -> list[str]:
    files = ", ".join(f"`{f['path']}`" for f in system["files"])
    lines = [
        f"## {system['name']}",
        "",
        f"Files: {files}",
        "",
        "| records | in | passed | modified | excluded | added | out |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    lines += [
        f"| {c['record_type']} | {c['in']} | {c['passed_through']} | {c['modified']} | "
        f"{c['excluded']} | {c['added']} | {c['out']} |"
        for c in system["counts"]
    ]
    return [*lines, ""]
