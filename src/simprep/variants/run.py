"""`simprep variants` at the I/O edge (TASK-005).

Pass 1 (variant snapshot missing or stale): evaluate every rotamer candidate, write
``manifest.json`` with the variant_build findings and stop for decisions. Pass 2 (snapshot
current): check the rotamer decisions and write ``wt/`` (the prepared wild type, with its
prep record), one directory per variant, ``variant_record.json`` and ``variant_report.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from simprep.canonical import sha256_canonical
from simprep.manifest import attach_variant_snapshot, config_from_manifest, write_json
from simprep.prep.apply import apply_plan
from simprep.prep.plan import build_plan
from simprep.prep.record import system_counts
from simprep.prep.run import RECORD_FILE as PREP_RECORD_FILE
from simprep.prep.run import PrepInputs, PrepRequest, load_inputs, write_prep
from simprep.prep.write import write_system
from simprep.provenance import sha256_file, simprep_provenance, utc_now
from simprep.schemas import validate
from simprep.severity import apply_context
from simprep.structure.model import ResidueId, Structure
from simprep.variants.apply import apply_mutations, mutated_residue
from simprep.variants.findings import Evaluation, evaluate_site, finding_id, site_finding
from simprep.variants.relaxation import (
    Relaxed,
    protonation_items,
    relax_system,
    relaxation_entry,
    residual_findings,
    site_label,
    union_shell,
)
from simprep.variants.report import render_variant_report
from simprep.variants.validate import Site, VariantError, resolve_sites

MANIFEST_FILE = "manifest.json"
RECORD_FILE = "variant_record.json"
REPORT_FILE = "variant_report.md"
WILD_TYPE_DIR = "wt"
UNRELAXED_DIR = "unrelaxed"
SYSTEM_NAME = "system"
CONTACTS_RECORDED = 5
CHOOSE_ROTAMER = "choose_rotamer"
APPLIED_OPTIONS = (CHOOSE_ROTAMER,)  # apply options of the variant_build family


class VariantDecisionError(VariantError):
    """Rotamer decisions are missing, not final, or name an unknown rotamer."""


@dataclass(frozen=True)
class VariantOutcome:
    """``needs_decisions`` with the refreshed manifest, or ``built`` with the record."""

    status: str
    manifest_path: Path | None = None
    record: dict | None = None


@dataclass(frozen=True)
class SiteResult:
    site: Site
    evaluations: list[Evaluation]


def run_variants(request: PrepRequest) -> VariantOutcome:
    inputs = load_inputs(request)
    variants = inputs.manifest.get("variants") or []
    if not variants:
        raise VariantError("the manifest lists no variants; add them under `variants` first")
    plan = build_plan(inputs.structure, inputs.manifest, inputs.ruleset)
    wild_type = _wild_type(inputs, plan)
    context = (inputs.structure, wild_type, plan)
    sites = resolve_sites(variants, context, inputs.ruleset.side_chains)
    results = {
        finding_id(s): SiteResult(s, evaluate_site(s, wild_type, inputs.ruleset)) for s in sites
    }
    findings = _finding_dicts(inputs, wild_type, results)
    if not _snapshot_current(inputs.manifest, findings):
        return _needs_decisions(inputs, findings, request.out_dir)
    choices = chosen_candidates(inputs.manifest, results)
    wt_record = write_prep(inputs, plan, request.out_dir / WILD_TYPE_DIR)
    built = [_build(variant, (wild_type, results, choices)) for variant in variants]
    protocol = inputs.manifest.get("relaxation") or inputs.ruleset.relaxation
    context = Context(inputs, wild_type, wt_record, protocol, request.out_dir)
    relaxed_wts, relaxed = _relax_all(built, context) if protocol["enabled"] else ([], {})
    entries = [_variant_entry(b, relaxed.get(b.name), context) for b in built]
    record = _record(context, entries, relaxed_wts)
    return VariantOutcome("built", record=record)


def _wild_type(inputs: PrepInputs, plan) -> Structure:
    systems = apply_plan(inputs.structure, plan)
    if len(systems) != 1:
        raise VariantError(
            "variants v0.1 build from one wild-type system; this manifest keeps an altloc "
            f"ensemble ({', '.join(s.name for s in systems)}). Choose one altloc first."
        )
    return systems[0].structure


def _finding_dicts(inputs: PrepInputs, wild_type: Structure, results: dict) -> list[dict]:
    raw = [site_finding(r.site, r.evaluations, inputs.ruleset) for r in results.values()]
    config = config_from_manifest(inputs.manifest, inputs.manifest_sha256)
    return [finding.to_dict() for finding in apply_context(raw, wild_type, config)]


def _snapshot_current(manifest: dict, findings: list[dict]) -> bool:
    snapshot = manifest.get("variant_snapshot")
    return (
        snapshot is not None
        and snapshot["variants_sha256"] == sha256_canonical(manifest["variants"])
        and snapshot["findings_sha256"] == sha256_canonical(findings)
    )


def _needs_decisions(inputs: PrepInputs, findings: list[dict], out_dir: Path) -> VariantOutcome:
    """Write the manifest with fresh candidates (and the default relaxation protocol when
    it has none, so the study sees and can edit it)."""
    manifest = {"relaxation": inputs.ruleset.relaxation, **inputs.manifest}
    updated = attach_variant_snapshot(manifest, findings, utc_now())
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / MANIFEST_FILE
    write_json(updated, path)
    return VariantOutcome("needs_decisions", manifest_path=path)


def chosen_candidates(manifest: dict, results: dict) -> dict[str, Evaluation]:
    """The decided candidate per site; raise VariantDecisionError listing every problem."""
    decisions = {d["finding_id"]: d for d in manifest["decisions"]}
    chosen, problems = {}, []
    for fid, result in results.items():
        evaluation, problem = _choice(fid, decisions.get(fid), result.evaluations)
        if problem:
            problems.append(problem)
        else:
            chosen[fid] = evaluation
    if problems:
        raise VariantDecisionError("rotamer decisions needed:\n  " + "\n  ".join(problems))
    return chosen


def _choice(fid: str, decision: dict | None, evaluations: list[Evaluation]):
    if decision is None:
        return None, f"{fid}: undecided"
    if decision["option_id"] != CHOOSE_ROTAMER:
        return None, f"{fid}: {decision['option_id']} is not a final decision"
    rotamer = decision.get("parameters", {}).get("rotamer")
    by_id = {e.candidate.rotamer_id: e for e in evaluations}
    if rotamer not in by_id:
        return (
            None,
            f"{fid}: parameters.rotamer must be one of {', '.join(by_id)} (got {rotamer!r})",
        )
    return by_id[rotamer], None


@dataclass(frozen=True)
class Built:
    """A variant built rigidly on the prepared wild type."""

    name: str
    sites: frozenset[ResidueId]
    structure: Structure
    mutations: list[dict]
    protonation: list[dict]


@dataclass(frozen=True)
class Context:
    inputs: PrepInputs
    wild_type: Structure
    wt_record: dict
    protocol: dict
    out_dir: Path


def _build(variant: dict, built: tuple) -> Built:
    wild_type, results, choices = built
    fids = [fid for fid, r in results.items() if r.site.variant == variant["name"]]
    residues = tuple(
        mutated_residue(results[fid].site.residue, results[fid].site.to, choices[fid].atoms)
        for fid in fids
    )
    return Built(
        variant["name"],
        frozenset(results[fid].site.residue.id for fid in fids),
        apply_mutations(wild_type, residues),
        [_mutation_entry(fid, results[fid].site, choices[fid]) for fid in fids],
        [_protonation_item(fid, results[fid].site) for fid in fids],
    )


def _relax_all(built: list[Built], context: Context) -> tuple[list[dict], dict]:
    """Per site: the relaxed wild type (record entries) and each relaxed variant."""
    groups: dict[frozenset, list[Built]] = {}
    for variant in built:
        groups.setdefault(variant.sites, []).append(variant)
    wild_types, relaxed = [], {}
    for sites, members in sorted(groups.items(), key=lambda item: sorted(item[0])):
        structures = [context.wild_type] + [m.structure for m in members]
        where = (sites, union_shell(structures, sites, context.protocol))
        wild_types.append(_relaxed_wild_type(sites, where, context))
        for member in members:
            relaxed[member.name] = relax_system(
                member.name, member.structure, where, (context.protocol, context.inputs.ruleset)
            )
    return wild_types, relaxed


def _relaxed_wild_type(sites: frozenset, where: tuple, context: Context) -> dict:
    name = f"{WILD_TYPE_DIR}_relaxed_{site_label(sites)}"
    result = relax_system(
        name, context.wild_type, where, (context.protocol, context.inputs.ruleset)
    )
    structure = result.outcome.structure
    return {
        "name": name,
        "directory": name,
        "sites": [rid.to_dict() for rid in sorted(sites)],
        "files": _write(structure, context.out_dir / name),
        "counts": system_counts(context.wild_type, structure),
        "relaxation": relaxation_entry(result),
        "findings": _residual(result, context),
        "work_order": context.wt_record["work_order"]
        + protonation_items(result, f"relaxation/{name}", frozenset()),
    }


def _variant_entry(variant: Built, result: Relaxed | None, context: Context) -> dict:
    directory = context.out_dir / variant.name
    structure = variant.structure if result is None else result.outcome.structure
    entry = {
        "name": variant.name,
        "directory": variant.name,
        "mutations": variant.mutations,
        "files": _write(structure, directory),
        "counts": system_counts(context.wild_type, structure),
        "work_order": context.wt_record["work_order"] + variant.protonation,
    }
    if result is not None:
        unrelaxed = _write(variant.structure, directory / UNRELAXED_DIR)
        entry["unrelaxed_files"] = [
            {**f, "path": f"{UNRELAXED_DIR}/{f['path']}"} for f in unrelaxed
        ]
        entry["relaxation"] = relaxation_entry(result)
        entry["findings"] = _residual(result, context)
        entry["work_order"] += protonation_items(
            result, variant.mutations[0]["finding_id"], variant.sites
        )
    return entry


def _write(structure: Structure, directory: Path) -> list[dict]:
    directory.mkdir(parents=True, exist_ok=True)
    return write_system(structure, directory, SYSTEM_NAME)


def _residual(result: Relaxed, context: Context) -> list[dict]:
    config = config_from_manifest(context.inputs.manifest, context.inputs.manifest_sha256)
    return residual_findings(result, result.outcome.structure, (context.inputs.ruleset, config))


def _mutation_entry(fid: str, site: Site, evaluation: Evaluation) -> dict:
    residue = site.residue
    return {
        "finding_id": fid,
        **residue.id.to_dict(),
        "from": residue.name,
        "to": site.to,
        "uniprot": site.uniprot,
        "rotamer": evaluation.candidate.rotamer_id,
        "chi_degree": list(evaluation.candidate.chi_degree),
        "frequency_percent": evaluation.candidate.frequency_percent,
        "closest_contacts": [
            {
                "atom_name": c.new_atom.name,
                "partner": {
                    **c.residue.id.to_dict(),
                    "res_name": c.residue.name,
                    "atom_name": c.atom.name,
                    "altloc": c.atom.altloc,
                },
                "distance_angstrom": round(c.distance_angstrom, 3),
                "overlap_angstrom": round(c.overlap_angstrom, 3),
            }
            for c in evaluation.contacts[:CONTACTS_RECORDED]
        ],
    }


def _protonation_item(fid: str, site: Site) -> dict:
    return {
        "stage": "protonation",
        "finding_id": fid,
        "option_id": CHOOSE_ROTAMER,
        "description": f"{site.label}: the mutated residue has no hydrogens (TASK-005 decision 2)",
        "residues": [site.residue.id.to_dict()],
    }


def _record(context: Context, entries: list[dict], relaxed_wild_types: list[dict]) -> dict:
    inputs, wt_record, out_dir = context.inputs, context.wt_record, context.out_dir
    record = {
        "schema_version": wt_record["schema_version"],
        "generated_at": utc_now(),
        "simprep": simprep_provenance(),
        "knowledge_base": wt_record["knowledge_base"],
        "input": inputs.source,
        "manifest_sha256": inputs.manifest_sha256,
        "relaxation_protocol": context.protocol,
        "wild_type": {
            "directory": WILD_TYPE_DIR,
            "prep_record_sha256": sha256_file(out_dir / WILD_TYPE_DIR / PREP_RECORD_FILE),
        },
        "relaxed_wild_types": relaxed_wild_types,
        "variants": entries,
    }
    validate(record, "variant_record")
    write_json(record, out_dir / RECORD_FILE)
    (out_dir / REPORT_FILE).write_text(render_variant_report(record))
    return record
