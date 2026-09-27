"""The protonation stage at the I/O edge (TASK-009, ADR-0009).

prep -> model -> variants -> relax -> protonate. Each system gets its own PROPKA
estimates (TASK-009 decision 5), states from states.py, and every hydrogen from OpenMM.
``simprep protonate`` protonates the wild type (prepared, or modelled when gaps are
decided model_loop); ``simprep variants`` protonates every system it wrote when the
manifest states a pH. Both write ``protonation_record.json`` and ``protonation_report.md``.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from importlib.metadata import version as installed_version
from pathlib import Path

from simprep.findings import residue_ref
from simprep.manifest import config_from_manifest, write_json
from simprep.model.run import MODELLED_DIR, model_wild_type, single_system, write_model
from simprep.prep.plan import build_plan
from simprep.prep.record import system_counts
from simprep.prep.run import PrepInputs, PrepRequest, load_inputs, write_prep
from simprep.prep.write import write_system
from simprep.protonate.hydrogens import add_hydrogens
from simprep.protonate.pka import estimate
from simprep.protonate.report import render_protonation_report
from simprep.protonate.states import VARIANTS, Context, ProtonationError, State, assign
from simprep.provenance import simprep_provenance, utc_now
from simprep.schemas import validate
from simprep.structure.model import ResidueClass, Structure

RECORD_FILE = "protonation_record.json"
REPORT_FILE = "protonation_report.md"
SYSTEM_NAME = "system"
WILD_TYPE = "wt"
PROTONATED_SUFFIX = "_protonated"
PROTONATION_STAGE = "protonation"


@dataclass(frozen=True)
class Protonated:
    """One protonated system: where it came from and went, its states and findings."""

    name: str
    source: str
    directory: str
    structure: Structure
    states: list[State]
    findings: list[dict]


def settings(inputs: PrepInputs) -> tuple[dict, dict]:
    """(the manifest's protonation, the method) or ProtonationError if the pH is missing."""
    chosen = inputs.manifest.get("protonation")
    if chosen is None:
        raise ProtonationError(
            "the manifest states no pH: add `protonation` with `ph` and `ph_rationale` "
            "(TASK-009 decision 1; simprep never defaults the pH)"
        )
    return chosen, chosen.get("method") or inputs.ruleset.protonation["method"]


def protonate(name: str, structure: Structure, where: tuple, inputs: PrepInputs) -> Protonated:
    """Protonate ``structure`` (system ``name``); ``where`` is (source, target directory)."""
    chosen, method = settings(inputs)
    estimates = estimate(structure, method)
    _check_covered(structure, {e.residue for e in estimates})
    decisions = {d["finding_id"]: d for d in inputs.manifest["decisions"]}
    config = config_from_manifest(inputs.manifest, inputs.manifest_sha256)
    context = Context(chosen["ph"], method, inputs.ruleset, decisions, name, config)
    states, findings = assign(structure, estimates, context)
    protonation = {**inputs.ruleset.protonation, "method": method}
    protonated = add_hydrogens(structure, states, protonation)
    source, directory = where
    return Protonated(
        name,
        source,
        directory,
        protonated,
        [_tautomer(s, protonated) for s in states],
        [f.to_dict() for f in findings],
    )


def _tautomer(state: State, structure: Structure) -> State:
    """Neutral His: record the tautomer OpenMM chose (from the hydrogens it added)."""
    if state.variant is not None:
        return state
    names = {a.name for a in structure.residue(state.residue).atoms}
    variant = "HID" if "HD1" in names else "HIE"
    return replace(state, variant=variant, basis=f"{state.basis}; tautomer by hydrogen bonding")


def _check_covered(structure: Structure, estimated: set) -> None:
    """Every titratable residue needs an estimate; a missing one would get OpenMM's default."""
    missing = [
        f"{r.name} {r.id.label()}"
        for r in structure.residues
        if r.residue_class is ResidueClass.POLYMER and r.name in VARIANTS and r.id not in estimated
    ]
    if missing:
        raise ProtonationError(f"PROPKA gave no pKa for {', '.join(missing)}")


def run_protonate(request: PrepRequest) -> dict:
    """``simprep protonate``: prepare (and model) the wild type, protonate it, write it."""
    inputs = load_inputs(request)
    settings(inputs)  # fail before any work when the pH is missing
    plan = build_plan(inputs.structure, inputs.manifest, inputs.ruleset)
    prepared = single_system(inputs, plan)
    modelled = model_wild_type(inputs, plan, prepared)
    wt_record = write_prep(inputs, plan, request.out_dir / WILD_TYPE)
    source, structure, work_order = WILD_TYPE, prepared, wt_record["work_order"]
    if modelled is not None:
        model_record = write_model(modelled, (inputs, wt_record), request.out_dir)
        source, structure, work_order = MODELLED_DIR, modelled.structure, model_record["work_order"]
    wild_type = protonate(WILD_TYPE, structure, (source, WILD_TYPE + PROTONATED_SUFFIX), inputs)
    return write_protonation([wild_type], (inputs, work_order), request.out_dir)


def write_protonation(systems: list[Protonated], context: tuple, out_dir: Path) -> dict:
    """Write every system and the record; ``systems[0]`` is the wild type the others are
    compared with; ``context`` is (inputs, the wild type's work order)."""
    inputs, work_order = context
    chosen, method = settings(inputs)
    entries = [_entry(s, (systems[0], inputs.structure), out_dir) for s in systems]
    record = {
        "schema_version": "0.1.0",
        "generated_at": utc_now(),
        "simprep": simprep_provenance(),
        "knowledge_base": {"version": inputs.ruleset.version, "sha256": inputs.ruleset.sha256},
        "input": inputs.source,
        "manifest_sha256": inputs.manifest_sha256,
        "ph": chosen["ph"],
        "ph_rationale": chosen["ph_rationale"],
        "method": method,
        "propka_version": installed_version("propka"),
        "systems": entries,
        "work_order": [w for w in work_order if w["stage"] != PROTONATION_STAGE],
    }
    validate(record, "protonation_record")
    write_json(record, out_dir / RECORD_FILE)
    (out_dir / REPORT_FILE).write_text(render_protonation_report(record))
    return record


def _entry(system: Protonated, reference: tuple, out_dir: Path) -> dict:
    wild_type, deposited = reference
    directory = out_dir / system.directory
    directory.mkdir(parents=True, exist_ok=True)
    return {
        "name": system.name,
        "directory": system.directory,
        "source": system.source,
        "files": write_system(system.structure, directory, SYSTEM_NAME),
        "counts": system_counts(deposited, system.structure),
        "states": [_state(s) for s in system.states],
        "findings": system.findings,
        "differs_from_wild_type": _differences(system, wild_type),
    }


def _state(state: State) -> dict:
    return {
        **residue_ref(state.residue, state.res_name),
        "variant": state.variant,
        "pka": None if state.pka is None else round(state.pka, 2),
        "basis": state.basis,
    }


def _differences(system: Protonated, wild_type: Protonated) -> list[dict]:
    reference = {s.residue: s for s in wild_type.states}
    return [
        {
            **residue_ref(s.residue, s.res_name),
            "wild_type_variant": reference[s.residue].variant,
            "variant": s.variant,
        }
        for s in system.states
        if s.residue in reference and reference[s.residue].variant != s.variant
    ]
