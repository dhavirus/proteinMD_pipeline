"""Prep orchestration at the I/O edge: read, gate, plan, apply, write, record."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from simprep.knowledge import load_ruleset
from simprep.manifest import check_input, load_manifest, write_json
from simprep.manifest.manifest import sha256_text
from simprep.paths import KNOWLEDGE_DIR
from simprep.prep.apply import apply_plan
from simprep.prep.plan import PrepPlan, build_plan
from simprep.prep.record import RecordContext, build_prep_record, render_prep_report, system_entry
from simprep.prep.write import write_system
from simprep.provenance import input_info, simprep_provenance, utc_now
from simprep.rules import RuleSet
from simprep.schemas import validate
from simprep.structure.model import Structure
from simprep.structure.parse import read_structure

RECORD_FILE = "prep_record.json"
REPORT_FILE = "prep_report.md"


@dataclass(frozen=True)
class PrepRequest:
    """Inputs of one prep run."""

    structure: Path
    manifest: Path
    out_dir: Path
    knowledge: Path = KNOWLEDGE_DIR


@dataclass(frozen=True)
class PrepInputs:
    """Everything a run reads, checked: the input file matches the manifest."""

    source: dict
    manifest: dict
    manifest_sha256: str
    ruleset: RuleSet
    structure: Structure


def load_inputs(request: PrepRequest) -> PrepInputs:
    source = input_info(request.structure)
    manifest = load_manifest(request.manifest)
    check_input(manifest, source["sha256"])
    return PrepInputs(
        source,
        manifest,
        sha256_text(request.manifest.read_text()),
        load_ruleset(request.knowledge),
        read_structure(request.structure),
    )


def run_prep(request: PrepRequest) -> dict:
    """Prepare the system(s) for ``request``; write coordinates, prep_record.json and
    prep_report.md into ``request.out_dir`` and return the (schema-valid) record."""
    inputs = load_inputs(request)
    plan = build_plan(inputs.structure, inputs.manifest, inputs.ruleset)
    return write_prep(inputs, plan, request.out_dir)


def write_prep(inputs: PrepInputs, plan: PrepPlan, out_dir: Path) -> dict:
    """Apply ``plan`` and write the system files, prep_record.json and prep_report.md."""
    out_dir.mkdir(parents=True, exist_ok=True)
    systems = [
        system_entry(inputs.structure, system, write_system(system.structure, out_dir, system.name))
        for system in apply_plan(inputs.structure, plan)
    ]
    context = RecordContext(
        inputs.source, inputs.manifest_sha256, inputs.ruleset, simprep_provenance(), utc_now()
    )
    record = build_prep_record(plan, systems, context)
    validate(record, "prep_record")
    write_json(record, out_dir / RECORD_FILE)
    (out_dir / REPORT_FILE).write_text(render_prep_report(record))
    return record
