"""Build a decided manifest for a structure (tests only): audit, snapshot, decisions."""

from __future__ import annotations

from pathlib import Path

from simprep.audit import ReportContext, build_findings_report, run_audit
from simprep.config import default_config
from simprep.knowledge import load_ruleset
from simprep.manifest import attach_snapshot, init_manifest, write_json
from simprep.provenance import input_info
from simprep.structure.parse import read_structure

FIXED_TIME = "2026-09-26T00:00:00Z"
TEST_PROVENANCE = {"version": "test", "git_commit": None}


def audited_manifest(structure_path: Path) -> dict:
    """A manifest for ``structure_path`` with its findings snapshot and no decisions."""
    ruleset = load_ruleset()
    source = input_info(structure_path)
    manifest = init_manifest(source, ruleset, (TEST_PROVENANCE, FIXED_TIME))
    structure = read_structure(structure_path)
    config = default_config(ruleset.audit_defaults)
    findings = run_audit(structure, ruleset, config)
    context = ReportContext(source, ruleset, config, TEST_PROVENANCE, FIXED_TIME)
    return attach_snapshot(manifest, build_findings_report(structure, findings, context))


def decision(finding_id: str, option_id: str, rationale: str = "test", parameters=None) -> dict:
    entry = {
        "finding_id": finding_id,
        "option_id": option_id,
        "rationale": rationale,
        "decided_by": "test",
        "timestamp": FIXED_TIME,
    }
    return entry | ({"parameters": parameters} if parameters else {})


def decided_manifest(structure_path: Path, decisions: dict, out_path: Path) -> Path:
    """Write the audited manifest with ``decisions`` ({finding_id: {option, rationale,
    parameters?}}) to ``out_path``."""
    manifest = audited_manifest(structure_path)
    manifest["decisions"] = [
        decision(fid, d["option"], d["rationale"], d.get("parameters"))
        for fid, d in decisions.items()
    ]
    write_json(manifest, out_path)
    return out_path
