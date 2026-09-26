"""Manifest v0.1: every decision about a system, with provenance."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from simprep.config import AuditConfig, region_from_dict, thresholds_from_dict
from simprep.rules import RuleSet
from simprep.schemas import SCHEMA_VERSION, validate


class ManifestError(ValueError):
    """The manifest is inconsistent with its input file or its own findings."""


def canonical_json(document: object) -> str:
    return json.dumps(document, indent=2, sort_keys=True) + "\n"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def write_json(document: object, path: Path) -> None:
    path.write_text(canonical_json(document))


def load_manifest(path: Path) -> dict:
    """Read and schema-validate a manifest; decisions are checked against its snapshot."""
    manifest = json.loads(path.read_text())
    validate(manifest, "manifest")
    check_decisions(manifest)
    return manifest


def init_manifest(input_info: dict, ruleset: RuleSet, provenance: tuple[dict, str]) -> dict:
    """A new manifest for ``input_info``: default thresholds, no regions, no decisions."""
    simprep, created_at = provenance
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at": created_at,
        "input": input_info,
        "knowledge_base": {"version": ruleset.version, "sha256": ruleset.sha256},
        "simprep": simprep,
        "regions": [],
        "severity_escalation": ruleset.audit_defaults["severity_escalation"],
        "findings_snapshot": None,
        "decisions": [],
    }
    validate(manifest, "manifest")
    return manifest


def check_input(manifest: dict, input_sha256: str) -> None:
    """Refuse to apply a manifest built for a different input file."""
    expected = manifest["input"]["sha256"]
    if expected != input_sha256:
        raise ManifestError(
            f"manifest was built for input SHA-256 {expected}, but the file given has "
            f"{input_sha256}. Use the manifest's own input file or run `simprep manifest init`."
        )


def config_from_manifest(manifest: dict, manifest_sha256: str) -> AuditConfig:
    return AuditConfig(
        regions=tuple(region_from_dict(region) for region in manifest["regions"]),
        thresholds=thresholds_from_dict(manifest["severity_escalation"]),
        manifest_sha256=manifest_sha256,
    )


def check_decisions(manifest: dict) -> None:
    """Every decision must name a snapshot finding and one of that finding's options."""
    snapshot = manifest["findings_snapshot"]
    if not manifest["decisions"]:
        return
    if snapshot is None:
        raise ManifestError("manifest has decisions but no findings snapshot to decide on")
    options = {f["id"]: {o["id"] for o in f["options"]} for f in snapshot["findings"]}
    for decision in manifest["decisions"]:
        finding_id = decision["finding_id"]
        if finding_id not in options:
            raise ManifestError(f"decision refers to unknown finding {finding_id!r}")
        if decision["option_id"] not in options[finding_id]:
            raise ManifestError(
                f"decision for {finding_id!r} picks unknown option {decision['option_id']!r}; "
                f"valid: {sorted(options[finding_id])}"
            )


def attach_snapshot(manifest: dict, report: dict) -> dict:
    """Copy of ``manifest`` carrying the findings of ``report`` as its snapshot."""
    findings = report["findings"]
    updated = {
        **manifest,
        "findings_snapshot": {
            "generated_at": report["generated_at"],
            "findings_sha256": sha256_text(canonical_json(findings)),
            "findings": findings,
        },
    }
    validate(updated, "manifest")
    check_decisions(updated)
    return updated
