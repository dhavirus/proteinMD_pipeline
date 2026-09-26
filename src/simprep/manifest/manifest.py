"""Manifest v0.1: every decision about a system, with provenance."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from simprep.canonical import sha256_canonical
from simprep.config import AuditConfig, region_from_dict, thresholds_from_dict
from simprep.rules import RuleSet
from simprep.schemas import SCHEMA_VERSION, validate


class ManifestError(ValueError):
    """The manifest is inconsistent with its input file or its own findings."""


def pretty_json(document: object) -> str:
    """Readable, stable file layout (hashes use :mod:`simprep.canonical`, not this)."""
    return json.dumps(document, indent=2, sort_keys=True) + "\n"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def write_json(document: object, path: Path) -> None:
    path.write_text(pretty_json(document))


def load_manifest(path: Path) -> dict:
    """Read and schema-validate a manifest; snapshot hash and decisions are checked too."""
    manifest = json.loads(path.read_text())
    validate(manifest, "manifest")
    check_snapshot(manifest)
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


def check_snapshot(manifest: dict) -> None:
    """Each snapshot's findings must hash to its recorded ``findings_sha256`` (RFC 8785)."""
    remedies = {
        "findings_snapshot": "simprep audit --manifest",
        "variant_snapshot": "simprep variants",
    }
    for key, command in remedies.items():
        snapshot = manifest.get(key)
        if snapshot is None:
            continue
        actual = sha256_canonical(snapshot["findings"])
        if actual != snapshot["findings_sha256"]:
            raise ManifestError(
                f"{key} was modified: recorded SHA-256 {snapshot['findings_sha256']}, "
                f"actual {actual}. Re-run `{command}` to refresh it."
            )


def check_decisions(manifest: dict) -> None:
    """Every decision must name a snapshot finding, once, with one of its options.

    All problems are reported together, so orphaned decisions after a re-audit are
    listed in full rather than one at a time.
    """
    snapshot = manifest["findings_snapshot"]
    if not manifest["decisions"]:
        return
    if snapshot is None:
        raise ManifestError("manifest has decisions but no findings snapshot to decide on")
    variant_snapshot = manifest.get("variant_snapshot") or {"findings": []}
    findings = snapshot["findings"] + variant_snapshot["findings"]
    options = {f["id"]: {o["id"] for o in f["options"]} for f in findings}
    problems = decision_problems(manifest["decisions"], options)
    if problems:
        raise ManifestError(
            f"{len(problems)} decision problem(s):\n  "
            + "\n  ".join(problems)
            + "\nRemove or re-decide these findings in the front end, then export again."
        )


def decision_problems(decisions: list[dict], options: dict[str, set[str]]) -> list[str]:
    """Human-readable problems: orphaned findings, unknown options, duplicate decisions."""
    problems, seen = [], set()
    for decision in decisions:
        finding_id, option_id = decision["finding_id"], decision["option_id"]
        who = f"(option {option_id!r}, decided by {decision['decided_by']})"
        if finding_id in seen:
            problems.append(f"{finding_id}: more than one decision {who}")
        seen.add(finding_id)
        if finding_id not in options:
            problems.append(f"{finding_id}: finding not in the snapshot (orphaned) {who}")
        elif option_id not in options[finding_id]:
            valid = ", ".join(sorted(options[finding_id]))
            problems.append(f"{finding_id}: unknown option {who}; valid: {valid}")
    return problems


def attach_snapshot(manifest: dict, report: dict) -> dict:
    """Copy of ``manifest`` carrying the findings of ``report`` as its snapshot."""
    findings = report["findings"]
    updated = {
        **manifest,
        "findings_snapshot": {
            "generated_at": report["generated_at"],
            "findings_sha256": sha256_canonical(findings),
            "findings": findings,
        },
    }
    validate(updated, "manifest")
    check_decisions(updated)
    return updated


def attach_variant_snapshot(manifest: dict, findings: list[dict], generated_at: str) -> dict:
    """Copy of ``manifest`` carrying variant_build ``findings`` for its current variants."""
    updated = {
        **manifest,
        "variant_snapshot": {
            "generated_at": generated_at,
            "variants_sha256": sha256_canonical(manifest.get("variants", [])),
            "findings_sha256": sha256_canonical(findings),
            "findings": findings,
        },
    }
    validate(updated, "manifest")
    check_decisions(updated)
    return updated
