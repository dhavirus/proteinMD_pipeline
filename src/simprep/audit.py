"""Audit orchestration: detectors -> unrecognized pass -> severity context -> report.

``run_audit`` and ``build_findings_report`` are pure; file access lives in the CLI.
"""

from __future__ import annotations

from dataclasses import dataclass

from simprep.config import AuditConfig
from simprep.detectors import load_detectors
from simprep.detectors.unrecognized import detect_unrecognized
from simprep.findings import Finding
from simprep.rules import RuleSet
from simprep.schemas import SCHEMA_VERSION
from simprep.severity import apply_context
from simprep.structure.model import ResidueClass, ResidueId, Structure

FAMILY_ORDER = ("metals", "nonstandard_residues", "altlocs", "missing_residues", "unrecognized")
SCOPE = (
    "Model 1 of the deposited asymmetric unit as parsed; no biological assembly is "
    "generated and no coordinates are modified by an audit."
)
ANNOTATION_CATEGORIES = (
    "struct_conn",
    "pdbx_struct_mod_residue",
    "pdbx_unobs_or_zero_occ_residues",
)


@dataclass(frozen=True)
class ReportContext:
    """Provenance and configuration stamped onto findings.json."""

    input: dict
    ruleset: RuleSet
    config: AuditConfig
    simprep: dict
    generated_at: str


def run_audit(structure: Structure, ruleset: RuleSet, config: AuditConfig) -> list[Finding]:
    """All findings for ``structure``, context-weighted and in deterministic order."""
    findings = [
        finding
        for detector in load_detectors().values()
        for finding in detector(structure, ruleset, config)
    ]
    claimed = frozenset(rid for finding in findings for rid in finding.claims)
    findings += detect_unrecognized(structure, ruleset, claimed)
    return sorted(apply_context(findings, structure, config), key=_sort_key)


def _sort_key(finding: Finding) -> tuple:
    residue = finding.locus.residue
    position = (residue.chain, residue.seq_num, residue.ins_code) if residue else ("", 0, "")
    return (FAMILY_ORDER.index(finding.rule.family), position, finding.id)


def record_counts(structure: Structure, findings: list[Finding]) -> list[dict]:
    """In / passed-through / flagged / excluded counts per record type.

    An audit applies no decisions, so ``excluded`` is 0 and every record is either
    passed through or flagged (covered by a finding awaiting a decision).
    """
    flagged = frozenset(
        rid for f in findings for rid in (*f.claims, *(r for r, _ in f.locus.extent))
    )
    referenced_links = {
        item["fields"]["id"]
        for f in findings
        for item in f.evidence
        if item["type"] == "source_record" and item["category"] == "struct_conn"
    }
    rows = [_residue_row(structure, flagged, cls) for cls in ResidueClass]
    rows.append(
        _row("atom records", [r.id in flagged for r in structure.residues for _ in r.atoms])
    )
    rows.append(
        _row(
            "unobserved polymer residues (annotation)",
            [u.residue in flagged for u in structure.unobserved_residues if u.is_polymer],
        )
    )
    rows.append(
        _row("struct_conn records", [link.conn_id in referenced_links for link in structure.links])
    )
    return rows


def _residue_row(structure: Structure, flagged: frozenset[ResidueId], cls: ResidueClass) -> dict:
    members = [r.id in flagged for r in structure.residues if r.residue_class is cls]
    return _row(f"{cls.value} residues", members)


def _row(record_type: str, flags: list[bool]) -> dict:
    return {
        "record_type": record_type,
        "in": len(flags),
        "passed_through": flags.count(False),
        "flagged": flags.count(True),
        "excluded": 0,
    }


def build_findings_report(
    structure: Structure, findings: list[Finding], context: ReportContext
) -> dict:
    """The findings.json document (validate with the findings_report schema)."""
    present = sorted(structure.annotation_categories)
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": context.generated_at,
        "simprep": context.simprep,
        "knowledge_base": {"version": context.ruleset.version, "sha256": context.ruleset.sha256},
        "input": context.input,
        "audit_config": context.config.to_dict(),
        "structure_summary": {
            "name": structure.name,
            "scope": SCOPE,
            "annotation_categories_present": present,
            "annotation_categories_absent": sorted(set(ANNOTATION_CATEGORIES) - set(present)),
        },
        "counts": record_counts(structure, findings),
        "findings": [finding.to_dict() for finding in findings],
    }
