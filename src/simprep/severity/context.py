"""Escalate severity by distance to manifest regions; resolve recommended options.

A finding's base severity comes from its rule. Its effective severity is raised by the
tightest configured threshold that its anchor residues fall within (heavy atom to heavy
atom; an anchor residue inside a region counts as distance 0). With no regions defined,
effective severity equals base severity.
"""

from __future__ import annotations

from simprep.config import AuditConfig, Region, Threshold
from simprep.findings import Finding, SeverityContext
from simprep.structure.geometry import min_distance
from simprep.structure.model import ResidueId, Structure


def apply_context(
    findings: list[Finding], structure: Structure, config: AuditConfig
) -> list[Finding]:
    """Findings with effective severity, severity context and recommendation filled in."""
    resolved = []
    for finding in findings:
        context = severity_context(finding, structure, config)
        resolved.append(finding.with_context(context, recommend(finding, context)))
    return resolved


def severity_context(
    finding: Finding, structure: Structure, config: AuditConfig
) -> SeverityContext:
    if not config.regions:
        return SeverityContext()
    nearest = min(
        (
            (_region_distance(finding.anchor_residues, region, structure), region.name)
            for region in config.regions
        ),
        key=lambda pair: (pair[0] is None, pair[0] if pair[0] is not None else 0.0, pair[1]),
    )
    distance_angstrom, region_name = nearest
    if distance_angstrom is None:
        return SeverityContext("far", None, None, 0)
    threshold = _tightest(config.thresholds, distance_angstrom)
    if distance_angstrom == 0.0:
        proximity = "inside"
    else:
        proximity = "near" if threshold is not None else "far"
    raise_levels = threshold.raise_levels if threshold is not None else 0
    return SeverityContext(proximity, region_name, distance_angstrom, raise_levels)


def _region_distance(
    anchors: tuple[ResidueId, ...], region: Region, structure: Structure
) -> float | None:
    """0 if an anchor is a region residue, else the minimum heavy-atom distance."""
    if any(anchor in region.residues for anchor in anchors):
        return 0.0
    return min_distance(
        _heavy_positions(anchors, structure),
        _heavy_positions(tuple(sorted(region.residues)), structure),
    )


def _heavy_positions(residue_ids: tuple[ResidueId, ...], structure: Structure) -> list[tuple]:
    index = structure.residue_index
    return [atom.position for rid in residue_ids if rid in index for atom in index[rid].heavy_atoms]


def _tightest(thresholds: tuple[Threshold, ...], distance_angstrom: float) -> Threshold | None:
    matching = [t for t in thresholds if distance_angstrom <= t.max_distance_angstrom]
    return min(matching, key=lambda t: t.max_distance_angstrom) if matching else None


def recommend(finding: Finding, context: SeverityContext) -> tuple[str, str]:
    """(option id, basis): the first rule condition matching the context and evidence."""
    recommendation = finding.rule.recommended_option
    for condition in recommendation["conditions"]:
        if _matches(condition["when"], finding, context):
            return condition["option"], f"condition {condition['name']}: {condition['reason']}"
    return recommendation["default"], "default"


def _matches(when: dict, finding: Finding, context: SeverityContext) -> bool:
    if "region_proximity" in when and context.proximity not in when["region_proximity"]:
        return False
    return all(
        _value_in(finding.evidence_value(key), allowed)
        for key, allowed in when.get("evidence", {}).items()
    )


def _value_in(value: object, allowed: list) -> bool:
    """Membership that does not conflate booleans with integers."""
    return any(type(value) is type(option) and value == option for option in allowed)
