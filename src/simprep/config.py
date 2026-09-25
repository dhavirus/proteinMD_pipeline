"""Audit configuration: regions of interest and severity-escalation thresholds."""

from __future__ import annotations

from dataclasses import dataclass

from simprep.structure.model import ResidueId


@dataclass(frozen=True)
class Region:
    """A named residue set from the manifest."""

    name: str
    description: str
    residues: frozenset[ResidueId]

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "description": self.description,
            "residues": [rid.to_dict() for rid in sorted(self.residues)],
        }


@dataclass(frozen=True)
class Threshold:
    max_distance_angstrom: float
    raise_levels: int


@dataclass(frozen=True)
class AuditConfig:
    """Everything context-dependent an audit uses; recorded verbatim in findings.json."""

    regions: tuple[Region, ...]
    thresholds: tuple[Threshold, ...]
    manifest_sha256: str | None = None

    def to_dict(self) -> dict:
        return {
            "manifest_sha256": self.manifest_sha256,
            "regions": [region.to_dict() for region in self.regions],
            "severity_escalation": {
                "thresholds": [
                    {
                        "max_distance_angstrom": t.max_distance_angstrom,
                        "raise_levels": t.raise_levels,
                    }
                    for t in self.thresholds
                ]
            },
        }


def region_from_dict(data: dict) -> Region:
    residues = frozenset(
        ResidueId(ref["chain"], ref["seq_num"], ref.get("ins_code", "")) for ref in data["residues"]
    )
    return Region(data["name"], data["description"], residues)


def thresholds_from_dict(severity_escalation: dict) -> tuple[Threshold, ...]:
    return tuple(
        Threshold(float(t["max_distance_angstrom"]), int(t["raise_levels"]))
        for t in severity_escalation["thresholds"]
    )


def default_config(audit_defaults: dict) -> AuditConfig:
    """No regions; thresholds from the knowledge base's audit defaults."""
    return AuditConfig((), thresholds_from_dict(audit_defaults["severity_escalation"]))
