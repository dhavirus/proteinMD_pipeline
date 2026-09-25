"""Finding model and typed evidence constructors.

Detectors build :class:`Finding` objects from a rule plus locus and evidence; the
severity step later fills the context-dependent fields. Serialization follows
``schema/finding.schema.json``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from simprep.rules import Rule
from simprep.structure.model import Atom, LinkPartner, Residue, ResidueId

DISTANCE_DECIMALS = 3
VALUE_DECIMALS = 2
SEVERITY_LEVELS = ("info", "warn", "blocking")


@dataclass(frozen=True)
class Locus:
    """Where a finding is: a primary residue (or none) and every residue it covers."""

    residue: ResidueId | None
    res_name: str | None
    atom_names: tuple[str, ...] = ()
    extent: tuple[tuple[ResidueId, str], ...] = ()

    def to_dict(self) -> dict:
        residue = self.residue
        return {
            "chain": residue.chain if residue else None,
            "seq_num": residue.seq_num if residue else None,
            "ins_code": residue.ins_code if residue else "",
            "res_name": self.res_name,
            "atom_names": list(self.atom_names),
            "extent": [residue_ref(rid, name) for rid, name in self.extent],
        }


@dataclass(frozen=True)
class SeverityContext:
    proximity: str = "no_regions"
    nearest_region: str | None = None
    min_distance_angstrom: float | None = None
    raise_levels: int = 0

    def to_dict(self) -> dict:
        return {
            "proximity": self.proximity,
            "nearest_region": self.nearest_region,
            "min_distance_angstrom": _round(self.min_distance_angstrom, DISTANCE_DECIMALS),
            "raise_levels": self.raise_levels,
        }


@dataclass(frozen=True)
class Finding:
    """One detected feature. ``claims`` lists residues whose chemistry it accounts for."""

    id: str
    rule: Rule
    title: str
    locus: Locus
    anchor_residues: tuple[ResidueId, ...]
    evidence: tuple[dict, ...]
    claims: frozenset[ResidueId] = frozenset()
    context: SeverityContext = field(default_factory=SeverityContext)
    effective_severity: str | None = None
    recommended_option: str | None = None
    recommendation_basis: str | None = None

    def with_context(self, context: SeverityContext, recommendation: tuple[str, str]) -> Finding:
        """Copy with context-dependent severity and recommendation filled in."""
        return replace(
            self,
            context=context,
            effective_severity=raise_severity(self.rule.base_severity, context.raise_levels),
            recommended_option=recommendation[0],
            recommendation_basis=recommendation[1],
        )

    def evidence_value(self, key: str) -> object | None:
        """Value of the first evidence item called ``key`` (None if absent)."""
        return next((item.get("value") for item in self.evidence if item["key"] == key), None)

    def to_dict(self) -> dict:
        if self.effective_severity is None:
            raise ValueError(f"{self.id}: severity context not applied before serialization")
        return {
            "id": self.id,
            "rule_id": self.rule.rule_id,
            "rule_family": self.rule.family,
            "title": self.title,
            "base_severity": self.rule.base_severity,
            "effective_severity": self.effective_severity,
            "severity_context": self.context.to_dict(),
            "locus": self.locus.to_dict(),
            "anchor_residues": [rid.to_dict() for rid in self.anchor_residues],
            "evidence": list(self.evidence),
            "options": list(self.rule.options),
            "recommended_option": self.recommended_option,
            "recommendation_basis": self.recommendation_basis,
            "rationale": self.rule.rationale,
            "references": list(self.rule.references),
            "verify_flags": list(self.rule.verify_flags),
        }


def raise_severity(severity: str, levels: int) -> str:
    """Raise ``severity`` by ``levels`` steps, saturating at ``blocking``."""
    index = min(SEVERITY_LEVELS.index(severity) + levels, len(SEVERITY_LEVELS) - 1)
    return SEVERITY_LEVELS[index]


def residue_locus(residue: Residue, atom_names: tuple[str, ...] = ()) -> Locus:
    return Locus(residue.id, residue.name, atom_names, ((residue.id, residue.name),))


def residue_ref(residue_id: ResidueId, res_name: str | None = None) -> dict:
    ref = residue_id.to_dict()
    if res_name:
        ref["res_name"] = res_name
    return ref


def atom_ref(residue: Residue, atom: Atom) -> dict:
    return {
        **residue.id.to_dict(),
        "res_name": residue.name,
        "atom_name": atom.name,
        "altloc": atom.altloc,
    }


def partner_ref(partner: LinkPartner) -> dict:
    return {
        **partner.residue.to_dict(),
        "res_name": partner.res_name,
        "atom_name": partner.atom_name,
        "altloc": partner.altloc,
    }


def _round(value: float | None, decimals: int) -> float | None:
    return None if value is None else round(value, decimals)


def _item(key: str, kind: str, note: str | None, **fields: object) -> dict:
    item = {"key": key, "type": kind, **fields}
    if note:
        item["note"] = note
    return item


def ev_distance(key: str, value: float, atoms: list[dict], note: str | None = None) -> dict:
    return _item(
        key, "distance", note, value=round(value, DISTANCE_DECIMALS), unit="angstrom", atoms=atoms
    )


def ev_angle(key: str, value: float, atoms: list[dict]) -> dict:
    return _item(key, "angle", None, value=round(value, VALUE_DECIMALS), unit="degree", atoms=atoms)


def ev_occupancy(key: str, value: float, note: str | None = None) -> dict:
    return _item(key, "occupancy", note, value=round(value, VALUE_DECIMALS))


def ev_b_factor(key: str, value: float, note: str | None = None) -> dict:
    return _item(key, "b_factor", note, value=round(value, VALUE_DECIMALS), unit="angstrom_squared")


def ev_count(key: str, value: int, note: str | None = None) -> dict:
    return _item(key, "count", note, value=value)


def ev_number(key: str, value: float, unit: str, note: str | None = None) -> dict:
    return _item(key, "number", note, value=round(value, DISTANCE_DECIMALS), unit=unit)


def ev_label(key: str, value: str, note: str | None = None) -> dict:
    return _item(key, "label", note, value=value)


def ev_source_record(key: str, category: str, fields: dict, note: str | None = None) -> dict:
    return _item(key, "source_record", note, category=category, fields=fields)


def ev_flag(key: str, value: bool, note: str | None = None) -> dict:
    return _item(key, "flag", note, value=value)
