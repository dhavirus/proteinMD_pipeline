"""Knowledge-base rule model (pure; file loading lives in :mod:`simprep.knowledge`)."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Rule:
    """One validated knowledge-base rule; ``matcher`` is interpreted by its family's detector."""

    rule_id: str
    family: str
    version: str
    title: str
    priority: int
    matcher: dict
    base_severity: str
    options: tuple[dict, ...]
    recommended_option: dict
    rationale: str
    references: tuple[dict, ...]
    verify_flags: tuple[str, ...]

    @classmethod
    def from_dict(cls, data: dict) -> Rule:
        return cls(
            **{
                **data,
                "options": tuple(data["options"]),
                "references": tuple(data["references"]),
                "verify_flags": tuple(data["verify_flags"]),
            }
        )

    def option_ids(self) -> tuple[str, ...]:
        return tuple(option["id"] for option in self.options)


@dataclass(frozen=True)
class RuleSet:
    """All rules of one knowledge-base version, with a content hash for provenance."""

    version: str
    sha256: str
    rules: tuple[Rule, ...]
    audit_defaults: dict
    residue_mappings: tuple[dict, ...] = ()
    side_chains: dict = field(default_factory=dict)
    relaxation: dict = field(default_factory=dict)
    modelling: dict = field(default_factory=dict)
    protonation: dict = field(default_factory=dict)

    def family(self, family: str) -> tuple[Rule, ...]:
        """Rules of ``family``, highest priority first, then by rule id (deterministic)."""
        members = [rule for rule in self.rules if rule.family == family]
        return tuple(sorted(members, key=lambda rule: (-rule.priority, rule.rule_id)))


class RuleConsistencyError(ValueError):
    """A rule is schema-valid but internally inconsistent (e.g. unknown option id)."""


def check_unique_ids(rules: tuple[Rule, ...]) -> None:
    counts = Counter(rule.rule_id for rule in rules)
    duplicates = sorted(rule_id for rule_id, count in counts.items() if count > 1)
    if duplicates:
        raise RuleConsistencyError(f"duplicate rule ids: {duplicates}")


def check_rule(rule: Rule, file_family: str) -> None:
    """Raise if the rule's family or option references are inconsistent."""
    if rule.family != file_family:
        raise RuleConsistencyError(f"{rule.rule_id}: family {rule.family} in a {file_family} file")
    known = set(rule.option_ids())
    recommendation = rule.recommended_option
    referenced = [recommendation["default"]] + [c["option"] for c in recommendation["conditions"]]
    unknown = sorted(set(referenced) - known)
    if unknown:
        raise RuleConsistencyError(f"{rule.rule_id}: recommends unknown option(s) {unknown}")
    explicit_only = {o["id"] for o in rule.options if o.get("requires_explicit_choice")}
    if explicit_only & set(referenced):
        raise RuleConsistencyError(
            f"{rule.rule_id}: recommends {sorted(explicit_only & set(referenced))}, "
            "which require an explicit human choice"
        )
