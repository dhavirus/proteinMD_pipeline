"""Decision coverage of a manifest: which findings are decided, which still block prep."""

from __future__ import annotations

from dataclasses import dataclass

from simprep.manifest.manifest import ManifestError

SEVERITIES = ("blocking", "warn", "info")


@dataclass(frozen=True)
class DecisionStatus:
    """Counts per effective severity plus the findings that still need a decision."""

    decided: dict[str, int]
    undecided: dict[str, int]
    undecided_findings: tuple[dict, ...]
    explicit_choice_decisions: tuple[dict, ...]

    @property
    def blocking_undecided(self) -> bool:
        return self.undecided["blocking"] > 0


def decision_status(manifest: dict) -> DecisionStatus:
    """Summarize decision coverage of a validated manifest (needs a findings snapshot)."""
    snapshot = manifest["findings_snapshot"]
    if snapshot is None:
        raise ManifestError(
            "manifest has no findings snapshot; run `simprep audit FILE --manifest M` first"
        )
    findings = snapshot["findings"]
    decided_ids = {decision["finding_id"] for decision in manifest["decisions"]}
    undecided = tuple(f for f in findings if f["id"] not in decided_ids)
    return DecisionStatus(
        decided=_count(f for f in findings if f["id"] in decided_ids),
        undecided=_count(undecided),
        undecided_findings=tuple(sorted(undecided, key=_severity_order)),
        explicit_choice_decisions=_explicit_choices(findings, manifest["decisions"]),
    )


def _count(findings) -> dict[str, int]:
    counts = dict.fromkeys(SEVERITIES, 0)
    for finding in findings:
        counts[finding["effective_severity"]] += 1
    return counts


def _severity_order(finding: dict) -> tuple[int, str]:
    return SEVERITIES.index(finding["effective_severity"]), finding["id"]


def _explicit_choices(findings: list[dict], decisions: list[dict]) -> tuple[dict, ...]:
    """Decisions whose option is marked requires_explicit_choice, kept visible at review."""
    explicit = {
        (f["id"], o["id"])
        for f in findings
        for o in f["options"]
        if o.get("requires_explicit_choice")
    }
    return tuple(d for d in decisions if (d["finding_id"], d["option_id"]) in explicit)


def render_status(status: DecisionStatus) -> str:
    """Plain-text report for the terminal."""
    lines = ["severity   decided  undecided"]
    lines += [f"{s:<10} {status.decided[s]:>7}  {status.undecided[s]:>9}" for s in SEVERITIES]
    if status.undecided_findings:
        lines.append("")
        lines.append("Undecided:")
        lines += [
            f"  [{f['effective_severity']}] {f['id']}  {f['title']}"
            for f in status.undecided_findings
        ]
    if status.explicit_choice_decisions:
        lines.append("")
        lines.append("Decisions using an explicit-choice option (review these):")
        lines += [
            f"  {d['finding_id']}: {d['option_id']} by {d['decided_by']}: {d['rationale']}"
            for d in status.explicit_choice_decisions
        ]
    lines.append("")
    lines.append(
        "BLOCKED: blocking findings are undecided."
        if status.blocking_undecided
        else "OK: every blocking finding has a decision."
    )
    return "\n".join(lines)
