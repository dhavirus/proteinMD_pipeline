"""Decision coverage of a manifest: which findings are decided, which still block prep."""

from __future__ import annotations

from dataclasses import dataclass

from simprep.errors import ManifestError

SEVERITIES = ("blocking", "warn", "info")


@dataclass(frozen=True)
class DecisionStatus:
    """Counts per effective severity plus the findings that still need a decision."""

    decided: dict[str, int]
    undecided: dict[str, int]
    undecided_findings: tuple[dict, ...]
    explicit_choice_decisions: tuple[dict, ...]
    unresolved_decisions: tuple[dict, ...] = ()

    @property
    def blocking_undecided(self) -> bool:
        return self.undecided["blocking"] > 0


def decision_status(manifest: dict, include_variants: bool = True) -> DecisionStatus:
    """Summarize decision coverage of a validated manifest (needs a findings snapshot).

    With ``include_variants`` the variant_build findings of ``variant_snapshot`` count
    too (``manifest status``); prep leaves them out, since the wild type does not
    depend on them."""
    snapshot = manifest["findings_snapshot"]
    if snapshot is None:
        raise ManifestError(
            "manifest has no findings snapshot; run `simprep audit FILE --manifest M` first"
        )
    findings = snapshot["findings"] + (variant_findings(manifest) if include_variants else [])
    unresolved = unresolved_decisions(findings, manifest["decisions"])
    open_ids = {decision["finding_id"] for decision in unresolved}
    decided_ids = {d["finding_id"] for d in manifest["decisions"]} - open_ids
    undecided = tuple(f for f in findings if f["id"] not in decided_ids)
    return DecisionStatus(
        decided=_count(f for f in findings if f["id"] in decided_ids),
        undecided=_count(undecided),
        undecided_findings=tuple(sorted(undecided, key=_severity_order)),
        explicit_choice_decisions=_explicit_choices(findings, manifest["decisions"]),
        unresolved_decisions=unresolved,
    )


def variant_findings(manifest: dict) -> list[dict]:
    snapshot = manifest.get("variant_snapshot")
    return snapshot["findings"] if snapshot else []


def unresolved_decisions(findings: list[dict], decisions: list[dict]) -> tuple[dict, ...]:
    """Decisions whose option is not final (prep_action "unresolved", e.g. expert_review);
    they count as undecided."""
    unresolved = {
        (f["id"], o["id"])
        for f in findings
        for o in f["options"]
        if o.get("prep_action") == "unresolved"
    }
    return tuple(d for d in decisions if (d["finding_id"], d["option_id"]) in unresolved)


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
    if status.unresolved_decisions:
        lines.append("")
        lines.append("Decisions that are not final (they count as undecided):")
        lines += [
            f"  {d['finding_id']}: {d['option_id']} by {d['decided_by']}"
            for d in status.unresolved_decisions
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
