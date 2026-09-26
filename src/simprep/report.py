"""Render findings.json as a human-readable Markdown report."""

from __future__ import annotations

SEVERITY_MARK = {"info": "info", "warn": "WARN", "blocking": "BLOCKING"}


def render_report(report: dict) -> str:
    """Markdown for one findings.json document."""
    sections = [_header(report), _counts(report), _summary(report)]
    sections += [_finding_section(finding) for finding in report["findings"]]
    return "\n\n".join(sections) + "\n"


def _header(report: dict) -> str:
    source, kb, tool = report["input"], report["knowledge_base"], report["simprep"]
    summary = report["structure_summary"]
    regions = report["audit_config"]["regions"]
    region_text = ", ".join(r["name"] for r in regions) if regions else "none (effective = base)"
    present = ", ".join(summary["annotation_categories_present"])
    absent = ", ".join(summary["annotation_categories_absent"])
    return "\n".join(
        [
            f"# simprep audit: {summary['name']}",
            "",
            f"- Input: `{source['path']}` ({source['format']}), SHA-256 `{source['sha256']}`",
            f"- Knowledge base: {kb['version']} (`{kb['sha256'][:12]}`)",
            f"- simprep: {tool['version']} (commit `{tool['git_commit']}`)",
            f"- Generated: {report['generated_at']}",
            f"- Regions of interest: {region_text}",
            f"- Scope: {summary['scope']}",
            f"- Annotation categories present: {present or 'none'}",
            f"- Annotation categories absent: {absent or 'none'}"
            + (" (findings that depend on them cannot be derived)" if absent else ""),
        ]
    )


def _counts(report: dict) -> str:
    lines = [
        "## Record accounting",
        "",
        "| record type | in | passed through | flagged (pending decision) | excluded |",
        "|---|---:|---:|---:|---:|",
    ]
    lines += [
        f"| {r['record_type']} | {r['in']} | {r['passed_through']} | {r['flagged']} | "
        f"{r['excluded']} |"
        for r in report["counts"]
    ]
    return "\n".join(lines)


def _summary(report: dict) -> str:
    findings = report["findings"]
    if not findings:
        return "## Findings\n\nNo findings."
    lines = [
        f"## Findings ({len(findings)})",
        "",
        "| id | severity (base -> effective) | recommended |",
        "|---|---|---|",
    ]
    lines += [
        f"| `{f['id']}` | {f['base_severity']} -> "
        f"{SEVERITY_MARK[f['effective_severity']]} | `{f['recommended_option']}` |"
        for f in findings
    ]
    return "\n".join(lines)


def _finding_section(finding: dict) -> str:
    context = finding["severity_context"]
    lines = [
        f"### {finding['title']}",
        "",
        f"`{finding['id']}`, rule `{finding['rule_id']}`; severity {finding['base_severity']}"
        f" -> **{finding['effective_severity']}** (proximity: {context['proximity']}"
        + (
            f", {context['min_distance_angstrom']} A to `{context['nearest_region']}`"
            if context["nearest_region"]
            else ""
        )
        + ")",
        "",
        "Evidence:",
    ]
    lines += [f"- {_evidence_text(item)}" for item in finding["evidence"]]
    lines += ["", "Options:"]
    for option in finding["options"]:
        marker = " **(recommended)**" if option["id"] == finding["recommended_option"] else ""
        explicit = " (explicit choice only)" if option.get("requires_explicit_choice") else ""
        lines.append(
            f"- `{option['id']}`: {option['label']}{marker}{explicit}; cost "
            f"{option['cost_hint']}. " + " ".join(option["trade_offs"])
        )
    lines += [
        "",
        f"Recommendation basis: {finding['recommendation_basis']}",
        "",
        f"Rationale: {finding['rationale']}",
    ]
    lines += [
        f"- Reference: {ref['citation']}" + (f" doi:{ref['doi']}" if "doi" in ref else "")
        for ref in finding["references"]
    ]
    lines += [f"- {flag}" for flag in finding["verify_flags"]]
    return "\n".join(lines)


def _evidence_text(item: dict) -> str:
    kind = item["type"]
    if kind == "source_record":
        fields = ", ".join(f"{k}={v}" for k, v in item["fields"].items())
        text = f"{item['key']} [{item['category']}]: {fields}"
    else:
        unit = f" {item['unit']}" if "unit" in item else ""
        text = f"{item['key']} ({kind}): {item['value']}{unit}"
    atoms = item.get("atoms")
    if atoms:
        text += " | " + " - ".join(
            f"{a['res_name']} {a['chain']}:{a['seq_num']}{a['ins_code']}.{a['atom_name']}"
            + (f"[{a['altloc']}]" if a["altloc"] else "")
            for a in atoms
        )
    if "note" in item:
        text += f" ({item['note']})"
    return text
