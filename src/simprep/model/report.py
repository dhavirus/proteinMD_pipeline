"""model_report.md: a short human-readable summary of model_record.json (pure)."""

from __future__ import annotations


def render_model_report(record: dict) -> str:
    minimization = record["minimization"]
    lines = [
        f"# Model report: {record['input']['path']} ({record['status']})",
        "",
        f"- input SHA-256 `{record['input']['sha256']}`",
        f"- manifest SHA-256 `{record['manifest_sha256']}`",
        f"- knowledge base {record['knowledge_base']['version']}, "
        f"simprep {record['simprep']['version']} ({record['simprep']['git_commit']})",
        f"- prepared wild type: `{record['wild_type']['directory']}/`; modelled: "
        + (f"`{record['system']['directory']}/`" if record["system"] else "not written"),
        f"- minimization: {minimization['moved_atoms']} heavy atoms moved (max "
        f"{minimization['max_displacement_angstrom']} A); loop clashes "
        f"{minimization['clashes_before']} -> {minimization['clashes_after']}",
        "",
    ]
    for loop in record["loops"]:
        lines += _loop_lines(loop)
    lines += [
        f"- **{f['effective_severity']}** {f['id']}: {f['title']}" for f in record["findings"]
    ]
    return "\n".join(lines) + "\n"


def _loop_lines(loop: dict) -> list[str]:
    sequence = "-".join(r["res_name"] for r in loop["residues"])
    lines = [
        f"## {loop['label']}: {sequence}",
        "",
        f"**{loop['note']}**",
        "",
        f"Source: {loop['source']}.",
        "",
        "| peptide bond | C-N (A) | CA-CA (A) | omega |",
        "|---|---:|---:|---:|",
    ]
    lines += [
        f"| {b['first']['res_name']}{b['first']['seq_num']}-{b['second']['res_name']}"
        f"{b['second']['seq_num']} | {b['c_n_angstrom']} | {b['ca_ca_angstrom']} | "
        f"{b['omega_degree']} |"
        for b in loop["peptide_bonds"]
    ]
    lines += ["", "| residue | phi | psi | CA improper |", "|---|---:|---:|---:|"]
    lines += [
        f"| {r['res_name']}{r['seq_num']} | {r['phi_degree']} | {r['psi_degree']} | "
        f"{'-' if r['ca_improper_degree'] is None else r['ca_improper_degree']} |"
        for r in loop["backbone"]
    ]
    return [*lines, ""]
