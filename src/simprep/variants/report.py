"""variant_report.md: a short human-readable summary of variant_record.json (pure)."""

from __future__ import annotations


def render_variant_report(record: dict) -> str:
    lines = [
        f"# Variant report: {record['input']['path']}",
        "",
        f"- input SHA-256 `{record['input']['sha256']}`",
        f"- manifest SHA-256 `{record['manifest_sha256']}`",
        f"- knowledge base {record['knowledge_base']['version']}, "
        f"simprep {record['simprep']['version']} ({record['simprep']['git_commit']})",
        f"- wild type: `{record['wild_type']['directory']}/` (prep record there)",
        "",
    ]
    for variant in record["variants"]:
        lines += _variant_lines(variant)
    return "\n".join(lines) + "\n"


def _variant_lines(variant: dict) -> list[str]:
    files = ", ".join(f"`{variant['directory']}/{f['path']}`" for f in variant["files"])
    lines = [f"## {variant['name']}", "", f"Files: {files}", ""]
    for mutation in variant["mutations"]:
        uniprot = mutation["uniprot"]
        where = (
            f"UniProt {uniprot['accession']} {uniprot['position']}"
            if uniprot
            else "no UniProt mapping"
        )
        chis = ", ".join(f"{chi:g}" for chi in mutation["chi_degree"])
        lines.append(
            f"- {mutation['chain']}:{mutation['seq_num']}{mutation['ins_code']} "
            f"{mutation['from']}>{mutation['to']} ({where}): rotamer {mutation['rotamer']} "
            f"(chi {chis}); closest contacts:"
        )
        lines += [
            f"  - {c['atom_name']} - {c['partner']['res_name']} {c['partner']['chain']}:"
            f"{c['partner']['seq_num']} {c['partner']['atom_name']}: {c['distance_angstrom']:.2f} A"
            for c in mutation["closest_contacts"]
        ]
    lines += [
        "",
        "| records (vs wild type) | in | passed | modified | excluded | added | out |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    lines += [
        f"| {c['record_type']} | {c['in']} | {c['passed_through']} | {c['modified']} | "
        f"{c['excluded']} | {c['added']} | {c['out']} |"
        for c in variant["counts"]
    ]
    return [*lines, ""]
