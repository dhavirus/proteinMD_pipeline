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
        _protocol_line(record["relaxation_protocol"]),
        "",
    ]
    for wild_type in record["relaxed_wild_types"]:
        lines += [f"## {wild_type['name']}", "", *_relaxation_lines(wild_type), ""]
    for variant in record["variants"]:
        lines += _variant_lines(variant)
    return "\n".join(lines) + "\n"


def _protocol_line(protocol: dict | None) -> str:
    if not protocol or not protocol["enabled"]:
        return "- relaxation: off (variants are rigid builds)"
    return (
        f"- relaxation: {protocol['engine']} {', '.join(protocol['force_field_files'])}; shell "
        f"{protocol['mobile_shell_angstrom']} A; restraint {protocol['restraint_kj_per_mol_nm2']} "
        f"kJ/mol/nm^2; platform {protocol['platform']}; seed {protocol['seed']}"
    )


def _relaxation_lines(system: dict) -> list[str]:
    relaxation = system.get("relaxation")
    if relaxation is None:
        return []
    lines = [
        f"Relaxation: {relaxation['moved_atoms']} heavy atoms in "
        f"{len(relaxation['moved_residues'])} residues moved (max "
        f"{relaxation['max_displacement_angstrom']} A, rms "
        f"{relaxation['rms_displacement_angstrom']} A); clashes at the site "
        f"{relaxation['clashes_before']} -> {relaxation['clashes_after']}.",
    ]
    if relaxation["frozen"]:
        frozen = ", ".join(
            f"{f['chain']}:{f['seq_num']} (near {f['near']})" for f in relaxation["frozen"]
        )
        lines.append(f"Held fixed (template-less chemistry within the cutoff): {frozen}.")
    if relaxation["left_out"]:
        left = ", ".join(
            f"{r['res_name']} {r['chain']}:{r['seq_num']}" for r in relaxation["left_out"]
        )
        lines.append(f"Left out of the calculation: {left}.")
    lines += [f"WARNING {f['id']}: {f['title']}" for f in system.get("findings", [])]
    return lines


def _variant_lines(variant: dict) -> list[str]:
    files = ", ".join(f"`{variant['directory']}/{f['path']}`" for f in variant["files"])
    lines = [f"## {variant['name']}", "", f"Files: {files}", "", *_relaxation_lines(variant), ""]
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
            f"(chi {chis}); closest contacts as built, before relaxation:"
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
