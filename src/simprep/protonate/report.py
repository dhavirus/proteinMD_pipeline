"""protonation_report.md: a short human-readable summary of protonation_record.json (pure)."""

from __future__ import annotations

from collections import Counter

PKA_BASIS = "pka"


def render_protonation_report(record: dict) -> str:
    lines = [
        f"# Protonation report: {record['input']['path']}",
        "",
        f"- pH {record['ph']}: {record['ph_rationale']}",
        f"- PROPKA {record['propka_version']}; hydrogens by OpenMM Modeller; ambiguity window "
        f"{record['method']['ambiguity_window_ph']} pH units; estimates within "
        f"{record['method']['unreliable_within_angstrom']} A of metal ions or non-standard "
        "residues are flagged",
        f"- knowledge base {record['knowledge_base']['version']}, "
        f"simprep {record['simprep']['version']} ({record['simprep']['git_commit']})",
        "",
    ]
    for system in record["systems"]:
        lines += _system_lines(system)
    return "\n".join(lines) + "\n"


def _system_lines(system: dict) -> list[str]:
    variants = Counter(s["variant"] for s in system["states"])
    lines = [
        f"## {system['name']} (`{system['directory']}/`, from `{system['source']}/`)",
        "",
        "States: " + ", ".join(f"{v} {n}" for v, n in sorted(variants.items())),
        "",
        "| residue | pKa | state | why |",
        "|---|---:|---|---|",
    ]
    lines += [
        f"| {s['res_name']} {s['chain']}:{s['seq_num']}{s['ins_code']} | {s['pka']} | "
        f"{s['variant']} | {s['basis']} |"
        for s in system["states"]
        if s["basis"] != PKA_BASIS
    ]
    lines += [
        "",
        *(f"- **{f['effective_severity']}** {f['id']}: {f['title']}" for f in system["findings"]),
    ]
    if system["differs_from_wild_type"]:
        lines += ["", "Differs from the wild type:"]
        lines += [
            f"- {d['res_name']} {d['chain']}:{d['seq_num']}: "
            f"{d['wild_type_variant']} -> {d['variant']}"
            for d in system["differs_from_wild_type"]
        ]
    return [*lines, ""]
