"""Derive expected findings for a panel structure from the file's own annotations.

Independent of the detectors by construction: this script reads raw mmCIF categories
with gemmi.cif only and imports nothing from simprep. It encodes the *expected* mapping
(which element or component belongs to which rule) in its own tables, so a change in
the knowledge base that alters findings shows up as a regression-test mismatch.

Usage: python tests/panel/derive_expected.py tests/panel/5FQL.cif.gz

Writes tests/panel/expected_findings/<ID>.yaml with ``reviewed: false``. Hand-written
``curated_families`` (families that annotations cannot give, such as covalent_contacts)
and ``curated_checks`` blocks already present in that file are preserved.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import gemmi
import yaml

METAL_RULES = {
    "ZN": "metals.zinc",
    "CA": "metals.calcium",
    "MG": "metals.magnesium",
    **{e: "metals.transition" for e in ("MN", "FE", "CO", "NI", "CU", "CD")},
    "NA": "metals.alkali",
    "K": "metals.alkali",
}
STANDARD = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
COMPONENT_RULES = {
    "ALS": "nonstandard_residues.formylglycine",
    "DDZ": "nonstandard_residues.formylglycine",
}
GROUP_CONN_TYPES = {"covale", "covale_base", "covale_phosphate", "covale_sugar"}
LINK_CONN_TYPES = GROUP_CONN_TYPES | {"disulf"}


def rows(block: gemmi.cif.Block, category: str) -> list[dict]:
    table = block.find_mmcif_category(category)
    tags = [t[len(category) :] for t in table.tags]
    return [
        {
            tag: (None if row[i] in ("?", ".") else gemmi.cif.as_string(row[i]))
            for i, tag in enumerate(tags)
        }
        for row in table
    ]


def label(chain: str, seq: str, icode: str | None) -> str:
    return f"{chain}:{seq}{icode or ''}"


def first_model_atoms(block: gemmi.cif.Block) -> list[dict]:
    atoms = rows(block, "_atom_site.")
    first = atoms[0].get("pdbx_PDB_model_num")
    return [a for a in atoms if a.get("pdbx_PDB_model_num") == first]


def entity_types(block: gemmi.cif.Block) -> dict[str, str]:
    return {e["id"]: e["type"] for e in rows(block, "_entity.")}


def residues(block: gemmi.cif.Block) -> dict[str, dict]:
    """Residue label -> {name, entity type, label_seq, atoms: [atom rows]} (file order)."""
    types = entity_types(block)
    found: dict[str, dict] = {}
    for atom in first_model_atoms(block):
        key = label(atom["auth_asym_id"], atom["auth_seq_id"], atom.get("pdbx_PDB_ins_code"))
        entry = found.setdefault(
            key,
            {
                "name": atom["label_comp_id"],
                "atoms": [],
                "type": types[atom["label_entity_id"]],
                "label_seq": atom.get("label_seq_id"),
                "chain": atom["auth_asym_id"],
            },
        )
        entry["atoms"].append(atom)
    return found


def link_label(link: dict, index: int) -> str:
    p = f"ptnr{index}_"
    return label(
        link[f"{p}auth_asym_id"], link[f"{p}auth_seq_id"], link.get(f"pdbx_{p}PDB_ins_code")
    )


def atom_element(found: dict, residue: str, atom_name: str) -> str | None:
    atoms = found.get(residue, {}).get("atoms", [])
    return next((a["type_symbol"].upper() for a in atoms if a["label_atom_id"] == atom_name), None)


def metals(block, found) -> list[dict]:
    links = rows(block, "_struct_conn.")
    expected = []
    for key, res in found.items():
        elements = {a["type_symbol"].upper() for a in res["atoms"]}
        if res["type"] != "non-polymer" or len(elements) != 1:
            continue
        (element,) = elements
        if element not in METAL_RULES:
            continue
        ligands = []
        for link in links:
            if link["conn_type_id"] != "metalc":
                continue
            for i, j in ((1, 2), (2, 1)):
                if link_label(link, i) == key:
                    ligands.append(
                        {
                            "residue": link_label(link, j),
                            "atom": link[f"ptnr{j}_label_atom_id"],
                            "element": atom_element(
                                found, link_label(link, j), link[f"ptnr{j}_label_atom_id"]
                            ),
                            "distance_angstrom": float(link["pdbx_dist_value"]),
                        }
                    )
        expected.append(
            {"id": f"metals/{key}", "rule_id": METAL_RULES[element], "annotated_ligands": ligands}
        )
    return expected


def nonstandard(block, found) -> list[dict]:
    # An annotation whose component is its own parent (glycosylated ASN) marks an
    # attachment site, not a non-standard residue.
    annotated = {
        label(m["auth_asym_id"], m["auth_seq_id"], m.get("PDB_ins_code"))
        for m in rows(block, "_pdbx_struct_mod_residue.")
        if m.get("auth_comp_id") != m.get("parent_comp_id")
    }
    expected = []
    for key, res in found.items():
        if (res["type"] == "polymer" and res["name"] not in STANDARD) or key in annotated:
            rule = COMPONENT_RULES.get(res["name"], "nonstandard_residues.generic")
            expected.append({"id": f"nonstandard_residues/{key}", "rule_id": rule})
    return expected


def altlocs(found) -> list[dict]:
    expected, waters = [], 0
    for key, res in found.items():
        if not any(a.get("label_alt_id") for a in res["atoms"]):
            continue
        if res["type"] == "water":
            waters += 1
        else:
            expected.append({"id": f"altlocs/{key}"})
    if waters:
        expected.append({"id": "altlocs/water", "residues_with_altlocs": waters})
    return expected


def missing(block, found) -> list[dict]:
    observed: dict[str, list[int]] = defaultdict(list)
    for res in found.values():
        if res["type"] == "polymer" and res["label_seq"]:
            observed[res["chain"]].append(int(res["label_seq"]))
    unobserved: dict[str, list[dict]] = defaultdict(list)
    for row in rows(block, "_pdbx_unobs_or_zero_occ_residues."):
        if row.get("polymer_flag") == "Y" and row.get("occupancy_flag") == "1":
            unobserved[row["auth_asym_id"]].append(row)
    expected = []
    for chain, members in sorted(unobserved.items()):
        members.sort(key=lambda r: int(r["label_seq_id"]))
        runs: list[list[dict]] = []
        for member in members:
            if runs and int(member["label_seq_id"]) == int(runs[-1][-1]["label_seq_id"]) + 1:
                runs[-1].append(member)
            else:
                runs.append([member])
        for run in runs:
            expected.append(missing_entry(chain, run, observed[chain]))
    return expected


def missing_entry(chain: str, run: list[dict], observed: list[int]) -> dict:
    first, last = int(run[0]["label_seq_id"]), int(run[-1]["label_seq_id"])
    if not observed:
        position = "whole_chain"
    elif last < min(observed):
        position = "n_terminal"
    elif first > max(observed):
        position = "c_terminal"
    else:
        position = "internal"
    start = label(chain, run[0]["auth_seq_id"], run[0].get("PDB_ins_code"))
    end = f"{run[-1]['auth_seq_id']}{run[-1].get('PDB_ins_code') or ''}"
    return {
        "id": f"missing_residues/{start}-{end}",
        "rule_id": f"missing_residues.{position}",
        "gap_length": len(run),
    }


def unrecognized(block, found) -> list[dict]:
    claimed_metals = {m["id"].split("/", 1)[1] for m in metals(block, found)}
    claimed_polymer = {n["id"].split("/", 1)[1] for n in nonstandard(block, found)}
    candidates = [
        key
        for key, res in found.items()
        if res["type"] != "water"
        and key not in claimed_metals
        and key not in claimed_polymer
        and not (res["type"] == "polymer" and res["name"] in STANDARD)
    ]
    parent = {key: key for key in candidates}

    def root(key):
        while parent[key] != key:
            key = parent[key]
        return key

    links = rows(block, "_struct_conn.")
    for link in links:
        a, b = link_label(link, 1), link_label(link, 2)
        if link["conn_type_id"] in GROUP_CONN_TYPES and a in parent and b in parent:
            ra, rb = sorted((root(a), root(b)), key=_sort_key)
            parent[rb] = ra
    groups: dict[str, list[str]] = defaultdict(list)
    for key in candidates:
        groups[root(key)].append(key)
    expected = []
    for members in groups.values():
        members.sort(key=_sort_key)
        expected.append(
            {
                "id": f"unrecognized/group/{members[0]}",
                "components": "-".join(found[m]["name"] for m in members),
            }
        )
    grouped = set(candidates)
    for link in links:
        if (
            link["conn_type_id"] in LINK_CONN_TYPES
            and not (link_label(link, 1) in grouped or link_label(link, 2) in grouped)
            and not standard_link(link, found)
        ):
            expected.append({"id": f"unrecognized/link/{link['id']}"})
    return expected


def _sort_key(key: str) -> tuple:
    chain, rest = key.split(":")
    digits = rest.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
    return chain, int(digits), rest[len(digits) :]


def standard_link(link: dict, found: dict) -> bool:
    """Disulfide between CYS SG atoms, or a C-N bond between two polymer residues."""
    names = (link["ptnr1_label_comp_id"], link["ptnr2_label_comp_id"])
    atoms = (link["ptnr1_label_atom_id"], link["ptnr2_label_atom_id"])
    if link["conn_type_id"] == "disulf":
        return names == ("CYS", "CYS") and atoms == ("SG", "SG")
    types = [found.get(link_label(link, i), {}).get("type") for i in (1, 2)]
    return (
        link["conn_type_id"] == "covale"
        and sorted(atoms) == ["C", "N"]
        and types == ["polymer", "polymer"]
    )


def sequon_asparagines(block, found) -> dict[str, list[str]]:
    """Per chain, observed ASN residues starting an N-X-S/T sequon (X not PRO), from
    _entity_poly_seq; used to check covalent_contacts candidates independently."""
    sequences: dict[str, dict[int, set[str]]] = defaultdict(lambda: defaultdict(set))
    for row in rows(block, "_entity_poly_seq."):
        sequences[row["entity_id"]][int(row["num"])].add(row["mon_id"])
    positions: list[tuple[str, int, str]] = []
    for key, res in found.items():
        if res["type"] == "polymer" and res["label_seq"]:
            positions.append((res["atoms"][0]["label_entity_id"], int(res["label_seq"]), key))
    result: dict[str, list[str]] = defaultdict(list)
    for entity, number, key in sorted(positions, key=lambda item: _sort_key(item[2])):
        seq = sequences.get(entity, {})
        if (
            "ASN" in seq.get(number, set())
            and "PRO" not in seq.get(number + 1, {"PRO"})
            and seq.get(number + 2, set()) & {"SER", "THR"}
        ):
            result[key.split(":")[0]].append(key)
    return dict(result)


def derive(path: Path) -> dict:
    block = gemmi.cif.read(str(path)).sole_block()
    found = residues(block)
    return {
        "pdb_id": block.name.upper(),
        "file": path.name,
        "reviewed": False,
        "derived_by": "tests/panel/derive_expected.py: raw mmCIF categories (atom_site, "
        "entity, struct_conn, pdbx_struct_mod_residue, "
        "pdbx_unobs_or_zero_occ_residues, entity_poly_seq); no simprep code",
        "families": {
            "metals": metals(block, found),
            "nonstandard_residues": nonstandard(block, found),
            "altlocs": altlocs(found),
            "missing_residues": missing(block, found),
            "unrecognized": unrecognized(block, found),
        },
        "sequon_asparagines": sequon_asparagines(block, found),
    }


def main(path: Path) -> None:
    expected = derive(path)
    out = path.parent / "expected_findings" / f"{expected['pdb_id']}.yaml"
    if out.exists():
        previous = yaml.safe_load(out.read_text())
        for key in ("curated_families", "curated_checks"):
            if key in previous:
                expected[key] = previous[key]
    out.write_text(yaml.safe_dump(expected, sort_keys=False, width=100))
    print(f"wrote {out}")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
