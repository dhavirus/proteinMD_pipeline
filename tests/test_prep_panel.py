"""Prep on the committed panel (TASK-004): fixtures decide every finding; expected outputs
come from the fixture files (written by hand) and from the raw mmCIF atom_site table read
with gemmi.cif, never from prep's own output."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from functools import cache
from pathlib import Path

import gemmi
import pytest
import yaml

from simprep.prep.run import RECORD_FILE, PrepRequest, run_prep
from simprep.structure.model import ResidueClass
from simprep.structure.parse import read_structure
from tests.prep_helpers import audited_manifest, decided_manifest

PANEL_DIR = Path(__file__).parent / "panel"
FIXTURES = sorted((PANEL_DIR / "prep_fixtures").glob("*.yaml"))
PARAMS = [pytest.param(path.stem, id=path.stem) for path in FIXTURES]


def fixture(name: str) -> dict:
    return yaml.safe_load((PANEL_DIR / "prep_fixtures" / f"{name}.yaml").read_text())


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    """Prep each panel structure once with its fixture manifest: {name: (record, out_dir)}."""
    root = tmp_path_factory.mktemp("prep")

    @cache
    def run(name: str) -> tuple[dict, Path]:
        entry = fixture(name)
        structure = PANEL_DIR / entry["structure"]
        manifest = decided_manifest(structure, entry["decisions"], root / f"{name}.json")
        out_dir = root / name
        return run_prep(PrepRequest(structure, manifest, out_dir)), out_dir

    return run


def raw_atom_rows(path: Path) -> list[dict]:
    """atom_site rows (model 1) straight from the mmCIF, without simprep."""
    table = (
        gemmi.cif.read(str(path))
        .sole_block()
        .find(
            "_atom_site.",
            [
                "auth_asym_id",
                "auth_seq_id",
                "?pdbx_PDB_ins_code",
                "label_alt_id",
                "occupancy",
                "pdbx_PDB_model_num",
            ],
        )
    )
    rows = []
    for row in table:
        if row[5] != "1":
            continue
        icode = row[2] if row.has(2) and row[2] not in "?." else ""
        altloc = "" if row[3] in "?." else row[3]
        rows.append(
            {"residue": f"{row[0]}:{row[1]}{icode}", "altloc": altloc, "occupancy": float(row[4])}
        )
    return rows


def expected_atoms_out(rows: list[dict], entry: dict, ensemble_altloc: str | None) -> int:
    """Atom records a system must keep: everything outside excluded residues, and for a
    residue with altlocs only the chosen one (fixture altloc, ensemble altloc, or the
    highest mean occupancy with ties by id, which is what the other fixture decisions ask)."""
    excluded = excluded_residue_labels(entry, rows)
    by_residue = defaultdict(list)
    for row in rows:
        by_residue[row["residue"]].append(row)
    kept = 0
    for label, atoms in by_residue.items():
        if label in excluded:
            continue
        altloc = chosen_altloc(label, atoms, entry, ensemble_altloc)
        kept += sum(1 for a in atoms if a["altloc"] in ("", altloc))
    return kept


def chosen_altloc(label: str, atoms: list[dict], entry: dict, ensemble_altloc) -> str | None:
    expected = entry["expected"]
    if label in expected.get("ensemble_residues", []):
        return ensemble_altloc
    if label in expected["collapsed_altloc_residues"]:
        return expected["collapsed_altloc_residues"][label]
    occupancies = defaultdict(list)
    for a in atoms:
        if a["altloc"]:
            occupancies[a["altloc"]].append(a["occupancy"])
    if not occupancies:
        return None
    return min(occupancies, key=lambda alt: (-sum(occupancies[alt]) / len(occupancies[alt]), alt))


def excluded_residue_labels(entry: dict, rows: list[dict]) -> set[str]:
    spec = entry["expected"]["excluded_residues"]
    chains = {row["residue"] for row in rows if row["residue"].split(":")[0] in spec["chains"]}
    return chains | set(spec["residues"])


def counts(system: dict) -> dict:
    return {row["record_type"]: row for row in system["counts"]}


def test_every_panel_structure_has_a_prep_fixture():
    structures = {p.name.split(".")[0] for p in PANEL_DIR.glob("*.cif.gz")}
    assert {p.stem for p in FIXTURES} == structures


@pytest.mark.parametrize("name", PARAMS)
def test_fixture_decides_every_finding(name):
    entry = fixture(name)
    manifest = audited_manifest(PANEL_DIR / entry["structure"])
    assert set(entry["decisions"]) == {f["id"] for f in manifest["findings_snapshot"]["findings"]}


@pytest.mark.parametrize("name", PARAMS)
def test_accounting_equation_holds(runs, name):
    record, _ = runs(name)
    assert [s["name"] for s in record["systems"]] == fixture(name)["expected"]["systems"]
    for system in record["systems"]:
        for row in system["counts"]:
            assert row["in"] == row["passed_through"] + row["modified"] + row["excluded"], row
            assert row["out"] == row["passed_through"] + row["modified"] + row["added"], row


@pytest.mark.parametrize("name", PARAMS)
def test_excluded_counts_match_fixture(runs, name):
    record, _ = runs(name)
    expected = fixture(name)["expected"]["excluded_by_class"]
    for system in record["systems"]:
        rows = counts(system)
        assert {cls: rows[f"{cls} residues"]["excluded"] for cls in expected} == expected


@pytest.mark.parametrize("name", PARAMS)
def test_written_mmcif_rereads_to_the_record_counts(runs, name):
    record, out_dir = runs(name)
    for system in record["systems"]:
        rows = counts(system)
        structure = read_structure(out_dir / f"{system['name']}.cif")
        classes = Counter(r.residue_class for r in structure.residues)
        for cls in ResidueClass:
            assert classes[cls] == rows[f"{cls.value} residues"]["out"], cls
        assert sum(len(r.atoms) for r in structure.residues) == rows["atom records"]["out"]
        assert len(structure.links) == rows["struct_conn records"]["out"]
        assert len(structure.modified_residues) == rows["pdbx_struct_mod_residue records"]["out"]
        unobserved = [u for u in structure.unobserved_residues if u.is_polymer]
        assert len(unobserved) == rows["unobserved polymer residues (annotation)"]["out"]


@pytest.mark.parametrize("name", PARAMS)
def test_written_pdb_has_the_same_residues_and_atoms(runs, name):
    record, out_dir = runs(name)
    for system in record["systems"]:
        cif = read_structure(out_dir / f"{system['name']}.cif")
        pdb = read_structure(out_dir / f"{system['name']}.pdb")
        assert [(r.id, r.name, r.atoms) for r in pdb.residues] == [
            (r.id, r.name, r.atoms) for r in cif.residues
        ]


@pytest.mark.parametrize("name", PARAMS)
def test_atoms_out_match_raw_atom_site(runs, name):
    record, _ = runs(name)
    entry = fixture(name)
    rows = raw_atom_rows(PANEL_DIR / entry["structure"])
    for system in record["systems"]:
        assert counts(system)["atom records"]["in"] == len(rows)
        expected = expected_atoms_out(rows, entry, system["altloc"])
        assert counts(system)["atom records"]["out"] == expected


@pytest.mark.parametrize("name", PARAMS)
def test_excluded_residues_are_gone_and_nothing_else(runs, name):
    record, out_dir = runs(name)
    entry = fixture(name)
    rows = raw_atom_rows(PANEL_DIR / entry["structure"])
    expected = {row["residue"] for row in rows} - excluded_residue_labels(entry, rows)
    for system in record["systems"]:
        structure = read_structure(out_dir / f"{system['name']}.cif")
        assert {r.id.label() for r in structure.residues} == expected


@pytest.mark.parametrize("name", PARAMS)
def test_no_altlocs_left_outside_ensembles(runs, name):
    record, out_dir = runs(name)
    for system in record["systems"]:
        structure = read_structure(out_dir / f"{system['name']}.cif")
        assert [r.id.label() for r in structure.residues if r.altlocs] == []


@pytest.mark.parametrize("name", PARAMS)
def test_added_links_are_in_the_output_struct_conn(runs, name):
    record, out_dir = runs(name)
    wanted = {
        (link["ligand"], link["ligand_atom"], link["polymer"], link["polymer_atom"])
        for link in fixture(name)["expected"]["added_links"]
    }
    for system in record["systems"]:
        structure = read_structure(out_dir / f"{system['name']}.cif")
        covalent = {
            (
                k.partner1.residue.label(),
                k.partner1.atom_name,
                k.partner2.residue.label(),
                k.partner2.atom_name,
            )
            for k in structure.links
            if k.conn_type == "covale"
        }
        assert wanted <= covalent
        assert counts(system)["struct_conn records"]["added"] == len(wanted)


@pytest.mark.parametrize("name", PARAMS)
def test_work_order_stages_match_fixture(runs, name):
    record, _ = runs(name)
    stages = Counter(item["stage"] for item in record["work_order"])
    assert dict(stages) == fixture(name)["expected"]["work_order"]


def test_ensemble_systems_differ_only_at_the_ensemble_residue(runs):
    record, out_dir = runs("1FO8")
    first, second = (read_structure(out_dir / f"{s['name']}.cif") for s in record["systems"])
    differing = {
        a.id.label() for a, b in zip(first.residues, second.residues, strict=True) if a != b
    }
    assert differing == set(fixture("1FO8")["expected"]["ensemble_residues"])


def test_prep_is_deterministic(runs, tmp_path):
    record, out_dir = runs("5FQL")
    entry = fixture("5FQL")
    structure = PANEL_DIR / entry["structure"]
    manifest = decided_manifest(structure, entry["decisions"], tmp_path / "m.json")
    again = run_prep(PrepRequest(structure, manifest, tmp_path / "again"))
    for system in record["systems"]:
        for file in system["files"]:
            assert (out_dir / file["path"]).read_bytes() == (
                tmp_path / "again" / file["path"]
            ).read_bytes()
    strip = lambda doc: {k: v for k, v in doc.items() if k != "generated_at"}  # noqa: E731
    assert strip(json.loads((out_dir / RECORD_FILE).read_text())) == strip(again)
