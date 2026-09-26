"""`simprep variants` on 5FQL (TASK-005): WT vs R468Q vs R468W, rigid builds (the
manifest switches relaxation off; test_relaxation_panel.py covers TASK-006). Expected
outcomes come from tests/panel/variant_fixtures/5FQL.yaml (written by hand, with
placements made independently of simprep's builder)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from simprep.cli import main
from simprep.knowledge import load_ruleset
from simprep.prep.apply import apply_plan
from simprep.prep.plan import build_plan
from simprep.structure.model import Atom, ResidueClass, ResidueId
from simprep.structure.parse import read_structure
from simprep.variants.check import ClashCriterion, clashes, contacts, neighbours
from tests.prep_helpers import decided_manifest

PANEL_DIR = Path(__file__).parent / "panel"
FIXTURE = yaml.safe_load((PANEL_DIR / "variant_fixtures" / "5FQL.yaml").read_text())
STRUCTURE = PANEL_DIR / FIXTURE["structure"]
PREP_DECISIONS = (
    yaml.safe_load((PANEL_DIR / "prep_fixtures" / FIXTURE["prep_fixture"]).read_text())["decisions"]
    | FIXTURE["prep_overrides"]
)
SITE = ResidueId("A", 468)


def cli(*args) -> int:
    return main([str(a) for a in args])


def manifest_with(variants: list[dict], path: Path, relaxation: bool = False) -> Path:
    decided_manifest(STRUCTURE, PREP_DECISIONS, path)
    document = json.loads(path.read_text())
    document["variants"] = variants
    document["relaxation"] = {**load_ruleset().relaxation, "enabled": relaxation}
    path.write_text(json.dumps(document))
    return path


def variant(name: str) -> dict:
    return next(v for v in FIXTURE["variants"] if v["name"] == name)


def decide(manifest_path: Path, choices: dict, out: Path) -> Path:
    """Add rotamer decisions ({variant: (option, rotamer)}) to a pass-1 manifest."""
    document = json.loads(manifest_path.read_text())
    for name, (option, rotamer) in choices.items():
        decision = {
            "finding_id": f"variant_build/{name}/A:468",
            "option_id": option,
            "rationale": "test",
            "decided_by": "test",
            "timestamp": "2026-09-26T00:00:00Z",
        }
        if rotamer:
            decision["parameters"] = {"rotamer": rotamer}
        document["decisions"].append(decision)
    out.write_text(json.dumps(document))
    return out


@pytest.fixture(scope="module")
def pass_one(tmp_path_factory):
    root = tmp_path_factory.mktemp("variants")
    manifest = manifest_with(FIXTURE["variants"], root / "m.json")
    code = cli("variants", STRUCTURE, "--manifest", manifest, "--out", root / "pass1")
    return code, root / "pass1" / "manifest.json", root


def variant_findings(manifest_path: Path) -> dict:
    snapshot = json.loads(manifest_path.read_text())["variant_snapshot"]
    return {f["id"].split("/")[1]: f for f in snapshot["findings"]}


def evidence(finding: dict, key: str):
    return next(item["value"] for item in finding["evidence"] if item["key"] == key)


def test_pass_one_stops_for_rotamer_decisions(pass_one):
    code, manifest_path, _ = pass_one
    assert code == 3
    assert set(variant_findings(manifest_path)) == {"R468Q", "R468W"}


@pytest.mark.parametrize("name", ["R468Q", "R468W"])
def test_findings_match_the_hand_written_expectations(pass_one, name):
    finding = variant_findings(pass_one[1])[name]
    expected, site = FIXTURE["expected"][name], FIXTURE["expected"]["site"]
    assert finding["effective_severity"] == "blocking"
    assert evidence(finding, "all_candidates_clash") is expected["all_candidates_clash"]
    assert finding["recommended_option"] == expected["recommended_option"]
    uniprot = site["uniprot"]
    assert evidence(finding, "uniprot_position") == f"{uniprot['accession']} {uniprot['position']}"
    assert evidence(finding, "mutation").startswith(f"{site['residue']} {site['res_name']}>")


@pytest.fixture(scope="module")
def wild_type(pass_one):
    document = json.loads(pass_one[1].read_text())
    structure = read_structure(STRUCTURE)
    plan = build_plan(structure, document, load_ruleset())
    (system,) = apply_plan(structure, plan)
    return system.structure


@pytest.mark.parametrize(("residue_type", "clash_expected"), [("TRP", True), ("GLN", False)])
def test_independent_pdbfixer_placements_are_judged_as_expected(
    wild_type, residue_type, clash_expected
):
    criterion = ClashCriterion.from_matcher(load_ruleset().family("variant_build")[0].matcher)
    placement = FIXTURE["pdbfixer_placements"][residue_type]
    atoms = tuple(Atom(name, name[0], tuple(xyz), 1.0, 60.0) for name, xyz in placement.items())
    found = clashes(contacts(atoms, neighbours(wild_type, SITE), criterion), criterion)
    assert bool(found) is clash_expected
    if clash_expected:
        assert ResidueId("A", 470) in {c.residue.id for c in found}


def test_clashing_rotamers_are_built_and_open_choices_refused(pass_one, capsys):
    """Clashes only flag (maintainer, TASK-005 review): a clashing rotamer is built."""
    _, manifest_path, root = pass_one
    forced = decide(
        manifest_path,
        {"R468Q": ("choose_rotamer", "mm-40"), "R468W": ("choose_rotamer", "m95")},
        root / "forced.json",
    )
    assert cli("variants", STRUCTURE, "--manifest", forced, "--out", root / "forced") == 0
    assert (root / "forced" / "R468W" / "system.cif").is_file()
    opened = decide(
        manifest_path,
        {"R468Q": ("choose_rotamer", "mm-40"), "R468W": ("expert_review", None)},
        root / "open.json",
    )
    assert cli("variants", STRUCTURE, "--manifest", opened, "--out", root / "open") == 3
    assert "expert_review is not a final decision" in capsys.readouterr().err
    assert cli("manifest", "status", opened) == 3


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """R468Q alone: pass 1, decide the recommended rotamer, pass 2."""
    root = tmp_path_factory.mktemp("r468q")
    manifest = manifest_with([variant("R468Q")], root / "m.json")
    assert cli("variants", STRUCTURE, "--manifest", manifest, "--out", root / "p1") == 3
    finding = variant_findings(root / "p1" / "manifest.json")["R468Q"]
    rotamer = evidence(finding, "recommended_rotamer")
    decided = decide(
        root / "p1" / "manifest.json", {"R468Q": ("choose_rotamer", rotamer)}, root / "decided.json"
    )
    assert cli("variants", STRUCTURE, "--manifest", decided, "--out", root / "out") == 0
    record = json.loads((root / "out" / "variant_record.json").read_text())
    return record, root / "out", decided, root


def test_variant_differs_from_wild_type_only_at_the_mutated_residue(built):
    _, out, _, _ = built
    wt = read_structure(out / "wt" / "system.cif")
    mutant = read_structure(out / "R468Q" / "system.cif")
    differing = {a.id for a, b in zip(wt.residues, mutant.residues, strict=True) if a != b}
    assert differing == {SITE}
    residue = mutant.residue(SITE)
    assert residue.name == "GLN" and all(atom.is_heavy for atom in residue.atoms)
    sequence = dict(mutant.polymer_sequences)["A"]
    assert sequence[residue.label_seq - 1] == "GLN"
    assert dict(wt.polymer_sequences)["A"][residue.label_seq - 1] == "ARG"


def test_accounting_holds_and_matches_the_written_files(built):
    record, out, _, _ = built
    (entry,) = record["variants"]
    rows = {row["record_type"]: row for row in entry["counts"]}
    for row in entry["counts"]:
        assert row["in"] == row["passed_through"] + row["modified"] + row["excluded"], row
        assert row["out"] == row["passed_through"] + row["modified"] + row["added"], row
    mutant = read_structure(out / "R468Q" / "system.cif")
    for cls in ResidueClass:
        count = sum(1 for r in mutant.residues if r.residue_class is cls)
        assert count == rows[f"{cls.value} residues"]["out"]
    assert sum(len(r.atoms) for r in mutant.residues) == rows["atom records"]["out"]
    assert rows["atom records"]["added"] == 2  # OE1, NE2: names Arg does not have
    assert rows["polymer residues"]["modified"] == 1


def test_record_carries_rotamer_uniprot_and_protonation_item(built):
    record, _, _, _ = built
    (mutation,) = record["variants"][0]["mutations"]
    assert mutation["uniprot"] == FIXTURE["expected"]["site"]["uniprot"]
    assert mutation["from"] == "ARG" and mutation["to"] == "GLN"
    stages = [item["stage"] for item in record["variants"][0]["work_order"]]
    assert stages.count("protonation") == 1


def test_wild_type_equals_plain_prep_and_runs_are_deterministic(built, tmp_path):
    record, out, decided, _ = built
    assert cli("prep", STRUCTURE, "--manifest", decided, "--out", tmp_path / "prep") == 0
    for name in ("system.cif", "system.pdb"):
        assert (out / "wt" / name).read_bytes() == (tmp_path / "prep" / name).read_bytes()
    assert cli("variants", STRUCTURE, "--manifest", decided, "--out", tmp_path / "again") == 0
    for name in ("system.cif", "system.pdb"):
        assert (out / "R468Q" / name).read_bytes() == (
            tmp_path / "again" / "R468Q" / name
        ).read_bytes()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ({"from": "LYS", "to": "GLN"}, "the structure has ARG there"),
        ({"from": "ARG", "to": "ALA"}, "no side-chain data for ALA"),
    ],
)
def test_bad_mutations_are_refused_with_actionable_messages(tmp_path, capsys, mutation, message):
    bad = {
        **variant("R468Q"),
        "name": "BAD",
        "mutations": [{"chain": "A", "seq_num": 468, "ins_code": "", **mutation}],
    }
    manifest = manifest_with([bad], tmp_path / "m.json")
    assert cli("variants", STRUCTURE, "--manifest", manifest, "--out", tmp_path / "out") == 2
    assert message in capsys.readouterr().err
