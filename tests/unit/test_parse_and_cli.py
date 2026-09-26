import gzip
import json

import gemmi
import pytest

from simprep.cli import main
from simprep.schemas import validation_errors
from simprep.structure.model import ResidueClass, ResidueId
from simprep.structure.parse import StructureReadError, read_structure
from tests.unit.mini_cif import mini_cif_text


@pytest.fixture
def mini_path(tmp_path):
    path = tmp_path / "mini.cif.gz"
    path.write_bytes(gzip.compress(mini_cif_text().encode()))
    return path


def test_parser_builds_model_and_annotations(mini_path):
    structure = read_structure(mini_path)
    assert structure.name == "MINI"
    classes = {r.id: r.residue_class for r in structure.residues}
    assert classes[ResidueId("A", 3)] is ResidueClass.POLYMER
    assert classes[ResidueId("A", 901)] is ResidueClass.NONPOLYMER
    assert classes[ResidueId("A", 1002)] is ResidueClass.WATER
    assert structure.residue(ResidueId("A", 1002)).altlocs == ("A", "B")
    assert structure.residue(ResidueId("A", 6)).label_seq == 6
    assert [link.conn_id for link in structure.links] == ["metalc1"]
    assert structure.modified_residues[0].parent_res_name == "MET"
    assert [u.label_seq for u in structure.unobserved_residues] == [1, 5]
    assert structure.annotation_categories == {
        "struct_conn",
        "pdbx_struct_mod_residue",
        "pdbx_unobs_or_zero_occ_residues",
    }


def test_multi_model_input_is_refused(tmp_path):
    path = tmp_path / "two.cif"
    path.write_text(mini_cif_text(models=2))
    with pytest.raises(StructureReadError, match="2 models"):
        read_structure(path)


def test_pdb_input_reports_missing_annotation_categories(mini_path, tmp_path):
    pdb_path = tmp_path / "mini.pdb"
    gemmi.read_structure(str(mini_path)).write_pdb(str(pdb_path))
    structure = read_structure(pdb_path)
    assert "pdbx_unobs_or_zero_occ_residues" not in structure.annotation_categories
    assert structure.residue(ResidueId("A", 901)).residue_class is ResidueClass.NONPOLYMER


def run_cli(*args):
    return main([str(a) for a in args])


def test_audit_cli_writes_valid_findings_and_report(mini_path, tmp_path):
    out = tmp_path / "out"
    assert run_cli("audit", mini_path, "--out", out) == 0
    report = json.loads((out / "findings.json").read_text())
    assert validation_errors(report, "findings_report") == []
    ids = {f["id"] for f in report["findings"]}
    assert ids == {
        "metals/A:901",
        "nonstandard_residues/A:3",
        "altlocs/water",
        "missing_residues/A:1-1",
        "missing_residues/A:5-5",
        "unrecognized/group/A:902",
    }
    for row in report["counts"]:
        assert row["in"] == row["passed_through"] + row["flagged"] + row["excluded"]
    text = (out / "report.md").read_text()
    assert "# simprep audit: MINI" in text and "Record accounting" in text


def test_manifest_round_trip_and_escalation(mini_path, tmp_path):
    manifest_path = tmp_path / "manifest.json"
    assert run_cli("manifest", "init", mini_path, "--out", manifest_path) == 0
    manifest = json.loads(manifest_path.read_text())
    manifest["regions"] = [
        {
            "name": "zinc_site",
            "description": "test",
            "residues": [{"chain": "A", "seq_num": 901, "ins_code": ""}],
        }
    ]
    manifest_path.write_text(json.dumps(manifest))
    out = tmp_path / "out"
    assert run_cli("audit", mini_path, "--manifest", manifest_path, "--out", out) == 0
    report = json.loads((out / "findings.json").read_text())
    metal = next(f for f in report["findings"] if f["id"] == "metals/A:901")
    assert metal["effective_severity"] == "blocking"
    assert report["audit_config"]["manifest_sha256"] is not None
    updated = json.loads((out / "manifest.json").read_text())
    assert validation_errors(updated, "manifest") == []
    assert len(updated["findings_snapshot"]["findings"]) == len(report["findings"])


def test_manifest_for_other_input_is_refused(mini_path, tmp_path, capsys):
    manifest_path = tmp_path / "manifest.json"
    run_cli("manifest", "init", mini_path, "--out", manifest_path)
    other = tmp_path / "other.cif"
    other.write_text(mini_cif_text().replace("MINI", "MINJ"))
    assert run_cli("audit", other, "--manifest", manifest_path, "--out", tmp_path / "o") == 2
    assert "built for input SHA-256" in capsys.readouterr().err


def test_audit_is_deterministic(mini_path, tmp_path):
    reports = []
    for name in ("a", "b"):
        run_cli("audit", mini_path, "--out", tmp_path / name)
        reports.append(json.loads((tmp_path / name / "findings.json").read_text()))
    assert reports[0]["findings"] == reports[1]["findings"]
    assert reports[0]["counts"] == reports[1]["counts"]
