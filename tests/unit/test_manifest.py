import copy
import json

import pytest

from simprep.manifest import ManifestError, check_decisions, check_input, check_snapshot
from simprep.manifest.status import decision_status, render_status
from simprep.paths import SCHEMA_DIR

EXAMPLE = SCHEMA_DIR / "examples" / "manifest" / "valid_with_decision.json"


@pytest.fixture
def manifest():
    return json.loads(EXAMPLE.read_text())


def test_example_decisions_are_consistent(manifest):
    check_decisions(manifest)


def test_decision_on_unknown_finding_is_refused(manifest):
    broken = copy.deepcopy(manifest)
    broken["decisions"][0]["finding_id"] = "metals/Z:1"
    with pytest.raises(ManifestError, match="orphaned"):
        check_decisions(broken)


def test_decision_with_unknown_option_is_refused(manifest):
    broken = copy.deepcopy(manifest)
    broken["decisions"][0]["option_id"] = "magic"
    with pytest.raises(ManifestError, match="unknown option"):
        check_decisions(broken)


def test_decisions_without_snapshot_are_refused(manifest):
    broken = copy.deepcopy(manifest)
    broken["findings_snapshot"] = None
    with pytest.raises(ManifestError, match="no findings snapshot"):
        check_decisions(broken)


def test_input_hash_mismatch_is_refused(manifest):
    with pytest.raises(ManifestError, match="built for input"):
        check_input(manifest, "f" * 64)


def test_all_orphaned_decisions_are_listed_together(manifest):
    broken = copy.deepcopy(manifest)
    orphan = dict(broken["decisions"][0], finding_id="metals/Z:1")
    broken["decisions"] = [orphan, dict(orphan, finding_id="metals/Z:2")]
    with pytest.raises(ManifestError) as error:
        check_decisions(broken)
    message = str(error.value)
    assert "2 decision problem(s)" in message
    assert "metals/Z:1" in message and "metals/Z:2" in message


def test_duplicate_decisions_are_refused(manifest):
    broken = copy.deepcopy(manifest)
    broken["decisions"].append(dict(broken["decisions"][0]))
    with pytest.raises(ManifestError, match="more than one decision"):
        check_decisions(broken)


def test_example_snapshot_hash_is_valid(manifest):
    check_snapshot(manifest)


def test_modified_snapshot_is_refused(manifest):
    broken = copy.deepcopy(manifest)
    broken["findings_snapshot"]["findings"][0]["title"] = "edited"
    with pytest.raises(ManifestError, match="snapshot was modified"):
        check_snapshot(broken)


def test_status_counts_and_blocking_flag(manifest):
    status = decision_status(manifest)
    assert status.blocking_undecided
    decided = {d["finding_id"] for d in manifest["decisions"]}
    assert all(f["id"] not in decided for f in status.undecided_findings)
    assert sum(status.decided.values()) == len(decided)
    assert "BLOCKED" in render_status(status)


def test_status_ok_when_every_blocking_finding_decided(manifest):
    complete = copy.deepcopy(manifest)
    template = complete["decisions"][0]
    complete["decisions"] = [
        dict(template, finding_id=f["id"], option_id=f["recommended_option"])
        for f in complete["findings_snapshot"]["findings"]
        if f["effective_severity"] == "blocking"
    ]
    status = decision_status(complete)
    assert not status.blocking_undecided
    assert render_status(status).endswith("OK: every blocking finding has a decision.")


def test_explicit_choice_decisions_are_listed(manifest):
    explicit = copy.deepcopy(manifest)
    finding = next(
        f
        for f in explicit["findings_snapshot"]["findings"]
        if any(o.get("requires_explicit_choice") for o in f["options"])
    )
    option = next(o["id"] for o in finding["options"] if o.get("requires_explicit_choice"))
    explicit["decisions"] = [
        dict(manifest["decisions"][0], finding_id=finding["id"], option_id=option)
    ]
    status = decision_status(explicit)
    assert [d["finding_id"] for d in status.explicit_choice_decisions] == [finding["id"]]


def test_status_needs_a_snapshot(manifest):
    empty = copy.deepcopy(manifest)
    empty["findings_snapshot"], empty["decisions"] = None, []
    with pytest.raises(ManifestError, match="no findings snapshot"):
        decision_status(empty)
