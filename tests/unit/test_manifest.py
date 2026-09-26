import copy
import json

import pytest

from simprep.manifest import ManifestError, check_decisions, check_input
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
    with pytest.raises(ManifestError, match="unknown finding"):
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
