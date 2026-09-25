import pytest

from simprep.rules import load_ruleset


@pytest.fixture(scope="session")
def ruleset():
    return load_ruleset()
