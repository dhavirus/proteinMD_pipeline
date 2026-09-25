import pytest

from simprep.knowledge import load_ruleset


@pytest.fixture(scope="session")
def ruleset():
    return load_ruleset()
