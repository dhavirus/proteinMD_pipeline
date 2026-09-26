"""Every schema, schema example and knowledge-base rule file is validated here (CI)."""

import json
from dataclasses import replace
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator

from simprep.knowledge import load_rule_file
from simprep.paths import KNOWLEDGE_DIR, SCHEMA_DIR
from simprep.rules import RuleConsistencyError, check_rule
from simprep.schemas import SCHEMA_NAMES, SCHEMA_VERSION, load_schema, validation_errors

EXAMPLES = SCHEMA_DIR / "examples"


def load_document(path: Path):
    return (
        yaml.safe_load(path.read_text()) if path.suffix == ".yaml" else json.loads(path.read_text())
    )


def examples(prefix):
    return sorted(
        pytest.param(path, name, id=f"{name}/{path.name}")
        for name in SCHEMA_NAMES
        for path in (EXAMPLES / name).glob(f"{prefix}_*")
    )


@pytest.mark.parametrize("name", SCHEMA_NAMES)
def test_schema_is_valid_draft_2020_12_with_version(name):
    schema = load_schema(name)
    Draft202012Validator.check_schema(schema)
    assert schema["x-schema-version"] == SCHEMA_VERSION


@pytest.mark.parametrize("name", SCHEMA_NAMES)
def test_every_schema_has_valid_and_invalid_examples(name):
    files = [p.name for p in (EXAMPLES / name).iterdir()]
    assert any(f.startswith("valid_") for f in files)
    assert any(f.startswith("invalid_") for f in files)


@pytest.mark.parametrize(("path", "name"), examples("valid"))
def test_valid_examples_pass(path, name):
    assert validation_errors(load_document(path), name) == []


@pytest.mark.parametrize(("path", "name"), examples("invalid"))
def test_invalid_examples_fail(path, name):
    assert validation_errors(load_document(path), name) != []


def rule_files():
    kb = yaml.safe_load((KNOWLEDGE_DIR / "kb.yaml").read_text())
    return [KNOWLEDGE_DIR / name for name in kb["rule_files"]]


def test_kb_lists_every_rule_file():
    listed = {path.name for path in rule_files()}
    not_rules = {
        "kb.yaml",
        "audit_defaults.yaml",
        "residue_mappings.yaml",
        "side_chains.yaml",
        "relaxation.yaml",
        "modelling.yaml",
    }
    present = {p.name for p in KNOWLEDGE_DIR.glob("*.yaml")} - not_rules
    assert listed == present


@pytest.mark.parametrize("path", rule_files(), ids=lambda p: p.name)
def test_rule_file_validates(path):
    assert validation_errors(yaml.safe_load(path.read_text()), "rule") == []
    assert load_rule_file(path)


def test_ruleset_loads_with_unique_ids(ruleset):
    assert len({rule.rule_id for rule in ruleset.rules}) == len(ruleset.rules)
    assert len(ruleset.sha256) == 64


def test_every_detected_family_has_rules(ruleset):
    for family in ("metals", "nonstandard_residues", "altlocs", "missing_residues", "unrecognized"):
        assert ruleset.family(family), family


def test_unrecognized_rules_are_blocking_with_expert_default(ruleset):
    for rule in ruleset.family("unrecognized"):
        assert rule.base_severity == "blocking"
        assert rule.recommended_option["default"] == "expert_review"


def test_rules_citing_literature_have_doi_or_verify_flag(ruleset):
    for rule in ruleset.rules:
        for reference in rule.references:
            assert "doi" in reference or rule.verify_flags, rule.rule_id


def test_consistency_check_rejects_unknown_recommended_option(ruleset):
    rule = ruleset.family("altlocs")[0]
    broken = replace(rule, recommended_option={"default": "nope", "conditions": []})
    with pytest.raises(RuleConsistencyError, match="unknown option"):
        check_rule(broken, "altlocs")


def test_consistency_check_rejects_recommending_explicit_only_option(ruleset):
    rule = ruleset.family("nonstandard_residues")[-1]
    broken = replace(rule, recommended_option={"default": "revert_to_parent", "conditions": []})
    with pytest.raises(RuleConsistencyError, match="explicit human choice"):
        check_rule(broken, "nonstandard_residues")
