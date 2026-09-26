"""Load and validate the knowledge base (versioned YAML rule files) from disk."""

from __future__ import annotations

import hashlib
from pathlib import Path

import yaml

from simprep.paths import KNOWLEDGE_DIR, require_dir
from simprep.rules import Rule, RuleConsistencyError, RuleSet, check_rule, check_unique_ids
from simprep.schemas import validate

KB_MANIFEST = "kb.yaml"
AUDIT_DEFAULTS = "audit_defaults.yaml"
RESIDUE_MAPPINGS = "residue_mappings.yaml"
SIDE_CHAINS = "side_chains.yaml"
DATA_FILES = (RESIDUE_MAPPINGS, SIDE_CHAINS)


def load_rule_file(path: Path) -> tuple[Rule, ...]:
    """Parse, schema-validate and consistency-check one rule file."""
    data = yaml.safe_load(path.read_text())
    validate(data, "rule")
    rules = tuple(Rule.from_dict(rule) for rule in data["rules"])
    for rule in rules:
        check_rule(rule, data["family"])
    return rules


def load_ruleset(knowledge_dir: Path = KNOWLEDGE_DIR) -> RuleSet:
    """Load every rule file listed in ``kb.yaml`` plus the residue mappings; rule ids
    must be unique. The knowledge-base hash covers all of these files."""
    require_dir(knowledge_dir)
    kb = yaml.safe_load((knowledge_dir / KB_MANIFEST).read_text())
    rule_files = [knowledge_dir / name for name in kb["rule_files"]]
    rules = tuple(rule for path in rule_files for rule in load_rule_file(path))
    check_unique_ids(rules)
    fixed = [knowledge_dir / KB_MANIFEST, knowledge_dir / AUDIT_DEFAULTS]
    return RuleSet(
        version=kb["kb_version"],
        sha256=_hash_files([*fixed, *rule_files, *(knowledge_dir / name for name in DATA_FILES)]),
        rules=rules,
        audit_defaults=yaml.safe_load((knowledge_dir / AUDIT_DEFAULTS).read_text()),
        residue_mappings=load_residue_mappings(knowledge_dir / RESIDUE_MAPPINGS),
        side_chains=load_side_chains(knowledge_dir / SIDE_CHAINS),
    )


def load_residue_mappings(path: Path) -> tuple[dict, ...]:
    """Schema-validated revert_to_parent atom mappings; one entry per source component."""
    data = yaml.safe_load(path.read_text())
    validate(data, "residue_mappings")
    sources = [mapping["from"] for mapping in data["mappings"]]
    duplicates = sorted({name for name in sources if sources.count(name) > 1})
    if duplicates:
        raise RuleConsistencyError(f"{path.name}: more than one mapping for {duplicates}")
    return tuple(data["mappings"])


def load_side_chains(path: Path) -> dict:
    """Schema-validated side-chain building data per residue type (variants)."""
    data = yaml.safe_load(path.read_text())
    validate(data, "side_chains")
    return data["residues"]


def _hash_files(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.name.encode() + b"\0" + path.read_bytes() + b"\0")
    return digest.hexdigest()
