"""Load and validate the knowledge base (versioned YAML rule files) from disk."""

from __future__ import annotations

import hashlib
from pathlib import Path

import yaml

from simprep.paths import KNOWLEDGE_DIR, require_dir
from simprep.rules import Rule, RuleSet, check_rule, check_unique_ids
from simprep.schemas import validate

KB_MANIFEST = "kb.yaml"
AUDIT_DEFAULTS = "audit_defaults.yaml"


def load_rule_file(path: Path) -> tuple[Rule, ...]:
    """Parse, schema-validate and consistency-check one rule file."""
    data = yaml.safe_load(path.read_text())
    validate(data, "rule")
    rules = tuple(Rule.from_dict(rule) for rule in data["rules"])
    for rule in rules:
        check_rule(rule, data["family"])
    return rules


def load_ruleset(knowledge_dir: Path = KNOWLEDGE_DIR) -> RuleSet:
    """Load every rule file listed in ``kb.yaml``; rule ids must be unique."""
    require_dir(knowledge_dir)
    kb = yaml.safe_load((knowledge_dir / KB_MANIFEST).read_text())
    files = [knowledge_dir / KB_MANIFEST, knowledge_dir / AUDIT_DEFAULTS] + [
        knowledge_dir / name for name in kb["rule_files"]
    ]
    rules = tuple(rule for path in files[2:] for rule in load_rule_file(path))
    check_unique_ids(rules)
    return RuleSet(
        version=kb["kb_version"],
        sha256=_hash_files(files),
        rules=rules,
        audit_defaults=yaml.safe_load((knowledge_dir / AUDIT_DEFAULTS).read_text()),
    )


def _hash_files(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.name.encode() + b"\0" + path.read_bytes() + b"\0")
    return digest.hexdigest()
