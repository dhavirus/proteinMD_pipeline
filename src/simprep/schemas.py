"""Load the JSON Schemas and validate documents against them."""

from __future__ import annotations

import json
from functools import cache
from pathlib import Path

from jsonschema import Draft202012Validator
from referencing import Registry, Resource

from simprep.paths import SCHEMA_DIR, require_dir

SCHEMA_NAMES = ("finding", "rule", "manifest", "findings_report")
SCHEMA_VERSION = "0.1.0"


class SchemaValidationError(ValueError):
    """A document does not conform to its schema."""


def schema_path(name: str) -> Path:
    return require_dir(SCHEMA_DIR) / f"{name}.schema.json"


@cache
def load_schema(name: str) -> dict:
    """Return the parsed schema called ``name`` (one of :data:`SCHEMA_NAMES`)."""
    if name not in SCHEMA_NAMES:
        raise KeyError(f"unknown schema {name!r}; expected one of {SCHEMA_NAMES}")
    return json.loads(schema_path(name).read_text())


@cache
def _registry() -> Registry:
    resources = [
        (load_schema(name)["$id"], Resource.from_contents(load_schema(name)))
        for name in SCHEMA_NAMES
    ]
    return Registry().with_resources(resources)


def validator(name: str) -> Draft202012Validator:
    return Draft202012Validator(load_schema(name), registry=_registry())


def validation_errors(document: object, name: str) -> list[str]:
    """Human-readable validation errors of ``document`` against schema ``name``."""
    errors = sorted(validator(name).iter_errors(document), key=lambda error: list(error.path))
    return [f"{'/'.join(map(str, error.path)) or '<root>'}: {error.message}" for error in errors]


def validate(document: object, name: str) -> None:
    """Raise :class:`SchemaValidationError` listing every error if ``document`` is invalid."""
    errors = validation_errors(document, name)
    if errors:
        raise SchemaValidationError(f"invalid {name} document:\n  " + "\n  ".join(errors))
