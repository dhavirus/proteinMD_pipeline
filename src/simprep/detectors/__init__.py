"""Detectors: one module per rule family, each a pure function
``(Structure, RuleSet, AuditConfig) -> list[Finding]``, registered by family name.

The ``unrecognized`` pass is not a family detector: it runs last and reports whatever
the family detectors did not claim.
"""

from collections.abc import Callable

from simprep.config import AuditConfig
from simprep.findings import Finding
from simprep.rules import RuleSet
from simprep.structure.model import Structure

Detector = Callable[[Structure, RuleSet, AuditConfig], list[Finding]]

DETECTORS: dict[str, Detector] = {}


def register(family: str) -> Callable[[Detector], Detector]:
    """Decorator registering ``detector`` for ``family``."""

    def decorate(detector: Detector) -> Detector:
        if family in DETECTORS:
            raise ValueError(f"detector for family {family!r} already registered")
        DETECTORS[family] = detector
        return detector

    return decorate


def load_detectors() -> dict[str, Detector]:
    """Import every family module (which registers itself) and return the registry."""
    from simprep.detectors import (  # noqa: F401
        altlocs,
        metals,
        missing_residues,
        nonstandard_residues,
    )

    return dict(DETECTORS)
