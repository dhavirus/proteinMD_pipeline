"""Errors shared by the pure core and the I/O edge (no dependencies)."""


class ManifestError(ValueError):
    """The manifest is inconsistent with its input file or its own findings."""
