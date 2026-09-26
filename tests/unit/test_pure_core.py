"""The detector core must import without compiled or I/O-only dependencies (ADR-0002)."""

import subprocess
import sys

BLOCKED = ("gemmi", "yaml", "jsonschema", "referencing")
PURE_MODULES = (
    "simprep.audit",
    "simprep.report",
    "simprep.severity",
    "simprep.rules",
    "simprep.decisions",
    "simprep.prep.plan",
    "simprep.prep.apply",
    "simprep.prep.record",
    "simprep.variants.build",
    "simprep.variants.check",
    "simprep.variants.findings",
    "simprep.variants.validate",
    "simprep.variants.apply",
    "simprep.variants.report",
    "simprep.variants.relaxation",
    "simprep.relax.shell",
    "simprep.relax.openmm_run",
)

PROBE = f"""
import builtins
real_import = builtins.__import__
def guarded(name, *args, **kwargs):
    if name.split(".")[0] in {BLOCKED!r}:
        raise ImportError("pure core imported " + name)
    return real_import(name, *args, **kwargs)
builtins.__import__ = guarded
import importlib
for module in {PURE_MODULES!r}:
    importlib.import_module(module)
from simprep.detectors import load_detectors
load_detectors()
"""


def test_detector_core_has_no_compiled_dependencies():
    result = subprocess.run([sys.executable, "-c", PROBE], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
