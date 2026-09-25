"""Compare findings computed under Pyodide with a CPython run on the same file.

Usage: python check_parity.py STRUCTURE PYODIDE_FINDINGS.json
"""

import json
import sys
from pathlib import Path

from simprep.audit import run_audit
from simprep.config import default_config
from simprep.knowledge import load_ruleset
from simprep.structure.parse import read_structure

ruleset = load_ruleset()
cpython = [
    f.to_dict()
    for f in run_audit(
        read_structure(Path(sys.argv[1])), ruleset, default_config(ruleset.audit_defaults)
    )
]
pyodide = json.loads(Path(sys.argv[2]).read_text())
print(f"findings: cpython {len(cpython)}, pyodide {len(pyodide)}; identical: {cpython == pyodide}")
sys.exit(0 if cpython == pyodide else 1)
