"""CPython side of the ADR-0002 spike: parse with gemmi, export plain JSON for Pyodide.

Usage: python export_inputs.py STRUCTURE OUT.json
"""

import json
import sys
from dataclasses import asdict
from pathlib import Path

from simprep.knowledge import load_ruleset
from simprep.structure.parse import read_structure


def export(structure_path: Path) -> dict:
    structure = read_structure(structure_path)
    ruleset = load_ruleset()
    return {
        "structure": {
            "name": structure.name,
            "residues": [
                asdict(r) | {"residue_class": r.residue_class.value} for r in structure.residues
            ],
            "links": [asdict(link) for link in structure.links],
            "modified_residues": [asdict(m) for m in structure.modified_residues],
            "unobserved_residues": [asdict(u) for u in structure.unobserved_residues],
            "annotation_categories": sorted(structure.annotation_categories),
        },
        "ruleset": {
            "version": ruleset.version,
            "sha256": ruleset.sha256,
            "rules": [asdict(rule) for rule in ruleset.rules],
            "audit_defaults": ruleset.audit_defaults,
        },
    }


if __name__ == "__main__":
    Path(sys.argv[2]).write_text(json.dumps(export(Path(sys.argv[1]))))
