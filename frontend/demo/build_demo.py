"""Generate (or check) the demo findings file for the review page.

    python frontend/demo/build_demo.py          # regenerate frontend/demo/5FQL.findings.json
    python frontend/demo/build_demo.py --check  # exit 1 if it is stale

The demo must never be hand-edited: it is the output of `simprep audit` on the committed
panel file, run from the repository root so the recorded input path is repo-relative.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from simprep.audit import ReportContext, build_findings_report, run_audit
from simprep.cli import input_info
from simprep.config import default_config
from simprep.knowledge import load_ruleset
from simprep.paths import REPO_ROOT
from simprep.provenance import simprep_provenance, utc_now
from simprep.schemas import validate
from simprep.structure.parse import read_structure

STRUCTURE = Path("tests/panel/5FQL.cif.gz")
OUTPUT = REPO_ROOT / "frontend" / "demo" / "5FQL.findings.json"
# Fields that must match for the demo to be current (provenance and timestamps may differ).
COMPARED = (
    "schema_version",
    "knowledge_base",
    "input",
    "audit_config",
    "structure_summary",
    "counts",
    "findings",
)


def build() -> dict:
    ruleset = load_ruleset()
    config = default_config(ruleset.audit_defaults)
    structure = read_structure(REPO_ROOT / STRUCTURE)
    source = input_info(REPO_ROOT / STRUCTURE) | {"path": STRUCTURE.as_posix()}
    context = ReportContext(source, ruleset, config, simprep_provenance(), utc_now())
    report = build_findings_report(structure, run_audit(structure, ruleset, config), context)
    validate(report, "findings_report")
    return report


def stale_fields(current: dict, committed: dict) -> list[str]:
    return [key for key in COMPARED if current[key] != committed.get(key)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="fail if the committed demo is stale")
    report = build()
    if parser.parse_args().check:
        stale = stale_fields(report, json.loads(OUTPUT.read_text()))
        if stale:
            print(
                f"{OUTPUT} is stale ({', '.join(stale)}); run {Path(__file__).name}",
                file=sys.stderr,
            )
            return 1
        return 0
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
