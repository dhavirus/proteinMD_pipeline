"""Colab back end for variants (TASK-005): the run bookkeeping of prep_module around
:func:`simprep.variants.run.run_variants`.

One run directory holds ``wt/``, one directory per variant, ``variant_record.json`` and
``variant_report.md``. A first run on a manifest without current rotamer decisions
returns status ``needs_decisions`` and leaves ``manifest.json`` (with the candidates)
in the run directory; decide them, then start a new run with that manifest.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from modules.prep_module import RunResult, Stage, files_match, run_stage
from simprep.prep.run import RECORD_FILE as PREP_RECORD_FILE
from simprep.prep.run import PrepRequest
from simprep.variants.run import RECORD_FILE, VariantDecisionError, run_variants


def _variants(request: PrepRequest) -> str:
    try:
        return "complete" if run_variants(request).status == "built" else "needs_decisions"
    except VariantDecisionError as error:
        request.out_dir.mkdir(parents=True, exist_ok=True)
        (request.out_dir / "decisions_needed.txt").write_text(f"{error}\n")
        return "needs_decisions"


def outputs_complete(run_dir: Path) -> bool:
    """True when variant_record.json exists and every file of the wild type (from its
    prep record) and of each variant is present and matches its recorded SHA-256."""
    record_path = run_dir / RECORD_FILE
    if not record_path.exists():
        return False
    record = json.loads(record_path.read_text())
    wild_type = run_dir / record["wild_type"]["directory"]
    if not (wild_type / PREP_RECORD_FILE).exists():
        return False
    prep = json.loads((wild_type / PREP_RECORD_FILE).read_text())
    wt_files = [f for system in prep["systems"] for f in system["files"]]
    return files_match(wild_type, wt_files) and all(
        files_match(run_dir / variant["directory"], variant["files"])
        for variant in record["variants"]
    )


def run_variants_job(config: dict, now: datetime | None = None) -> RunResult:
    """Run (or resume) one variants run described by the notebook's config dict."""
    return run_stage(config, Stage(_variants, outputs_complete), now)
