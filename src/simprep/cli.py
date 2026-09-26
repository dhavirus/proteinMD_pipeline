"""Command line: ``simprep audit``, ``simprep manifest init|status``, ``simprep prep`` and
``simprep variants``."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

from simprep.audit import ReportContext, build_findings_report, run_audit
from simprep.config import default_config
from simprep.knowledge import load_ruleset
from simprep.manifest import (
    attach_snapshot,
    check_input,
    config_from_manifest,
    init_manifest,
    load_manifest,
    write_json,
)
from simprep.manifest.manifest import sha256_text
from simprep.manifest.status import decision_status, render_status
from simprep.paths import KNOWLEDGE_DIR
from simprep.prep.run import RECORD_FILE, PrepRequest, run_prep
from simprep.prep.run import REPORT_FILE as PREP_REPORT_FILE
from simprep.provenance import input_info, sha256_file, simprep_provenance, utc_now
from simprep.report import render_report
from simprep.schemas import validate
from simprep.structure.parse import read_structure
from simprep.variants.run import VariantDecisionError, run_variants

FINDINGS_FILE = "findings.json"
REPORT_FILE = "report.md"
MANIFEST_FILE = "manifest.json"


def audit(args: argparse.Namespace) -> int:
    ruleset = load_ruleset(args.knowledge)
    source = input_info(args.file)
    config, manifest = default_config(ruleset.audit_defaults), None
    if args.manifest:
        manifest = load_manifest(args.manifest)
        check_input(manifest, source["sha256"])
        config = config_from_manifest(manifest, sha256_text(args.manifest.read_text()))
    structure = read_structure(args.file)
    findings = run_audit(structure, ruleset, config)
    context = ReportContext(source, ruleset, config, simprep_provenance(), utc_now())
    report = build_findings_report(structure, findings, context)
    validate(report, "findings_report")
    args.out.mkdir(parents=True, exist_ok=True)
    write_json(report, args.out / FINDINGS_FILE)
    (args.out / REPORT_FILE).write_text(render_report(report))
    if manifest is not None:
        write_json(attach_snapshot(manifest, report), args.out / MANIFEST_FILE)
    _print_summary(report, args.out)
    return 0


def _print_summary(report: dict, out: Path) -> None:
    severities = Counter(f["effective_severity"] for f in report["findings"])
    print(
        f"{len(report['findings'])} findings "
        f"(blocking {severities['blocking']}, warn {severities['warn']}, info {severities['info']})"
        f" -> {out / FINDINGS_FILE}, {out / REPORT_FILE}"
    )


def manifest_init(args: argparse.Namespace) -> int:
    if args.out.exists():
        raise FileExistsError(f"{args.out} exists; refusing to overwrite a manifest")
    ruleset = load_ruleset(args.knowledge)
    manifest = init_manifest(input_info(args.file), ruleset, (simprep_provenance(), utc_now()))
    write_json(manifest, args.out)
    print(f"manifest written to {args.out}; add regions, then run `simprep audit --manifest`")
    return 0


EXIT_BLOCKING_UNDECIDED = 3


def manifest_status(args: argparse.Namespace) -> int:
    """Exit 0 when every blocking finding is decided, 3 when any is undecided."""
    manifest = load_manifest(args.manifest)
    if args.structure:
        check_input(manifest, sha256_file(args.structure))
    status = decision_status(manifest)
    print(render_status(status))
    return EXIT_BLOCKING_UNDECIDED if status.blocking_undecided else 0


def prep(args: argparse.Namespace) -> int:
    """Apply a fully decided manifest: prepared system file(s) + prep record + report."""
    record = run_prep(PrepRequest(args.file, args.manifest, args.out, args.knowledge))
    for system in record["systems"]:
        files = ", ".join(f["path"] for f in system["files"])
        print(f"{system['name']}: {files}")
    print(
        f"{len(record['work_order'])} work-order items -> "
        f"{args.out / RECORD_FILE}, {args.out / PREP_REPORT_FILE}"
    )
    return 0


def variants(args: argparse.Namespace) -> int:
    """Pass 1: candidate findings into manifest.json (exit 3). Pass 2: build the variants."""
    try:
        outcome = run_variants(PrepRequest(args.file, args.manifest, args.out, args.knowledge))
    except VariantDecisionError as error:
        print(f"simprep: {error}", file=sys.stderr)
        return EXIT_BLOCKING_UNDECIDED
    if outcome.status == "needs_decisions":
        print(
            f"rotamer candidates written to {outcome.manifest_path} (variant_snapshot); decide "
            "them (option choose_rotamer with parameters.rotamer), then run again with it"
        )
        return EXIT_BLOCKING_UNDECIDED
    for variant in outcome.record["variants"]:
        print(f"{variant['name']}: {args.out / variant['directory']}")
    print(f"wild type: {args.out / 'wt'}; record: {args.out / 'variant_record.json'}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="simprep")
    parser.add_argument(
        "--knowledge",
        type=Path,
        default=KNOWLEDGE_DIR,
        help="knowledge-base directory (default: the repository's knowledge/)",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    audit_parser = commands.add_parser(
        "audit", help="audit a structure -> findings.json + report.md"
    )
    audit_parser.add_argument("file", type=Path)
    audit_parser.add_argument("--manifest", type=Path, help="manifest with regions and thresholds")
    audit_parser.add_argument("--out", type=Path, default=Path("out"), help="output directory")
    audit_parser.set_defaults(handler=audit)
    manifest_parser = commands.add_parser("manifest", help="manifest utilities")
    manifest_commands = manifest_parser.add_subparsers(dest="manifest_command", required=True)
    init_parser = manifest_commands.add_parser("init", help="create a manifest for a structure")
    init_parser.add_argument("file", type=Path)
    init_parser.add_argument("--out", type=Path, required=True)
    init_parser.set_defaults(handler=manifest_init)
    status_parser = manifest_commands.add_parser(
        "status",
        help="decision coverage; exit 0 = all blocking findings decided, 3 = some undecided",
    )
    status_parser.add_argument("manifest", type=Path)
    status_parser.add_argument("--structure", type=Path, help="check the input file's SHA-256")
    status_parser.set_defaults(handler=manifest_status)
    prep_parser = commands.add_parser(
        "prep", help="apply a decided manifest -> system file(s) + prep_record.json"
    )
    prep_parser.add_argument("file", type=Path)
    prep_parser.add_argument("--manifest", type=Path, required=True)
    prep_parser.add_argument("--out", type=Path, required=True, help="output directory")
    prep_parser.set_defaults(handler=prep)
    variants_parser = commands.add_parser(
        "variants", help="build the manifest's variants from the prepared wild type"
    )
    variants_parser.add_argument("file", type=Path)
    variants_parser.add_argument("--manifest", type=Path, required=True)
    variants_parser.add_argument("--out", type=Path, required=True, help="output directory")
    variants_parser.set_defaults(handler=variants)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.handler(args)
    except (OSError, ValueError, LookupError) as error:
        print(f"simprep: error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
