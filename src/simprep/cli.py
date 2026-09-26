"""Command line: ``simprep audit`` and ``simprep manifest init``."""

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
from simprep.paths import KNOWLEDGE_DIR
from simprep.provenance import sha256_file, simprep_provenance, utc_now
from simprep.report import render_report
from simprep.schemas import validate
from simprep.structure.parse import read_structure

FINDINGS_FILE = "findings.json"
REPORT_FILE = "report.md"
MANIFEST_FILE = "manifest.json"


FORMAT_BY_EXTENSION = {".cif": "mmcif", ".mmcif": "mmcif", ".pdb": "pdb", ".ent": "pdb"}


def input_info(path: Path) -> dict:
    """Path, SHA-256 and format (from the extension, ignoring a trailing .gz)."""
    suffixes = [s.lower() for s in path.suffixes if s.lower() != ".gz"]
    fmt = FORMAT_BY_EXTENSION.get(suffixes[-1] if suffixes else "")
    if fmt is None:
        raise ValueError(f"{path}: unsupported extension; expected .cif/.mmcif/.pdb/.ent[.gz]")
    return {"path": str(path), "sha256": sha256_file(path), "format": fmt}


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
