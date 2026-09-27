"""The modelling stage at the I/O edge (TASK-007, ADR-0007).

prep (select) -> model (add) -> variants (build) -> relax. ``model_wild_type`` builds
every gap decided ``model_loop`` on the prepared wild type: PDBFixer places the residues,
OpenMM minimizes the loops (free) and their flanks (restrained), and the result is
checked. ``write_model`` writes ``wt_modelled/`` (only if no check is blocking),
``model_record.json`` and ``model_report.md``; a blocking check then raises LoopRejected.
``simprep variants`` calls both in-process; ``simprep model`` runs them for the wild
type alone.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from importlib.metadata import version as installed_version
from pathlib import Path

from simprep.findings import residue_ref
from simprep.manifest import config_from_manifest, write_json
from simprep.model.checks import (
    NO_EXPERIMENTAL_SUPPORT,
    backbone,
    gap_findings,
    loop_clashes,
    peptide_bonds,
)
from simprep.model.checks import criterion as clash_criterion
from simprep.model.loops import (
    APPLIED_OPTIONS,
    Gap,
    ModelError,
    gaps_to_model,
    loop_shell,
    torsion_restraints,
    with_loops,
)
from simprep.model.placement import place_loops
from simprep.model.report import render_model_report
from simprep.prep.apply import apply_plan
from simprep.prep.plan import PrepPlan, build_plan
from simprep.prep.record import system_counts
from simprep.prep.run import RECORD_FILE as PREP_RECORD_FILE
from simprep.prep.run import PrepInputs, PrepRequest, load_inputs, write_prep
from simprep.prep.write import write_system
from simprep.provenance import sha256_file, simprep_provenance, utc_now
from simprep.relax.openmm_run import MOVED_ANGSTROM, RelaxOutcome, relax
from simprep.schemas import validate
from simprep.severity import apply_context
from simprep.structure.model import Structure
from simprep.variants.relaxation import Relaxed, relaxation_entry

MODELLED_DIR = "wt_modelled"
WILD_TYPE_DIR = "wt"
RECORD_FILE = "model_record.json"
REPORT_FILE = "model_report.md"
SYSTEM_NAME = "system"
DIGITS = 3


class LoopRejected(ModelError):
    """A modelled loop failed a blocking check; the record lists why, no system was written."""


@dataclass(frozen=True)
class Modelled:
    """The wild type with its loops built, and what the stage did."""

    structure: Structure
    gaps: tuple[Gap, ...]
    relaxed: Relaxed
    findings: list[dict]
    protocol: dict

    @property
    def rejected(self) -> bool:
        return any(f["effective_severity"] == "blocking" for f in self.findings)


def run_model(request: PrepRequest) -> dict:
    """``simprep model``: prepare the wild type, build its loops, write everything."""
    inputs = load_inputs(request)
    plan = build_plan(inputs.structure, inputs.manifest, inputs.ruleset)
    modelled = model_wild_type(inputs, plan, single_system(inputs, plan))
    if modelled is None:
        raise ModelError("the manifest decides no gap model_loop; there is nothing to model")
    wt_record = write_prep(inputs, plan, request.out_dir / WILD_TYPE_DIR)
    return write_model(modelled, (inputs, wt_record), request.out_dir)


def single_system(inputs: PrepInputs, plan: PrepPlan) -> Structure:
    systems = apply_plan(inputs.structure, plan)
    if len(systems) != 1:
        raise ModelError("modelling v0.1 builds on one wild-type system; choose one altloc first")
    return systems[0].structure


def model_wild_type(inputs: PrepInputs, plan: PrepPlan, wild_type: Structure) -> Modelled | None:
    """The wild type with every ``model_loop`` gap built, or None if there is none."""
    gaps = gaps_to_model(wild_type, plan)
    if not gaps:
        return None
    protocol = inputs.manifest.get("modelling") or inputs.ruleset.modelling
    placed = place_loops(wild_type, gaps, protocol["placement"])
    built = with_loops(wild_type, gaps, (placed, protocol["marking"]))
    relaxed = _minimize(built, gaps, (protocol, inputs.ruleset))
    structure = relaxed.outcome.structure
    config = config_from_manifest(inputs.manifest, inputs.manifest_sha256)
    raw = [
        f
        for gap in gaps
        for f in gap_findings(structure, gap, (inputs.ruleset, protocol["checks"]))
    ]
    findings = [f.to_dict() for f in apply_context(raw, structure, config)]
    return Modelled(structure, gaps, relaxed, findings, protocol)


def _minimize(built: Structure, gaps: tuple[Gap, ...], context: tuple) -> Relaxed:
    """Pass 1: restrained (bonded pre-stage, trans-omega and L-chirality torsion
    restraints). Pass 2: hydrogens added afresh to the corrected heavy atoms, no torsion
    restraints. Hydrogens placed on a D placement stay on the D side when the pre-stage
    flips CB, and pull the centre back once the restraints go (ADR-0007, amended)."""
    protocol, ruleset = context
    settings = protocol["minimization"]
    shell = loop_shell(built, gaps)
    torsions = torsion_restraints(built, gaps, settings["l_chirality_improper_degree"])
    restrained = relax(built, shell, settings, torsions)
    final = relax(restrained.structure, loop_shell(restrained.structure, gaps), settings)
    criterion = clash_criterion(ruleset)
    return Relaxed(
        MODELLED_DIR,
        _against(built, final),
        shell,
        tuple(c for gap in gaps for _, c in loop_clashes(built, gap, criterion)),
        tuple(c for gap in gaps for _, c in loop_clashes(final.structure, gap, criterion)),
    )


def _against(built: Structure, final: RelaxOutcome) -> RelaxOutcome:
    """``final`` with what moved measured from the placed loop (``built``), not from pass 1."""
    shifts, moved = [], []
    index = final.structure.residue_index
    for residue in built.residues:
        after = {a.name: a.position for a in index[residue.id].heavy_atoms}
        mine = [math.dist(a.position, after[a.name]) for a in residue.heavy_atoms]
        shifts += mine
        if any(d >= MOVED_ANGSTROM for d in mine):
            moved.append(residue.id)
    moving = [d for d in shifts if d >= MOVED_ANGSTROM]
    return replace(
        final,
        moved_residues=tuple(moved),
        moved_atoms=len(moving),
        max_displacement_angstrom=max(shifts, default=0.0),
        rms_displacement_angstrom=math.sqrt(sum(d * d for d in moving) / max(1, len(moving))),
    )


def write_model(modelled: Modelled, context: tuple, out_dir: Path) -> dict:
    """Write wt_modelled/ (if accepted), model_record.json and model_report.md; return
    the record, or raise LoopRejected after writing the record of a rejected model."""
    inputs, wt_record = context
    system = None
    if not modelled.rejected:
        directory = out_dir / MODELLED_DIR
        directory.mkdir(parents=True, exist_ok=True)
        system = {
            "directory": MODELLED_DIR,
            "files": write_system(modelled.structure, directory, SYSTEM_NAME),
            "counts": system_counts(inputs.structure, modelled.structure),
        }
    record = _record(modelled, (inputs, wt_record, system), out_dir)
    if modelled.rejected:
        blocking = [f["title"] for f in modelled.findings if f["effective_severity"] == "blocking"]
        raise LoopRejected(
            f"modelled loop rejected ({'; '.join(blocking)}); see {out_dir / RECORD_FILE}"
        )
    return record


def _record(modelled: Modelled, context: tuple, out_dir: Path) -> dict:
    inputs, wt_record, system = context
    outcome = modelled.relaxed.outcome
    source = (
        f"PDBFixer {installed_version('pdbfixer')} placement, then OpenMM "
        f"{outcome.openmm_version} minimization (modelling protocol)"
    )
    record = {
        "schema_version": wt_record["schema_version"],
        "generated_at": utc_now(),
        "simprep": simprep_provenance(),
        "knowledge_base": wt_record["knowledge_base"],
        "input": inputs.source,
        "manifest_sha256": inputs.manifest_sha256,
        "status": "rejected" if modelled.rejected else "modelled",
        "protocol": modelled.protocol,
        "wild_type": {
            "directory": WILD_TYPE_DIR,
            "prep_record_sha256": sha256_file(out_dir / WILD_TYPE_DIR / PREP_RECORD_FILE),
        },
        "loops": [_loop_entry(gap, modelled.structure, source) for gap in modelled.gaps],
        "minimization": relaxation_entry(modelled.relaxed),
        "system": system,
        "findings": modelled.findings,
        "work_order": _work_order(modelled, wt_record["work_order"]),
    }
    validate(record, "model_record")
    write_json(record, out_dir / RECORD_FILE)
    (out_dir / REPORT_FILE).write_text(render_model_report(record))
    return record


def _loop_entry(gap: Gap, structure: Structure, source: str) -> dict:
    residues = [structure.residue(rid) for rid in gap.residue_ids]
    return {
        "finding_id": gap.finding_id,
        "label": gap.label,
        "flanks": [residue_ref(rid, structure.residue(rid).name) for rid in gap.flanks],
        "source": source,
        "note": NO_EXPERIMENTAL_SUPPORT,
        "residues": [
            {
                **residue_ref(r.id, r.name),
                "label_seq": r.label_seq,
                "atoms": [a.name for a in r.atoms],
            }
            for r in residues
        ],
        "peptide_bonds": [
            {
                "first": residue_ref(b.first.id, b.first.name),
                "second": residue_ref(b.second.id, b.second.name),
                "c_n_angstrom": round(b.c_n_angstrom, DIGITS),
                "ca_ca_angstrom": round(b.ca_ca_angstrom, DIGITS),
                "omega_degree": round(b.omega_degree, 1),
            }
            for b in peptide_bonds(structure, gap)
        ],
        "backbone": [
            {
                **residue_ref(b.residue.id, b.residue.name),
                "phi_degree": round(b.phi_degree, 1),
                "psi_degree": round(b.psi_degree, 1),
                "ca_improper_degree": None
                if b.improper_degree is None
                else round(b.improper_degree, 1),
            }
            for b in backbone(structure, gap)
        ],
    }


def _work_order(modelled: Modelled, prep_items: list[dict]) -> list[dict]:
    """The prepared wild type's work order without the gaps built here, plus protonation
    items: modelled residues have no hydrogens, moved flanks lost theirs."""
    built = {gap.finding_id for gap in modelled.gaps}
    items = [
        item
        for item in prep_items
        if not (item["finding_id"] in built and item["option_id"] in APPLIED_OPTIONS)
    ]
    moved = set(modelled.relaxed.outcome.moved_residues)
    for gap in modelled.gaps:
        items.append(_protonation(gap, gap.residue_ids, "modelled residues have no hydrogens"))
        items += [
            _protonation(gap, (rid,), f"{rid.label()}: flank moved in the loop minimization")
            for rid in gap.flanks
            if rid in moved
        ]
    return items


def _protonation(gap: Gap, residues: tuple, description: str) -> dict:
    return {
        "stage": "protonation",
        "finding_id": gap.finding_id,
        "option_id": APPLIED_OPTIONS[0],
        "description": f"{gap.label}: {description}" if len(residues) > 1 else description,
        "residues": [rid.to_dict() for rid in residues],
    }
