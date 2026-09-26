"""Which gaps the modelling stage builds, and the structure with the loops in (pure;
TASK-007, ADR-0007).

A gap is modelled when its missing_residues finding is decided ``model_loop`` (an apply
option carried out by this stage). Only internal gaps are built; the sequence comes from
the unobserved-residue annotation and must agree with the entity sequence.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from simprep.detectors.missing_residues import Run, unobserved_runs
from simprep.prep.plan import PrepPlan
from simprep.relax.openmm_run import TorsionRestraint
from simprep.relax.shell import Shell
from simprep.structure.model import Atom, Residue, ResidueClass, ResidueId, Structure

MODEL_LOOP = "model_loop"
APPLIED_OPTIONS = (MODEL_LOOP,)  # apply options carried out by the modelling stage
INTERNAL = "internal"
GLYCINE = "GLY"
TRANS_OMEGA_DEGREE = 180.0


class ModelError(ValueError):
    """The decided gaps cannot be modelled as specified; the message says why."""


@dataclass(frozen=True)
class Gap:
    """An internal run of unobserved residues decided ``model_loop``."""

    finding_id: str
    run: Run

    @property
    def label(self) -> str:
        first, last = self.run.members[0].residue, self.run.members[-1].residue
        return f"{first.label()}-{last.seq_num}{last.ins_code}"

    @property
    def residue_ids(self) -> tuple[ResidueId, ...]:
        return tuple(member.residue for member in self.run.members)

    @property
    def flanks(self) -> tuple[ResidueId, ResidueId]:
        return self.run.before.id, self.run.after.id


def gaps_to_model(structure: Structure, plan: PrepPlan) -> tuple[Gap, ...]:
    """The gaps of ``structure`` (a prepared system) that ``plan`` decides ``model_loop``."""
    runs = {tuple(m.residue for m in run.members): run for run in unobserved_runs(structure)}
    gaps = []
    for action in plan.actions:
        if action.option_id != MODEL_LOOP or action.prep_action != "apply":
            continue
        run = runs.get(action.residues)
        if run is None or run.position != INTERNAL or None in (run.before, run.after):
            raise ModelError(
                f"{action.finding_id}: model_loop builds internal gaps between two observed "
                "residues of the prepared system only; this run is not one"
            )
        _check_sequence(structure, action.finding_id, run)
        gaps.append(Gap(action.finding_id, run))
    return tuple(gaps)


def _check_sequence(structure: Structure, finding_id: str, run: Run) -> None:
    sequence = structure.sequence(run.chain)
    for member in run.members:
        expected = sequence[member.label_seq - 1] if member.label_seq <= len(sequence) else None
        if expected != member.res_name:
            raise ModelError(
                f"{finding_id}: {member.residue.label()} is {member.res_name} in the "
                f"unobserved-residue annotation but {expected!r} in the entity sequence"
            )


def with_loops(structure: Structure, gaps: tuple[Gap, ...], built: tuple) -> Structure:
    """``structure`` with each gap's residues inserted after its first flank, their atoms
    marked (``built`` is (placed atoms per residue id, protocol ``marking``)), and the
    residues removed from the unobserved-residue annotation."""
    placed, marking = built
    new = {gap.run.before.id: new_residues(gap, placed, marking) for gap in gaps}
    residues = []
    for residue in structure.residues:
        residues += [residue, *new.get(residue.id, [])]
    modelled = {rid for gap in gaps for rid in gap.residue_ids}
    unobserved = tuple(u for u in structure.unobserved_residues if u.residue not in modelled)
    return replace(structure, residues=tuple(residues), unobserved_residues=unobserved)


def new_residues(gap: Gap, placed: dict, marking: dict) -> list[Residue]:
    return [
        Residue(
            member.residue,
            member.res_name,
            ResidueClass.POLYMER,
            tuple(
                replace(a, occupancy=marking["occupancy"], b_iso=marking["b_iso_angstrom2"])
                for a in placed[member.residue]
            ),
            label_seq=member.label_seq,
        )
        for member in gap.run.members
    ]


def segment(structure: Structure, gap: Gap) -> list[Residue]:
    """The flank before, the modelled residues and the flank after, in chain order."""
    index = structure.residue_index
    return [index[rid] for rid in (gap.flanks[0], *gap.residue_ids, gap.flanks[1])]


def loop_shell(structure: Structure, gaps: tuple[Gap, ...]) -> Shell:
    """Loop heavy atoms move freely; flank heavy atoms move under the positional restraint."""
    loops = frozenset(rid for gap in gaps for rid in gap.residue_ids)
    flanks = frozenset(rid for gap in gaps for rid in gap.flanks)
    mobile = frozenset(
        (rid, atom.name)
        for rid in loops | flanks
        for atom in structure.residue(rid).atoms
        if atom.is_heavy
    )
    return Shell(loops, loops | flanks, mobile, frozenset(k for k in mobile if k[0] in loops))


def torsion_restraints(structure: Structure, gaps: tuple, target_degree: float) -> tuple:
    """Trans omega for every peptide bond from flank to flank; ``target_degree`` for the
    CA-N-C-CB improper of every non-glycine modelled residue (L)."""
    restraints = []
    for gap in gaps:
        residues = segment(structure, gap)
        for first, second in zip(residues, residues[1:], strict=False):
            atoms = ((first, "CA"), (first, "C"), (second, "N"), (second, "CA"))
            restraints.append(TorsionRestraint(_keys(atoms), TRANS_OMEGA_DEGREE))
        for residue in residues[1:-1]:
            if residue.name != GLYCINE:
                atoms = tuple((residue, name) for name in ("CA", "N", "C", "CB"))
                restraints.append(TorsionRestraint(_keys(atoms), target_degree))
    return tuple(restraints)


def _keys(atoms: tuple) -> tuple:
    return tuple((residue.id, name) for residue, name in atoms)


def atom(residue: Residue, name: str) -> Atom:
    """The heavy atom ``name`` of ``residue``; ModelError if absent."""
    found = next((a for a in residue.atoms if a.name == name), None)
    if found is None:
        raise ModelError(f"{residue.name} {residue.id.label()} has no atom {name}")
    return found
