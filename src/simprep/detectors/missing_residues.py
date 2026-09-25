"""Unobserved polymer residues from the file's annotations, as terminal or internal runs."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from simprep.config import AuditConfig
from simprep.detectors import register
from simprep.detectors.common import primary_atom
from simprep.findings import (
    Finding,
    Locus,
    atom_ref,
    ev_count,
    ev_distance,
    ev_flag,
    ev_label,
    ev_number,
    ev_source_record,
)
from simprep.rules import Rule, RuleSet
from simprep.structure.geometry import distance
from simprep.structure.model import Atom, Residue, Structure, UnobservedResidue

UNOBSERVED_CATEGORY = "pdbx_unobs_or_zero_occ_residues"
FLANK_ATOM = "CA"


@dataclass(frozen=True)
class Run:
    """Consecutive unobserved residues of one chain and their observed neighbours."""

    chain: str
    members: tuple[UnobservedResidue, ...]
    position: str
    before: Residue | None
    after: Residue | None


@register("missing_residues")
def detect_missing_residues(
    structure: Structure, ruleset: RuleSet, config: AuditConfig
) -> list[Finding]:
    """One finding per run of consecutive unobserved polymer residues."""
    rules = {rule.matcher["position"]: rule for rule in ruleset.family("missing_residues")}
    findings = []
    for run in unobserved_runs(structure):
        if run.position not in rules:
            raise LookupError(f"no missing_residues rule for position {run.position!r}")
        findings.append(_finding(run, rules[run.position]))
    return findings


def unobserved_runs(structure: Structure) -> list[Run]:
    """Split unobserved polymer residues into runs by entity sequence position."""
    by_chain: dict[str, list[UnobservedResidue]] = defaultdict(list)
    for unobserved in structure.unobserved_residues:
        if unobserved.is_polymer:
            by_chain[unobserved.residue.chain].append(unobserved)
    runs = []
    for chain in sorted(by_chain):
        observed = {r.label_seq: r for r in structure.polymer_residues(chain) if r.label_seq}
        members = sorted(by_chain[chain], key=_sequence_position)
        runs += [_run(chain, group, observed) for group in _consecutive(members)]
    return runs


def _sequence_position(unobserved: UnobservedResidue) -> int:
    if unobserved.label_seq is None:
        raise ValueError(
            f"unobserved residue {unobserved.residue.label()} has no label_seq_id; "
            "cannot place it in the entity sequence"
        )
    return unobserved.label_seq


def _consecutive(members: list[UnobservedResidue]) -> list[list[UnobservedResidue]]:
    groups: list[list[UnobservedResidue]] = []
    for member in members:
        if groups and member.label_seq == groups[-1][-1].label_seq + 1:
            groups[-1].append(member)
        else:
            groups.append([member])
    return groups


def _run(chain: str, group: list[UnobservedResidue], observed: dict[int, Residue]) -> Run:
    first, last = group[0].label_seq, group[-1].label_seq
    if not observed:
        position = "whole_chain"
    elif last < min(observed):
        position = "n_terminal"
    elif first > max(observed):
        position = "c_terminal"
    else:
        position = "internal"
    return Run(chain, tuple(group), position, observed.get(first - 1), observed.get(last + 1))


def _finding(run: Run, rule: Rule) -> Finding:
    first, last = run.members[0], run.members[-1]
    span = f"{first.residue.label()}-{last.residue.seq_num}{last.residue.ins_code}"
    evidence = [
        ev_label("position", run.position),
        ev_count("gap_length", len(run.members)),
        ev_label("missing_sequence", "-".join(m.res_name for m in run.members)),
        ev_source_record(
            "unobserved_record",
            UNOBSERVED_CATEGORY,
            {
                "first": first.residue.label(),
                "last": last.residue.label(),
                "first_label_seq": first.label_seq,
                "last_label_seq": last.label_seq,
            },
        ),
    ]
    evidence += _flank_evidence(run, rule.matcher["peptide_rise_angstrom"])
    anchors = tuple(r.id for r in (run.before, run.after) if r is not None)
    return Finding(
        id=f"missing_residues/{span}",
        rule=rule,
        title=f"{rule.title}: {span} ({len(run.members)} residues)",
        locus=Locus(
            first.residue, first.res_name, (), tuple((m.residue, m.res_name) for m in run.members)
        ),
        anchor_residues=anchors,
        evidence=tuple(evidence),
    )


def _flank_evidence(run: Run, rise_angstrom: float) -> list[dict]:
    """Flanking residues and, for internal gaps, whether the gap can be bridged."""
    items = [
        ev_label(f"flank_{side}", f"{r.name}{r.id.label()}")
        for side, r in (("before", run.before), ("after", run.after))
        if r is not None
    ]
    if run.position != "internal" or run.before is None or run.after is None:
        return items
    atoms = [_flank_atom(run.before), _flank_atom(run.after)]
    if None in atoms:
        return items + [ev_label("flank_distance", "not computed: flanking CA missing")]
    d = distance(atoms[0].position, atoms[1].position)
    max_span = rise_angstrom * (len(run.members) + 1)
    return items + [
        ev_distance(
            "flank_ca_distance", d, [atom_ref(run.before, atoms[0]), atom_ref(run.after, atoms[1])]
        ),
        ev_number(
            "max_bridgeable_span",
            max_span,
            "angstrom",
            note=f"{rise_angstrom} A x (gap length + 1)",
        ),
        ev_flag("bridgeable", d <= max_span),
    ]


def _flank_atom(residue: Residue) -> Atom | None:
    candidates = tuple(atom for atom in residue.atoms if atom.name == FLANK_ATOM)
    return primary_atom(candidates) if candidates else None
