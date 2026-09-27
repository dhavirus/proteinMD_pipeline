"""Protonation states from pKa estimates (pure; TASK-009, ADR-0009).

Order of rules per titratable residue: a cysteine in a disulfide is CYX; a residue whose
side chain is linked to a metal in struct_conn is not protonated there (Asp/Glu/Cys
deprotonated, His HID when NE2 is bound, HIE when ND1 is); otherwise a decision on its
finding, else the finding's recommended option, else the pKa estimate at the pH.
Findings: residues within the method's distance of a metal ion or non-standard residue
(unreliable estimate), or within the ambiguity window of the pH (ambiguous state).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from simprep.config import AuditConfig
from simprep.findings import Finding, ev_label, ev_number, residue_locus
from simprep.rules import RuleSet
from simprep.severity import apply_context
from simprep.structure.model import Residue, ResidueClass, ResidueId, Structure

FAMILY = "protonation"
# (protonated, deprotonated) variant per residue type; None = neutral His, whose tautomer
# OpenMM chooses by hydrogen bonding. Tyr and Arg have no alternative in the force field.
VARIANTS = {
    "ASP": ("ASH", "ASP"),
    "GLU": ("GLH", "GLU"),
    "LYS": ("LYS", "LYN"),
    "CYS": ("CYS", "CYX"),
    "HIS": ("HIP", None),
}
METAL_BOUND = {"ASP": "ASP", "GLU": "GLU", "CYS": "CYX"}
HIS_BOUND = {"NE2": "HID", "ND1": "HIE"}
OPTION_STATE = {"predicted_state": "pka", "standard_state": "model_pka"}
APPLIED_OPTIONS = (*OPTION_STATE, "specific_state")  # carried out by this stage


class ProtonationError(ValueError):
    """States cannot be assigned as specified; the message says why."""


@dataclass(frozen=True)
class Estimate:
    """A pKa estimate for one titratable residue."""

    residue: ResidueId
    res_name: str
    pka: float
    model_pka: float


@dataclass(frozen=True)
class State:
    """The variant used for a residue (None: neutral His, tautomer by OpenMM) and why."""

    residue: ResidueId
    res_name: str
    variant: str | None
    basis: str
    pka: float | None


@dataclass(frozen=True)
class Context:
    ph: float
    method: dict
    ruleset: RuleSet
    decisions: dict  # finding id -> decision
    system: str
    config: AuditConfig  # regions and thresholds for the severity context


def assign(structure: Structure, estimates: list[Estimate], context: Context) -> tuple:
    """(states, findings) for every titratable residue with an estimate."""
    special = _special(structure)
    titratable = [e for e in estimates if e.res_name in VARIANTS]
    raw = [
        f
        for e in titratable
        if e.residue not in special and (f := _finding(structure, e, context)) is not None
    ]
    findings = {f.anchor_residues[0]: f for f in apply_context(raw, structure, context.config)}
    states = []
    for estimate in titratable:
        if estimate.residue in special:
            variant, basis = special[estimate.residue]
            states.append(State(estimate.residue, estimate.res_name, variant, basis, estimate.pka))
        else:
            states.append(_state(estimate, findings.get(estimate.residue), context))
    return states, list(findings.values())


def variant_at(estimate: Estimate, value: float, ph: float) -> str | None:
    protonated, deprotonated = VARIANTS[estimate.res_name]
    return protonated if value > ph else deprotonated


def _special(structure: Structure) -> dict:
    """Residue id -> (variant, basis) for disulfide cysteines and metal ligands."""
    special = {}
    for link in structure.links:
        for mine, other in ((link.partner1, link.partner2), (link.partner2, link.partner1)):
            if link.conn_type == "disulf" and mine.res_name == "CYS":
                special[mine.residue] = ("CYX", f"disulfide {link.conn_id}")
            if link.conn_type != "metalc":
                continue
            variant = METAL_BOUND.get(mine.res_name)
            if mine.res_name == "HIS":
                variant = HIS_BOUND.get(mine.atom_name)
            if variant is not None:
                where = f"{other.res_name} {other.residue.label()}"
                special[mine.residue] = (variant, f"metal rule: {mine.atom_name} bound to {where}")
    return special


def _state(estimate: Estimate, finding: Finding | None, context: Context) -> State:
    if finding is None:
        variant = variant_at(estimate, estimate.pka, context.ph)
        return State(estimate.residue, estimate.res_name, variant, "pka", estimate.pka)
    decision = context.decisions.get(finding.id)
    option = decision["option_id"] if decision else finding.recommended_option
    basis = f"{option} ({'decided' if decision else 'recommended, undecided'})"
    if option == "specific_state":
        variant = _specific(estimate, (decision or {}).get("parameters", {}), finding.id)
    elif option in OPTION_STATE:
        value = getattr(estimate, OPTION_STATE[option])
        variant = variant_at(estimate, value, context.ph)
    else:
        raise ProtonationError(f"{finding.id}: {option} is not a final protonation decision")
    return State(estimate.residue, estimate.res_name, variant, basis, estimate.pka)


def _specific(estimate: Estimate, parameters: dict, finding_id: str) -> str | None:
    allowed = {v for v in VARIANTS[estimate.res_name] if v} | (
        {"HID", "HIE"} if estimate.res_name == "HIS" else set()
    )
    state = parameters.get("state")
    if state not in allowed:
        raise ProtonationError(
            f"{finding_id}: specific_state needs parameters.state, one of {sorted(allowed)}"
            f" (got {state!r})"
        )
    return state


def _finding(structure: Structure, estimate: Estimate, context: Context) -> Finding | None:
    rules = {rule.matcher["kind"]: rule for rule in context.ruleset.family(FAMILY)}
    residue = structure.residue(estimate.residue)
    limit = context.method["unreliable_within_angstrom"]
    near = _nearest_blind(structure, residue, (limit, _blind_names(context.ruleset)))
    if near is not None:
        rule, extra = rules["unreliable_pka"], [ev_label("near", near[0]), _distance(near[1])]
    elif abs(estimate.pka - context.ph) < context.method["ambiguity_window_ph"]:
        rule, extra = rules["ambiguous_pka"], []
    else:
        return None
    return Finding(
        id=f"{FAMILY}/{context.system}/{estimate.residue.label()}",
        rule=rule,
        title=f"{rule.title}: {residue.name} {estimate.residue.label()} ({context.system})",
        locus=residue_locus(residue),
        anchor_residues=(estimate.residue,),
        evidence=(
            ev_number("pka", estimate.pka, "pH"),
            ev_number("model_pka", estimate.model_pka, "pH"),
            ev_number("ph", context.ph, "pH"),
            *extra,
        ),
    )


def _distance(value: float) -> dict:
    return ev_number("distance_to_ignored_chemistry", value, "angstrom")


def _blind_names(ruleset: RuleSet) -> tuple[frozenset[str], frozenset[str]]:
    """(metal elements of the metals rules, standard residues of the generic
    non-standard rule): what the knowledge base calls a metal and a standard residue."""
    metals = frozenset(e for r in ruleset.family("metals") for e in r.matcher.get("elements", ()))
    generic = next(r for r in ruleset.rules if r.rule_id == "nonstandard_residues.generic")
    return metals, frozenset(generic.matcher["standard_residues"])


def _nearest_blind(structure: Structure, residue: Residue, where: tuple) -> tuple | None:
    """(label, distance) of the closest metal ion or non-standard residue within the limit;
    ``where`` is (limit A, (metal elements, standard residues))."""
    limit, names = where
    best = None
    for other in structure.residues:
        if not _blind(other, names) or other.id == residue.id:
            continue
        distance = min(
            math.dist(a.position, b.position)
            for a in residue.heavy_atoms
            for b in other.heavy_atoms
        )
        if distance <= limit and (best is None or distance < best[1]):
            best = (f"{other.name} {other.id.label()}", distance)
    return best


def _blind(residue: Residue, names: tuple) -> bool:
    """Chemistry PROPKA does not include: metal ions and non-standard polymer residues."""
    metals, standard = names
    if residue.residue_class is ResidueClass.POLYMER:
        return residue.name not in standard
    atoms = residue.heavy_atoms
    return (
        residue.residue_class is ResidueClass.NONPOLYMER
        and len(atoms) == 1
        and atoms[0].element in metals
    )
