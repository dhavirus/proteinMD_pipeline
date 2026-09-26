"""Non-standard residues in polymer chains, from residue names and modification records."""

from __future__ import annotations

from simprep.config import AuditConfig
from simprep.detectors import register
from simprep.detectors.common import GENERIC_NONSTANDARD_KIND, standard_polymer_residues
from simprep.findings import (
    Finding,
    ev_b_factor,
    ev_count,
    ev_label,
    ev_occupancy,
    ev_source_record,
    residue_locus,
)
from simprep.rules import Rule, RuleSet
from simprep.structure.model import ModifiedResidue, Residue, ResidueClass, Structure


@register("nonstandard_residues")
def detect_nonstandard_residues(
    structure: Structure, ruleset: RuleSet, config: AuditConfig
) -> list[Finding]:
    """One finding per residue claimed by a component rule or the generic rule."""
    rules = ruleset.family("nonstandard_residues")
    standard = standard_polymer_residues(ruleset)
    annotations = {mod.residue: mod for mod in structure.modified_residues}
    findings = []
    for residue in structure.residues:
        sources = _detection_sources(residue, standard, annotations)
        rule = _matching_rule(residue, rules, (sources, annotations.get(residue.id)))
        if rule is not None:
            findings.append(_finding(residue, rule, (sources, annotations.get(residue.id))))
    return findings


def _detection_sources(
    residue: Residue, standard: frozenset[str], annotations: dict
) -> tuple[str, ...]:
    """Why the residue is non-standard: its name, its modification record, or both."""
    sources = []
    if residue.residue_class is ResidueClass.POLYMER and residue.name not in standard:
        sources.append("residue_name")
    if residue.id in annotations:
        sources.append("modification_annotation")
    return tuple(sources)


def _matching_rule(
    residue: Residue, rules: tuple[Rule, ...], sources_annotation: tuple
) -> Rule | None:
    sources, annotation = sources_annotation
    for rule in rules:
        matcher = rule.matcher
        if matcher["kind"] == "component_id" and (
            residue.name in matcher["comp_ids"]
            and residue.residue_class.value in matcher["residue_classes"]
        ):
            return rule
        if matcher["kind"] == GENERIC_NONSTANDARD_KIND and _generic_match(
            matcher, (sources, annotation)
        ):
            return rule
    return None


def _generic_match(matcher: dict, sources_annotation: tuple) -> bool:
    sources, annotation = sources_annotation
    if not matcher["use_modification_annotations"] or (
        matcher["skip_annotations_on_parent_residue"] and is_attachment_site(annotation)
    ):
        sources = tuple(s for s in sources if s != "modification_annotation")
    return bool(sources)


def is_attachment_site(annotation: ModifiedResidue | None) -> bool:
    """An annotation whose component is its own parent marks a site, not a new residue."""
    return annotation is not None and annotation.res_name == annotation.parent_res_name


def modeled_form(residue: Residue, matcher: dict) -> str:
    """Name of the first modeled form whose atom criteria the residue satisfies."""
    names = {atom.name for atom in residue.atoms}
    for form in matcher.get("modeled_forms", []):
        if (
            residue.name in form["comp_ids"]
            and set(form["required_atoms"]) <= names
            and not set(form.get("absent_atoms", [])) & names
        ):
            return form["form"]
    return "unclassified"


def _finding(
    residue: Residue, rule: Rule, sources_annotation: tuple[tuple, ModifiedResidue | None]
) -> Finding:
    sources, annotation = sources_annotation
    evidence = [
        ev_label("detected_by", ",".join(sources) or "component_id"),
        ev_count("atom_records", len(residue.atoms)),
        ev_label("atom_names", " ".join(sorted({atom.name for atom in residue.atoms}))),
        ev_occupancy("mean_occupancy", _mean([a.occupancy for a in residue.atoms])),
        ev_b_factor("mean_b_factor", _mean([a.b_iso for a in residue.atoms])),
    ]
    if "modeled_forms" in rule.matcher:
        evidence.append(
            ev_label(
                "modeled_form",
                modeled_form(residue, rule.matcher),
                note="from the atoms present in the deposited residue",
            )
        )
    if annotation is not None:
        evidence.append(
            ev_source_record(
                "modification_record",
                "pdbx_struct_mod_residue",
                {
                    "residue": residue.id.label(),
                    "res_name": annotation.res_name,
                    "parent_comp_id": annotation.parent_res_name,
                    "details": annotation.details,
                },
            )
        )
    return Finding(
        id=f"nonstandard_residues/{residue.id.label()}",
        rule=rule,
        title=f"{rule.title}: {residue.name} {residue.id.label()}",
        locus=residue_locus(residue),
        anchor_residues=(residue.id,),
        evidence=tuple(evidence),
        claims=frozenset({residue.id}),
    )


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0
