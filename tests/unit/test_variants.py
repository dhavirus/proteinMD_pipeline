"""Variant validation, application and rotamer decisions on hand-built structures
(TASK-005)."""

from dataclasses import replace

import pytest

from simprep.prep.plan import Action, PrepPlan
from simprep.structure.model import Link, LinkPartner, ResidueId, SequenceReference
from simprep.variants.apply import apply_mutations, mutated_residue
from simprep.variants.build import Candidate
from simprep.variants.findings import Evaluation
from simprep.variants.run import SiteResult, VariantDecisionError, chosen_candidates
from simprep.variants.validate import VariantError, resolve_sites, uniprot_position
from tests.unit import builders as b

SITE = ResidueId("A", 2)


def arg(num=2, with_cb=True):
    atoms = [
        b.atom("N", "N", (1.46, 0, 0)),
        b.atom("CA", "C", (0, 0, 0)),
        b.atom("C", "C", (-0.55, 1.42, 0)),
        b.atom("O", "O", (-0.2, 2.4, 0.5)),
        b.atom("CB", "C", (-0.53, -0.77, -1.2)),
        b.atom("CG", "C", (0.2, -1.9, -1.8)),
        b.atom("NH1", "N", (1.0, -3.0, -2.5)),
        b.atom("HA", "H", (0.3, 0.4, 0.9)),
    ]
    atoms = [a for a in atoms if with_cb or a.name != "CB"]
    return b.residue("A", num, "ARG", b.P, atoms, label_seq=num)


def wild_type(residue=None, links=()):
    gly = b.residue("A", 1, "GLY", b.P, [b.atom("CA", "C", (5, 0, 0))], label_seq=1)
    return b.structure(
        [gly, residue or arg()],
        links,
        polymer_sequences=(("A", ("GLY", "ARG")),),
        sequence_references=(SequenceReference("A", "UNP", "P00001", 1, 2, 11),),
    )


def plan(actions=()):
    return PrepPlan(tuple(actions), (), frozenset(), {}, frozenset(), (), (), {})


def variant(to="GLN", source="ARG", num=2, name="R2Q"):
    mutation = {"chain": "A", "seq_num": num, "ins_code": "", "from": source, "to": to}
    return {"name": name, "mutations": [mutation], "rationale": "t", "references": []}


def resolve(variants, structure=None, actions=(), side_chains=None):
    structure = structure or wild_type()
    return resolve_sites(
        variants, (structure, structure, plan(actions)), side_chains or {"GLN": {}, "TRP": {}}
    )


def test_valid_mutation_resolves_with_uniprot_position():
    (site,) = resolve([variant()])
    assert site.label == "A:2 ARG>GLN"
    assert site.uniprot == {"accession": "P00001", "position": 12}


def test_every_problem_is_listed_at_once():
    variants = [variant(source="LYS"), variant(to="ALA", name="R2A"), variant(num=9, name="X9Q")]
    with pytest.raises(VariantError) as error:
        resolve(variants)
    message = str(error.value)
    assert "R2Q: A:2 LYS>GLN: the structure has ARG there" in message
    assert "the entity sequence has ARG there" in message
    assert "R2A: A:2 ARG>ALA: no side-chain data for ALA" in message
    assert "X9Q: A:9 ARG>GLN: residue not in the prepared system" in message


def test_residue_without_cb_or_covered_by_a_deferred_finding_is_refused():
    deferred = Action(
        "metals/A:9", "metals.x", "restrained", "defer", "parameterization", (SITE,), "note"
    )
    with pytest.raises(VariantError, match="missing CB") as error:
        resolve([variant()], wild_type(arg(with_cb=False)), actions=[deferred])
    assert "covered by metals/A:9 (defer)" in str(error.value)


def test_duplicate_names_and_positions_are_refused():
    twice = variant()
    twice["mutations"] = twice["mutations"] * 2
    with pytest.raises(VariantError, match="mutated twice"):
        resolve([twice])
    with pytest.raises(VariantError, match="variant name used twice"):
        resolve([variant(), variant()])


def test_uniprot_position_is_none_outside_the_mapped_range():
    assert uniprot_position(wild_type(), ResidueId("A", 7)) is None


def test_mutated_residue_keeps_backbone_and_cb_and_drops_hydrogens():
    side_chain = (b.atom("CG", "C", (0.1, -1.9, -1.8)), b.atom("CD", "C", (0.9, -3.0, -2.4)))
    residue = mutated_residue(arg(), "GLN", side_chain)
    assert residue.name == "GLN"
    assert [a.name for a in residue.atoms] == ["N", "CA", "C", "O", "CB", "CG", "CD"]


def link(atom_name):
    return Link(
        f"l_{atom_name}",
        "covale",
        LinkPartner(SITE, "ARG", atom_name),
        LinkPartner(ResidueId("A", 1), "GLY", "CA"),
    )


def test_apply_updates_sequence_and_drops_links_to_removed_atoms():
    structure = wild_type(links=(link("CG"), link("C")))
    new = mutated_residue(arg(), "GLN", (b.atom("CG", "C", (0.1, -1.9, -1.8)),))
    variant_structure = apply_mutations(structure, (new,))
    assert dict(variant_structure.polymer_sequences)["A"] == ("GLY", "GLN")
    assert [k.conn_id for k in variant_structure.links] == ["l_C"]
    assert variant_structure.residues[0] == structure.residues[0]


def evaluation(rotamer, clash_count):
    return Evaluation(Candidate(rotamer, (0.0,), 10.0), (), (), (object(),) * clash_count)


def results():
    return {"variant_build/R2Q/A:2": SiteResult(None, [evaluation("a", 0), evaluation("b", 2)])}


def manifest(option_id, rotamer=None):
    decision = {"finding_id": "variant_build/R2Q/A:2", "option_id": option_id}
    if rotamer:
        decision["parameters"] = {"rotamer": rotamer}
    return {"decisions": [decision]}


@pytest.mark.parametrize("rotamer", ["a", "b"])
def test_a_decided_rotamer_is_chosen_even_when_it_clashes(rotamer):
    """Clashes only flag (maintainer, TASK-005 review); relaxation follows (TASK-006)."""
    chosen = chosen_candidates(manifest("choose_rotamer", rotamer), results())
    assert chosen["variant_build/R2Q/A:2"].candidate.rotamer_id == rotamer


@pytest.mark.parametrize(
    ("document", "message"),
    [
        ({"decisions": []}, "undecided"),
        (manifest("expert_review"), "expert_review is not a final decision"),
        (manifest("choose_rotamer", "zz"), "parameters.rotamer must be one of a, b"),
    ],
)
def test_missing_open_or_unknown_choices_are_refused(document, message):
    with pytest.raises(VariantDecisionError, match=message):
        chosen_candidates(document, results())


def test_sequence_reference_maps_positions():
    reference = SequenceReference("H", "UNP", "P0DOX5", 111, 478, 117)
    assert reference.db_position(111) == 117
    assert reference.db_position(110) is None
    assert replace(reference, db_begin=1).db_position(478) == 368
