from simprep.audit import run_audit
from simprep.detectors.covalent_contacts import compatible, detect_covalent_contacts, in_sequon
from simprep.structure.model import ResidueId
from tests.unit import builders as b

ASN = ResidueId("A", 80)
NAG = ResidueId("B", 1)
SEQUON = ("ASN", "GLY", "THR")


def glycan_site(nd2_c1_angstrom=2.5, sequence=SEQUON, altlocs=("", ""), links=()):
    """ASN (label_seq 1) with ND2 at the origin and a NAG whose C1 lies on the x axis."""
    asn = b.residue("A", 80, "ASN", b.P, [b.atom("ND2", "N", (0, 0, 0), altloc=altlocs[0])], 1)
    nag = b.residue(
        "B",
        1,
        "NAG",
        b.BR,
        [
            b.atom("C1", "C", (nd2_c1_angstrom, 0, 0), altloc=altlocs[1]),
            b.atom("O5", "O", (nd2_c1_angstrom + 1.4, 0, 0)),
        ],
    )
    return b.structure([asn, nag], links, polymer_sequences=(("A", sequence),))


def detect(structure, ruleset):
    return detect_covalent_contacts(structure, ruleset, b.config())


def test_unannotated_glycan_contact_is_reported(ruleset):
    (finding,) = detect(glycan_site(), ruleset)
    assert finding.id == "covalent_contacts/n_glycosylation/B:1-A:80"
    assert finding.evidence_value("pattern") == "n_glycosylation"
    assert finding.evidence_value("contact_distance") == 2.5
    assert finding.evidence_value("sequon") == "ASN-GLY-THR"
    assert finding.evidence_value("in_sequon") is True
    assert finding.anchor_residues == (NAG, ASN)


def test_annotated_link_suppresses_the_pair(ruleset):
    link = b.link("covale1", "covale", [(ASN, "ASN", "ND2"), (NAG, "NAG", "C1")], 1.45)
    assert detect(glycan_site(links=[link]), ruleset) == []


def test_contact_beyond_the_limit_is_ignored(ruleset):
    assert detect(glycan_site(nd2_c1_angstrom=3.2), ruleset) == []


def test_incompatible_altlocs_are_never_compared(ruleset):
    assert detect(glycan_site(altlocs=("A", "B")), ruleset) == []
    assert len(detect(glycan_site(altlocs=("A", "A")), ruleset)) == 1
    assert len(detect(glycan_site(altlocs=("A", "")), ruleset)) == 1


def test_non_pattern_atoms_are_ignored(ruleset):
    site = glycan_site()
    asn, nag = site.residues
    swapped = b.residue("A", 80, "ASN", b.P, [b.atom("OD1", "O", (0, 0, 0))], 1)
    assert (
        detect(b.structure([swapped, nag], polymer_sequences=site.polymer_sequences), ruleset) == []
    )


def contact(findings):
    return next(f for f in findings if f.rule.family == "covalent_contacts")


def test_sequon_drives_the_recommendation(ruleset):
    inside = run_audit(glycan_site(), ruleset, b.config())
    outside = run_audit(glycan_site(sequence=("ASN", "PRO", "THR")), ruleset, b.config())
    assert contact(inside).recommended_option == "add_link"
    assert contact(outside).evidence_value("in_sequon") is False
    assert contact(outside).recommended_option == "expert_review"
    assert contact(inside).effective_severity == "blocking"


def test_missing_sequence_is_reported_not_guessed(ruleset):
    (finding,) = detect(glycan_site(sequence=()), ruleset)
    assert finding.evidence_value("in_sequon") is None
    assert "not determinable" in finding.evidence_value("sequon")


def test_candidate_names_the_glycan_group_finding(ruleset):
    findings = run_audit(glycan_site(), ruleset, b.config())
    assert contact(findings).evidence_value("ligand_group_finding") == "unrecognized/group/B:1"


def test_cysteine_adduct_pattern(ruleset):
    cys = b.residue("A", 12, "CYS", b.P, [b.atom("SG", "S", (0, 0, 0))], 12)
    ligand = b.residue(
        "A", 900, "LIG", b.NP, [b.atom("C7", "C", (1.9, 0, 0)), b.atom("O1", "O", (0.5, 0, 0))]
    )
    (finding,) = detect(b.structure([cys, ligand]), ruleset)
    assert finding.evidence_value("pattern") == "cysteine_adduct"
    assert finding.locus.atom_names == ("C7",)
    assert finding.evidence_value("in_sequon") is None
    assert finding.rule.recommended_option["default"] == "expert_review"


def test_sequon_rules():
    assert in_sequon(("ASN", "GLY", "SER"))
    assert not in_sequon(("ASN", "PRO", "SER"))
    assert not in_sequon(("ASN", "GLY", "ALA"))
    assert in_sequon(("ASN", "GLY", "ALA,THR"))


def test_altloc_compatibility():
    a, b_, blank = (
        b.atom("X", "C", (0, 0, 0), altloc="A"),
        b.atom("X", "C", (0, 0, 0), altloc="B"),
        b.atom("X", "C", (0, 0, 0)),
    )
    assert compatible(a, blank) and compatible(blank, b_) and not compatible(a, b_)
