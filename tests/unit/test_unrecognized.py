from simprep.audit import run_audit
from simprep.structure.model import ModifiedResidue, ResidueId
from tests.unit import builders as b

ASN = ResidueId("A", 80)
NAG1 = ResidueId("A", 1001)
NAG2 = ResidueId("A", 1002)


def glycosylated():
    asn = b.residue("A", 80, "ASN", b.P, [b.atom("ND2", "N", (0, 0, 0))], 80)
    nag1 = b.residue(
        "A", 1001, "NAG", b.BR, [b.atom("C1", "C", (1.45, 0, 0)), b.atom("O4", "O", (3, 0, 0))]
    )
    nag2 = b.residue("A", 1002, "NAG", b.BR, [b.atom("C1", "C", (4.4, 0, 0))])
    so4 = b.residue("A", 1100, "SO4", b.NP, [b.atom("S", "S", (20, 0, 0))])
    links = [
        b.link("covale1", "covale", [(ASN, "ASN", "ND2"), (NAG1, "NAG", "C1")], 1.45),
        b.link("covale2", "covale", [(NAG1, "NAG", "O4"), (NAG2, "NAG", "C1")], 1.44),
    ]
    site = ModifiedResidue(ASN, "ASN", "ASN", "GLYCOSYLATION SITE")
    return b.structure([asn, nag1, nag2, so4], links, modified_residues=(site,))


def test_glycan_tree_is_one_blocking_group_with_attachment(ruleset):
    findings = {f.id: f for f in run_audit(glycosylated(), ruleset, b.config())}
    assert set(findings) == {"unrecognized/group/A:1001", "unrecognized/group/A:1100"}
    glycan = findings["unrecognized/group/A:1001"]
    assert glycan.effective_severity == "blocking"
    assert glycan.recommended_option == "expert_review"
    assert glycan.evidence_value("components") == "NAG-NAG"
    attachments = [e for e in glycan.evidence if e["key"] == "attachment"]
    assert [a["fields"]["id"] for a in attachments] == ["covale1"]
    sites = [e for e in glycan.evidence if e["key"] == "attachment_site_annotation"]
    assert sites[0]["fields"]["details"] == "GLYCOSYLATION SITE"


def cys_pair(link_type="disulf", atoms=("SG", "SG")):
    c1 = b.residue(
        "A", 10, "CYS", b.P, [b.atom("SG", "S", (0, 0, 0)), b.atom("CB", "C", (-1.5, 0, 0))], 10
    )
    c2 = b.residue(
        "A", 50, "CYS", b.P, [b.atom("SG", "S", (2.05, 0, 0)), b.atom("CB", "C", (3.5, 0, 0))], 50
    )
    connection = b.link("x1", link_type, [(c1.id, "CYS", atoms[0]), (c2.id, "CYS", atoms[1])])
    return b.structure([c1, c2], [connection])


def test_disulfide_is_recognized(ruleset):
    assert run_audit(cys_pair(), ruleset, b.config()) == []


def test_odd_covalent_link_between_standard_residues_is_blocking(ruleset):
    (finding,) = run_audit(cys_pair("covale", ("SG", "CB")), ruleset, b.config())
    assert finding.id == "unrecognized/link/x1"
    assert finding.effective_severity == "blocking"


def test_peptide_link_of_claimed_nonstandard_residue_is_recognized(ruleset):
    ala = b.residue("A", 4, "ALA", b.P, [b.atom("C", "C", (0, 0, 0))], 4)
    mse = b.residue("A", 5, "MSE", b.P, [b.atom("N", "N", (1.33, 0, 0))], 5)
    peptide = b.link("covale9", "covale", [(ala.id, "ALA", "C"), (mse.id, "MSE", "N")], 1.33)
    findings = run_audit(b.structure([ala, mse], [peptide]), ruleset, b.config())
    assert [f.rule.family for f in findings] == ["nonstandard_residues"]


def test_ligand_c_n_link_is_not_mistaken_for_peptide(ruleset):
    ala = b.residue("A", 4, "ALA", b.P, [b.atom("C", "C", (0, 0, 0))], 4)
    lig = b.residue("A", 900, "LIG", b.NP, [b.atom("N", "N", (1.4, 0, 0))])
    bond = b.link("covale3", "covale", [(ala.id, "ALA", "C"), (lig.id, "LIG", "N")], 1.4)
    (finding,) = run_audit(b.structure([ala, lig], [bond]), ruleset, b.config())
    assert finding.id == "unrecognized/group/A:900"
    assert finding.evidence_value("components") == "LIG"


def test_waters_and_standard_residues_need_no_claim(ruleset):
    site = b.structure(
        [
            b.residue("A", 1, "GLY", b.P, [b.atom("CA", "C", (0, 0, 0))], 1),
            b.water("A", 1000, (5, 5, 5)),
        ]
    )
    assert run_audit(site, ruleset, b.config()) == []
