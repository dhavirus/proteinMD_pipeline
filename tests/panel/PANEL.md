# Regression panel

Five structures, committed gzipped (mmCIF). Tests never fetch anything: these files are
the source for every regression test. Expected findings live in
`expected_findings/<ID>.yaml`, derived by `derive_expected.py` from the file's own
annotations (no simprep code) plus hand-written `publication_checks`, all with
`reviewed: false` until the maintainer reviews them.

RCSB (files.rcsb.org, data.rcsb.org) was blocked by this environment's network policy
during TASK-001, so no candidate could be checked against RCSB metadata online. Files
were uploaded by the maintainer; each choice below is verified against the records in
the committed file itself, and that is stated per entry.

| slot | ID | status |
|---|---|---|
| IDS instance (required) | 5FQL | committed; verified against file records and the primary publication |
| catalytic Zn metalloenzyme | - | **open**: awaiting file |
| glycoprotein, N-glycans resolved | - | **open**: awaiting file (5FQL itself carries 7 resolved N-glycans) |
| covalently bound ligand | - | **open**: awaiting file |
| altlocs + internal missing loop | - | **open**: awaiting file (5FQL has both: 2 altloc residues, gap 444-453) |

## 5FQL: human iduronate-2-sulfatase

- Source: `5FQL.cif.gz`, uploaded by the maintainer (SHA-256 recorded by every audit).
- File records: `_exptl.method` X-RAY DIFFRACTION, `_refine.ls_d_res_high` 2.300,
  one model, one polymer chain A (auth 34-550 observed), entities: polymer, two branched
  oligosaccharides (3 + 2 copies), CA, NAG (2), CL (3), water (213).
- Publication: Demydchuk M et al. Insights into Hunter syndrome from the structure of
  iduronate-2-sulfatase. Nat Commun 2017;8:15786, doi:10.1038/ncomms15786 (full text read
  via PubMed Central PMC5472762).
- What it exercises:
  - metals: Ca 1551, 6-coordinate (D45, D46, D334, H335, two sulfate O of residue 84);
    the file's `struct_conn` also lists the sulfate S-Ca contact (2.897 A) as `metalc`,
    which is not a donor and is excluded from the shell by the rule's donor elements.
  - nonstandard_residues: residue 84 is component **ALS**, "(3S)-3-(sulfooxy)-L-serine"
    (file `chem_comp`; `pdbx_struct_mod_residue` parent ALA), i.e. the FGly sulfate ester
    the paper calls "FGS84". The formylglycine rule matches ALS; FGS is not a component
    used in this file.
  - `pdbx_struct_mod_residue` also lists seven ASN residues with parent ASN and details
    "GLYCOSYLATION SITE". These are attachment sites, not non-standard residues: the
    generic rule skips annotations whose component equals its parent, and the records
    appear as evidence on the attached glycan findings.
  - altlocs: ARG 294 and ARG 421 (A/B).
  - missing_residues: 26-33 (N-terminal, 8) and 444-453 (internal, 10; the loop the
    paper could not model).
  - unrecognized: 5 branched glycans (NAG-NAG[-FUC]) and 2 single NAG, each attached to
    an ASN; 3 chloride ions. All blocking, as specified for TASK-001.

### Open question

The paper states eight putative N-glycosylation sites with at least one NAG built at
each; the deposited file attaches NAG to seven asparagines (115, 144, 246, 280, 325,
513, 537). The expected findings follow the file.
