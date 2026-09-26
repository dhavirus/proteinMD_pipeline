# Regression panel

Five structures, committed gzipped (mmCIF). Tests never fetch anything: these files are
the source for every regression test. Expected findings live in
`expected_findings/<ID>.yaml`, derived by `derive_expected.py` from the file's own
annotations (no simprep code) plus hand-written `curated_families` / `curated_checks`, all with
`reviewed: false` until the maintainer reviews them.

RCSB (files.rcsb.org, data.rcsb.org) was blocked by this environment's network policy
during TASK-001, so **no candidate was checked against RCSB metadata online**. The
maintainer uploaded the candidate files; every choice below was checked against the
records inside the committed file (`_exptl`, `_refine`, `_entity`, `_struct_conn`,
`_pdbx_unobs_or_zero_occ_residues`, `_atom_site`, `_citation`), and the candidates were
proposed from memory, not from a search.

| slot | ID | resolution | SHA-256 (file as committed) |
|---|---|---|---|
| IDS instance (required) | 5FQL | 2.30 A | `ac5ffc6fed697e7b2659f6747ad4d0638eb949220cd7ca7c8c3b372ce183f8d5` |
| catalytic Zn metalloenzyme | 3KS3 | 0.90 A | `ee750b2e2cdd2d2c27ce505593719a8850525dea53a4f5c9e863c283e8c876f8` |
| glycoprotein, N-glycans resolved | 1HZH | 2.70 A | `9bd5149e49d247ac6210fdad409682ca2f3b8eeb4fadf597c5aedcdbbddb5391` |
| covalently bound ligand | 6OIM | 1.65 A | `fb1594808cb5ce47fe60260b8b28ac3bb8fb1d10dfd6759ca96bade5d4b67241` |
| altlocs + internal missing loop | 1FO8 | 1.40 A | `030af970cb46b6bb39b348e0d54f29ef9ab4787c2b3b4d13d6ecc21c17815fa3` |

All five: X-RAY DIFFRACTION, one model.

## Candidates considered

| slot | candidates proposed | outcome |
|---|---|---|
| catalytic Zn | 3KS3, 2CBA (CA II), 8TLN (thermolysin) | 3KS3 uploaded and chosen; others not uploaded, not checked |
| glycoprotein | 1HZH, 3AVE, 1OGS/2V3F | 1HZH uploaded and chosen; others not uploaded, not checked |
| covalent ligand | 6OIM, 5P9J | 6OIM uploaded and chosen; 5P9J not uploaded, not checked |
| altlocs + internal gap | no confident proposal; maintainer uploaded 1FO8, 11GN, 1KMZ | **1FO8** chosen: one chain, 11 non-water altloc residues, one internal 13-residue gap (318-330). 11GN (STING, 1.29 A) also qualifies (32 altloc residues, several internal gaps) but is larger and multi-chain; 1KMZ has no altlocs and does not qualify. 11GN and 1KMZ were removed from the panel; they remain in git history. |

## Entries

### 5FQL: human iduronate-2-sulfatase

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

#### Resolved (maintainer, TASK-001 review): 7 vs 8 glycosylation sites

The paper states eight putative N-glycosylation sites with at least one NAG built at
each; the deposited file attaches NAG to seven asparagines (115, 144, 246, 280, 325,
513, 537). Decision: rely on the deposited structure; the expected findings follow the
file.

### 3KS3: human carbonic anhydrase II, 0.9 A

- File citation: "A short, strong hydrogen bond in the active site of human carbonic
  anhydrase II", doi:10.1021/bi902007b (PubMed 20000378). Not read in this session.
- Exercises the zinc rule's open-site criterion: Zn 262 coordinated by His94 NE2, His96
  NE2, His119 ND1 and a water (1.90 A) -> `likely_catalytic`. A second water (HOH 707,
  occupancy 0.29) lies at 2.57 A, inside the 2.8 A cutoff, giving coordination number 5;
  the evidence labels its partial occupancy. 13 non-water altloc residues; glycerol
  (GOL) as unrecognized chemistry; N-terminal residues 1-3 unobserved.

### 1HZH: intact human IgG1 b12, 2.7 A

- File citation: "Crystal structure of a neutralizing human IGG against HIV-1: a
  template for vaccine design", doi:10.1126/science.1061692 (PubMed 11498595). Not read.
- Two complex biantennary N-glycans (branched entities, chains A and B; 9 residues each,
  NAG/BMA/MAN/GAL/FUC), each one `unrecognized` group. Internal gaps in heavy chain K
  (130-136, 236-238) and a C-terminal gap. No altlocs.
- **The file has no ASN-NAG `struct_conn` record**: the glycans are reported without an
  attachment. Measured: NAG 1 C1 lies 2.64 A (chain A glycan) and 2.45 A (chain B glycan)
  from ND2 of ASN 314 in heavy chains H and K, i.e. the Fc glycosylation site, but at a
  non-bonding distance and without annotation. Since TASK-003 the `covalent_contacts`
  family reports both as blocking `n_glycosylation` candidates (recommended `add_link`:
  ASN 314 starts an N-S-T sequon in the file's own entity sequence). They are listed by
  hand in `curated_families` of `expected_findings/1HZH.yaml`, because no annotation
  can give them.

### 6OIM: KRAS G12C with AMG 510 (sotorasib), 1.65 A

- File citation: "The clinical KRAS(G12C) inhibitor AMG 510 drives anti-tumour
  immunity", doi:10.1038/s41586-019-1694-1 (PubMed 31666701). Not read.
- AMG 510 is component **MOV** ("AMG 510 (bound form)"), covalently linked to Cys12 SG
  (`struct_conn` covale1): reported as one blocking `unrecognized` group with the Cys12
  attachment; Cys12 itself stays a standard residue. Mg 301 is octahedral (Ser17 OG,
  GDP O2B, four waters) and is classified `ambiguous`: by maintainer decision a bound
  ligand does not make Mg catalytic, and the file's crystallization conditions
  ("1mM MgCl2, 0.1M MES pH6.5, 30% PEG4000") are shown as evidence. GDP is unrecognized chemistry. Gaps: N-terminal tag (-13 to -1)
  and internal 105-107.

### 1FO8: rabbit N-acetylglucosaminyltransferase I, 1.4 A

- File citation: "X-ray crystal structure of rabbit N-acetylglucosaminyltransferase I:
  catalytic mechanism and a new protein superfamily", doi:10.1093/emboj/19.20.5269
  (PubMed 11032794). Not read.
- One internal 13-residue gap (318-330), 11 non-water altloc residues plus one water.
- Methylmercury (MMC, a heavy-atom derivative) bound to Cys123 and backbone O: a metal
  inside a multi-atom residue is not treated as an ion; it is one `unrecognized` group
  whose `metalc` records appear as attachments.
