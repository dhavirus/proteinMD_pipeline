"""Generate the OpenMM force-field XML for DDZ (gem-diol formylglycine; TASK-010, ADR-0010).

Run once with AmberTools (antechamber, parmchk2, tleap, ParmEd) and RDKit available, e.g.
in a conda env: ``conda create -p ENV -c conda-forge ambertools=26.0 rdkit``; then
``ENV/bin/python tools/parameterize_ddz.py OPENMM_DATA_DIR``, where OPENMM_DATA_DIR holds
``amber14/protein.ff14SB.xml`` (e.g. ``python -c "import openmm.app, os;
print(os.path.join(os.path.dirname(openmm.app.__file__), 'data'))"`` in simprep's env).
Writes ``knowledge/forcefield/ddz.xml`` and ``knowledge/forcefield/ddz_provenance.yaml``.
Not run in CI; the outputs are committed data.

Recipe (TASK-010 decisions 2-3): AM1-BCC charges and GAFF2 types on ACE-DDZ-NME; backbone
atoms (N H CA HA C O) take ff14SB serine types and charges; side-chain charges are
shifted evenly to a neutral residue; every bonded term that touches a side-chain atom
takes GAFF2 parameters (parmchk2 for missing ones); side-chain atoms get unique types, so
those terms apply to DDZ only.
"""

from __future__ import annotations

import math
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import parmed
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem

SMILES = "CC(=O)N[C@@H](C(O)O)C(=O)NC"  # ACE-DDZ-NME; DDZ from the CCD canonical SMILES
SEED = 20260927
# Heavy atoms in SMILES order -> names; hydrogens are named by their heavy atom below.
HEAVY = ["CH3", "C1", "O1", "N", "CA", "CB", "OG1", "OG2", "C", "O", "N2", "CH32"]
HYDROGEN_NAMES = {
    "N": ["H"],
    "CA": ["HA"],
    "CB": ["HB"],
    "OG1": ["HG1"],
    "OG2": ["HG2"],
    "CH3": ["HH31", "HH32", "HH33"],
    "N2": ["HN2"],
    "CH32": ["HH34", "HH35", "HH36"],
}
BACKBONE = {
    "N": "protein-N",
    "H": "protein-H",
    "CA": "protein-CX",
    "HA": "protein-H1",
    "C": "protein-C",
    "O": "protein-O",
}
# Cap atoms stand in for the neighbouring residues' atoms in the chain.
CAP_TYPES = {"C1": "protein-C", "O1": "protein-O", "N2": "protein-N", "HN2": "protein-H"}
SIDE_CHAIN = {
    "CB": "ddz-CB",
    "HB": "ddz-HB",
    "OG1": "ddz-OG",
    "OG2": "ddz-OG",
    "HG1": "ddz-HG",
    "HG2": "ddz-HG",
}
DDZ_ATOMS = [*BACKBONE, *SIDE_CHAIN]
KCAL_TO_KJ = 4.184
ANGSTROM_TO_NM = 0.1
ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "knowledge" / "forcefield"


def run(*command: str, cwd: Path) -> None:
    subprocess.run(command, cwd=cwd, check=True, capture_output=True, text=True)


def fragment(work: Path) -> Path:
    """ACE-DDZ-NME with hydrogens, embedded with a fixed seed, as a named PDB."""
    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))
    if AllChem.EmbedMolecule(mol, randomSeed=SEED) != 0:
        raise RuntimeError("RDKit could not embed the fragment")
    AllChem.MMFFOptimizeMolecule(mol)
    names, counters = {}, {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() > 1:
            names[atom.GetIdx()] = HEAVY[atom.GetIdx()]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            parent = names[atom.GetNeighbors()[0].GetIdx()]
            index = counters.get(parent, 0)
            names[atom.GetIdx()] = HYDROGEN_NAMES[parent][index]
            counters[parent] = index + 1
    position = mol.GetConformer().GetAtomPosition
    lines = [
        f"HETATM{i + 1:5d} {names[i]:<4s} DDZ A   1    {position(i).x:8.3f}{position(i).y:8.3f}"
        f"{position(i).z:8.3f}  1.00  0.00          {a.GetSymbol():>2s}"
        for i, a in enumerate(mol.GetAtoms())
    ]
    path = work / "fragment.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")
    check_l(mol, names)
    return path


def check_l(mol, names) -> None:
    """Positive CA-N-C-CB improper = L (the convention of simprep's modelling checks)."""
    index = {name: i for i, name in names.items()}
    conf = mol.GetConformer()
    improper = Chem.rdMolTransforms.GetDihedralDeg(
        conf, *(index[n] for n in ("CA", "N", "C", "CB"))
    )
    if improper <= 0:
        raise RuntimeError(f"fragment is not L (CA-N-C-CB improper {improper:.1f})")


def amber_parameters(pdb: Path, work: Path) -> tuple[parmed.Structure, dict]:
    run(
        "antechamber",
        "-i",
        pdb.name,
        "-fi",
        "pdb",
        "-o",
        "fragment.mol2",
        "-fo",
        "mol2",
        "-c",
        "bcc",
        "-nc",
        "0",
        "-at",
        "gaff2",
        "-rn",
        "DDZ",
        "-s",
        "0",
        cwd=work,
    )
    run(
        "parmchk2",
        "-i",
        "fragment.mol2",
        "-f",
        "mol2",
        "-o",
        "fragment.frcmod",
        "-s",
        "gaff2",
        cwd=work,
    )
    (work / "leap.in").write_text(
        "source leaprc.gaff2\nloadamberparams fragment.frcmod\nm = loadmol2 fragment.mol2\n"
        "saveamberparm m fragment.prmtop fragment.inpcrd\nquit\n"
    )
    run("tleap", "-f", "leap.in", cwd=work)
    structure = parmed.load_file(str(work / "fragment.prmtop"), str(work / "fragment.inpcrd"))
    charges = {a.name: a.charge for a in structure.atoms}
    return structure, charges


def openmm_type(name: str) -> str | None:
    return SIDE_CHAIN.get(name) or BACKBONE.get(name) or CAP_TYPES.get(name)


def touches_side_chain(atoms) -> bool:
    return any(a.name in SIDE_CHAIN for a in atoms)


def terms(structure: parmed.Structure) -> dict:
    """Bonded terms touching a side-chain atom, keyed by OpenMM type tuples."""
    found = {"bonds": {}, "angles": {}, "torsions": {}}
    for b in structure.bonds:
        if touches_side_chain((b.atom1, b.atom2)):
            key = tuple(openmm_type(a.name) for a in (b.atom1, b.atom2))
            add(
                found["bonds"],
                key,
                (b.type.req * ANGSTROM_TO_NM, 2 * b.type.k * KCAL_TO_KJ / ANGSTROM_TO_NM**2),
            )
    for a in structure.angles:
        if touches_side_chain((a.atom1, a.atom2, a.atom3)):
            key = tuple(openmm_type(x.name) for x in (a.atom1, a.atom2, a.atom3))
            add(found["angles"], key, (math.radians(a.type.theteq), 2 * a.type.k * KCAL_TO_KJ))
    for d in structure.dihedrals:
        atoms = (d.atom1, d.atom2, d.atom3, d.atom4)
        if d.improper or not touches_side_chain(atoms):
            continue
        key = tuple(openmm_type(x.name) for x in atoms)
        if None in key:
            continue  # reaches a cap methyl: that atom is not in the chain
        found["torsions"].setdefault(key, set()).add(
            (d.type.per, round(math.radians(d.type.phase), 6), round(d.type.phi_k * KCAL_TO_KJ, 6))
        )
    return found


def add(table: dict, key: tuple, value: tuple) -> None:
    """One entry per type tuple (either direction); the same key must carry the same value."""
    canonical = min(key, key[::-1])
    rounded = tuple(round(v, 6) for v in value)
    if canonical in table and table[canonical] != rounded:
        raise RuntimeError(f"{canonical}: two different parameters {table[canonical]} {rounded}")
    table[canonical] = rounded


def serine(openmm_data: Path) -> dict:
    """ff14SB serine backbone charges from OpenMM's protein.ff14SB.xml."""
    tree = ET.parse(openmm_data / "amber14" / "protein.ff14SB.xml")
    residue = next(r for r in tree.iter("Residue") if r.get("name") == "SER")
    return {
        a.get("name"): float(a.get("charge"))
        for a in residue.iter("Atom")
        if a.get("name") in BACKBONE
    }


def ddz_charges(am1bcc: dict, backbone: dict) -> tuple[dict, float]:
    """Backbone from ff14SB serine; side chain AM1-BCC shifted evenly to a neutral residue."""
    side = [n for n in SIDE_CHAIN]
    shift = (-sum(backbone.values()) - sum(am1bcc[n] for n in side)) / len(side)
    return {**backbone, **{n: am1bcc[n] + shift for n in side}}, shift


def lj(structure: parmed.Structure) -> dict:
    """GAFF2 sigma (nm) / epsilon (kJ/mol) per side-chain OpenMM type."""
    table = {}
    for atom in structure.atoms:
        if atom.name in SIDE_CHAIN:
            sigma = 2 * atom.atom_type.rmin * ANGSTROM_TO_NM / 2 ** (1 / 6)
            add(table, (SIDE_CHAIN[atom.name],), (sigma, atom.atom_type.epsilon * KCAL_TO_KJ))
    return table


def xml(structure, found: dict, charges: dict) -> str:
    gaff = {a.name: a.type for a in structure.atoms}
    masses = {a.name: a.mass for a in structure.atoms}
    elements = {a.name: parmed.periodic_table.Element[a.atomic_number] for a in structure.atoms}
    root = ET.Element("ForceField")
    info = ET.SubElement(root, "Info")
    ET.SubElement(info, "Source").text = "tools/parameterize_ddz.py (see ddz_provenance.yaml)"
    types = ET.SubElement(root, "AtomTypes")
    for name in ("CB", "HB", "OG1", "HG1"):
        ET.SubElement(
            types,
            "Type",
            name=SIDE_CHAIN[name],
            **{"class": f"ddz-{gaff[name]}"},
            element=elements[name],
            mass=f"{masses[name]:.3f}",
        )
    residues = ET.SubElement(root, "Residues")
    residue = ET.SubElement(residues, "Residue", name="DDZ")
    for name in DDZ_ATOMS:
        ET.SubElement(
            residue, "Atom", name=name, type=openmm_type(name), charge=f"{charges[name]:.6f}"
        )
    for a, b in (
        ("N", "H"),
        ("N", "CA"),
        ("CA", "HA"),
        ("CA", "CB"),
        ("CA", "C"),
        ("C", "O"),
        ("CB", "HB"),
        ("CB", "OG1"),
        ("CB", "OG2"),
        ("OG1", "HG1"),
        ("OG2", "HG2"),
    ):
        ET.SubElement(residue, "Bond", atomName1=a, atomName2=b)
    ET.SubElement(residue, "ExternalBond", atomName="N")
    ET.SubElement(residue, "ExternalBond", atomName="C")
    bonds = ET.SubElement(root, "HarmonicBondForce")
    for (t1, t2), (length, k) in sorted(found["bonds"].items()):
        ET.SubElement(bonds, "Bond", type1=t1, type2=t2, length=str(length), k=str(k))
    angles = ET.SubElement(root, "HarmonicAngleForce")
    for (t1, t2, t3), (angle, k) in sorted(found["angles"].items()):
        ET.SubElement(angles, "Angle", type1=t1, type2=t2, type3=t3, angle=str(angle), k=str(k))
    torsions = ET.SubElement(root, "PeriodicTorsionForce")
    for key, parts in sorted(_torsions(found["torsions"]).items()):
        attributes = {f"type{i + 1}": t for i, t in enumerate(key)}
        for i, (per, phase, k) in enumerate(sorted(parts), start=1):
            attributes |= {f"periodicity{i}": str(per), f"phase{i}": str(phase), f"k{i}": str(k)}
        ET.SubElement(torsions, "Proper", **attributes)
    nonbonded = ET.SubElement(
        root, "NonbondedForce", coulomb14scale="0.8333333333333334", lj14scale="0.5"
    )
    ET.SubElement(nonbonded, "UseAttributeFromResidue", name="charge")
    for (type_,), (sigma, epsilon) in sorted(lj(structure).items()):
        ET.SubElement(nonbonded, "Atom", type=type_, sigma=str(sigma), epsilon=str(epsilon))
    ET.indent(root)
    return ET.tostring(root, encoding="unicode") + "\n"


def _torsions(raw: dict) -> dict:
    """Merge both directions of each torsion; the parameter sets must agree."""
    merged = {}
    for key, parts in raw.items():
        canonical = min(key, key[::-1])
        if canonical in merged and merged[canonical] != parts:
            raise RuntimeError(f"{canonical}: different torsion terms by direction")
        merged[canonical] = parts
    return merged


def versions() -> dict:
    import rdkit

    meta = sorted((Path(sys.prefix) / "conda-meta").glob("ambertools-*.json"))
    ambertools = meta[0].name.split("-")[1] if meta else "unknown"
    return {"ambertools": ambertools, "parmed": parmed.__version__, "rdkit": rdkit.__version__}


def main(openmm_data: Path) -> None:
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        structure, am1bcc = amber_parameters(fragment(work), work)
    backbone = serine(openmm_data)
    charges, shift = ddz_charges(am1bcc, backbone)
    found = terms(structure)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "ddz.xml").write_text(xml(structure, found, charges))
    provenance = {
        "generated_by": "tools/parameterize_ddz.py",
        "fragment_smiles": SMILES,
        "rdkit_seed": SEED,
        "tools": versions(),
        "am1bcc_charges": {n: round(am1bcc[n], 6) for n in DDZ_ATOMS},
        "side_chain_shift": round(shift, 6),
        "backbone_from": "OpenMM amber14/protein.ff14SB.xml, residue SER",
        "gaff2_types": {a.name: a.type for a in structure.atoms if a.name in DDZ_ATOMS},
        "terms": {k: len(v) for k, v in found.items()},
    }
    (OUT_DIR / "ddz_provenance.yaml").write_text(yaml.safe_dump(provenance, sort_keys=False))
    print(yaml.safe_dump(provenance, sort_keys=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
