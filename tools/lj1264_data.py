"""Extract the 12-6-4 Ca2+ model data from AmberTools, and ParmEd reference C4 values
(TASK-010, ADR-0010).

Run once with AmberTools on PATH (AMBERHOME set): ``ENV/bin/python tools/lj1264_data.py
OPENMM_DATA_DIR``. Writes ``knowledge/forcefield/lj1264.yaml`` (Ca2+ 12-6-4 TIP3P terms,
C4 with water, water polarizability, polarizability per OpenMM atom class, with sources)
and ``tests/panel/parameter_fixtures/lj1264_reference.yaml`` (C4 between Ca2+ and each
atom type of a small tleap system after ParmEd's add12_6_4: the independent reference
simprep's C4 term is tested against). Not run in CI; the outputs are committed data.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import parmed
import yaml
from parmed.tools import add12_6_4
from parmed.tools.add1264 import DEFAULT_C4_PARAMS, WATER_POL

ROOT = Path(__file__).resolve().parents[1]
ION_FILE = "frcmod.ions234lm_1264_tip3p"
POL_FILE = "lj_1264_pol.dat"
ION = "Ca2+"
C4_KEY = "Ca2"
WATER_MODEL = "TIP3P"
TUNING_FACTOR = 1.0
# OpenMM classes whose Amber type has another name; every other class is its Amber type.
CLASS_TO_AMBER = {
    "tip3p-O": "OW",
    "tip3p-H": "HW",
    "tip3p_standard-Cl-": "Cl-",
    "tip3p_standard-Ca2+": "Ca2+",
    "ddz-c3": "c3",
    "ddz-h2": "h2",
    "ddz-oh": "oh",
    "ddz-ho": "ho",
}
LEAP = """source leaprc.protein.ff14SB
source leaprc.water.tip3p
loadamberparams {ion_file}
p = sequence {{ ACE ASP HIE LYS SER CYS NME }}
m = combine {{ p CA CL WAT }}
saveamberparm m small.prmtop small.inpcrd
quit
"""


def parm_dir() -> Path:
    return Path(os.environ["AMBERHOME"]) / "dat" / "leap" / "parm"


def ion_terms() -> dict:
    """Rmin/2 and epsilon of Ca2+ from the 12-6-4 TIP3P frcmod, with its comment."""
    for line in (parm_dir() / ION_FILE).read_text().splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0] == ION and "." in fields[1]:
            return {
                "rmin_half_angstrom": float(fields[1]),
                "epsilon_kcal_per_mol": float(fields[2]),
                "source": f"AmberTools {ION_FILE}: {' '.join(fields[3:])}",
            }
    raise RuntimeError(f"{ION} not found in {ION_FILE}")


def polarizabilities(classes: list[str]) -> dict:
    table = {}
    for line in (parm_dir() / POL_FILE).read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2:
            table.setdefault(fields[0], (float(fields[1]), " ".join(fields[2:])))
    missing = [c for c in classes if CLASS_TO_AMBER.get(c, c) not in table]
    if missing:
        raise RuntimeError(f"no polarizability for {missing}")
    return {
        c: {
            "amber_type": CLASS_TO_AMBER.get(c, c),
            "value_cubic_angstrom": table[CLASS_TO_AMBER.get(c, c)][0],
            "source": table[CLASS_TO_AMBER.get(c, c)][1],
        }
        for c in classes
    }


def openmm_classes(openmm_data: Path) -> list[str]:
    classes = set()
    for path in (
        openmm_data / "amber14" / "protein.ff14SB.xml",
        ROOT / "knowledge" / "forcefield" / "ddz.xml",
    ):
        classes |= {t.get("class") for t in ET.parse(path).iter("Type")}
    return sorted(classes | {"tip3p-O", "tip3p-H", "tip3p_standard-Cl-", "tip3p_standard-Ca2+"})


def reference() -> dict:
    """C4 (kcal/mol/A^4) between Ca2+ and each atom type of a tleap system, by ParmEd."""
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        (work / "leap.in").write_text(LEAP.format(ion_file=ION_FILE))
        subprocess.run(["tleap", "-f", "leap.in"], cwd=work, check=True, capture_output=True)
        parm = parmed.amber.AmberParm(str(work / "small.prmtop"), str(work / "small.inpcrd"))
    add12_6_4(parm, "@%Ca2+", "watermodel", WATER_MODEL, "tunfactor", str(TUNING_FACTOR)).execute()
    ntypes = parm.ptr("ntypes")
    ion_index = next(a.nb_idx for a in parm.atoms if a.type == ION)
    c4 = {}
    for atom in parm.atoms:
        i, j = sorted((ion_index, atom.nb_idx))
        index = parm.parm_data["NONBONDED_PARM_INDEX"][ntypes * (i - 1) + j - 1] - 1
        value = parm.parm_data["LENNARD_JONES_CCOEF"][index]
        if atom.type in c4 and abs(c4[atom.type] - value) > 1e-9:
            raise RuntimeError(f"type {atom.type}: two C4 values")
        c4[atom.type] = value
    return {t: round(v, 6) for t, v in sorted(c4.items())}


def main(openmm_data: Path) -> None:
    classes = openmm_classes(openmm_data)
    meta = sorted((Path(sys.prefix) / "conda-meta").glob("ambertools-*.json"))
    data = {
        "schema_version": "0.1.0",
        "generated_by": "tools/lj1264_data.py",
        "ambertools": meta[0].name.split("-")[1] if meta else "unknown",
        "parmed": parmed.__version__,
        "ions": {
            "CA": {
                "openmm_class": "tip3p_standard-Ca2+",
                **ion_terms(),
                "c4_water_kcal_per_mol_a4": DEFAULT_C4_PARAMS[WATER_MODEL][C4_KEY],
                "c4_source": (
                    f"ParmEd {parmed.__version__} parmed/tools/add1264.py, {WATER_MODEL} table"
                ),
            }
        },
        "water_polarizability_cubic_angstrom": WATER_POL,
        "tuning_factor": TUNING_FACTOR,
        "polarizabilities": polarizabilities(classes),
    }
    out = ROOT / "knowledge" / "forcefield" / "lj1264.yaml"
    out.write_text(yaml.safe_dump(data, sort_keys=False))
    fixture = ROOT / "tests" / "panel" / "parameter_fixtures"
    fixture.mkdir(parents=True, exist_ok=True)
    ref = {
        "note": "C4 between Ca2+ and each Amber atom type of a tleap system "
        "(ACE-ASP-HIE-LYS-SER-CYS-NME, Ca2+, Cl-, one TIP3P water) after ParmEd add12_6_4 "
        f"(watermodel {WATER_MODEL}, tunfactor {TUNING_FACTOR}); written by tools/lj1264_data.py, "
        "independent of simprep's implementation.",
        "c4_kcal_per_mol_a4": reference(),
    }
    (fixture / "lj1264_reference.yaml").write_text(yaml.safe_dump(ref, sort_keys=False))
    print(yaml.safe_dump(ref, sort_keys=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
