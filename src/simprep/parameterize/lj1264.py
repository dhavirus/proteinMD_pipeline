"""The 12-6-4 ion model (TASK-010, ADR-0010).

C4 between the ion and an atom j = C4(ion, water O) / polarizability(water) x
polarizability(class of j) x tuning factor, as ParmEd's add12_6_4 computes it
(knowledge/forcefield/lj1264.yaml, from AmberTools). The ion's 12-6 terms are the 12-6-4
set's. ``c4_by_class`` is pure; ``apply`` edits an OpenMM System.
"""

from __future__ import annotations

KCAL_TO_KJ = 4.184
A4_TO_NM4 = 1e-4
ANGSTROM_TO_NM = 0.1
C4_ENERGY = "-(ion1*c42 + ion2*c41)/r^4"


class IonModelError(ValueError):
    """The 12-6-4 model cannot be applied as specified; the message says why."""


def c4_by_class(data: dict, ion: str) -> dict[str, float]:
    """C4 (kJ/mol nm^4) between ``ion`` and an atom of each OpenMM class."""
    scale = (
        data["ions"][ion]["c4_water_kcal_per_mol_a4"] / data["water_polarizability_cubic_angstrom"]
    )
    return {
        cls: scale * p["value_cubic_angstrom"] * data["tuning_factor"] * KCAL_TO_KJ * A4_TO_NM4
        for cls, p in data["polarizabilities"].items()
    }


def sigma_epsilon(data: dict, ion: str) -> tuple[float, float]:
    """12-6 terms of the 12-6-4 set as OpenMM sigma (nm) and epsilon (kJ/mol)."""
    terms = data["ions"][ion]
    sigma = 2 * terms["rmin_half_angstrom"] * ANGSTROM_TO_NM / 2 ** (1 / 6)
    return sigma, terms["epsilon_kcal_per_mol"] * KCAL_TO_KJ


def apply(system, particles: tuple, data: dict) -> int:
    """Give the ion particles (indices) the 12-6-4 model; ``particles`` is (ion residue name,
    ion particle indices, OpenMM class per particle). Returns the number of C4 pairs."""
    import openmm

    ion, indices, classes = particles
    c4 = c4_by_class(data, ion)
    missing = sorted({c for c in classes if c not in c4})
    if missing:
        raise IonModelError(f"no polarizability for classes {missing}; add them to lj1264.yaml")
    nonbonded = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    sigma, epsilon = sigma_epsilon(data, ion)
    for index in indices:
        charge, _, _ = nonbonded.getParticleParameters(index)
        nonbonded.setParticleParameters(index, charge, sigma, epsilon)
    force = openmm.CustomNonbondedForce(C4_ENERGY)
    force.addPerParticleParameter("ion")
    force.addPerParticleParameter("c4")
    ions = set(indices)
    for i, cls in enumerate(classes):
        force.addParticle([1.0 if i in ions else 0.0, c4[cls]])
    for exception in range(nonbonded.getNumExceptions()):
        first, second, *_ = nonbonded.getExceptionParameters(exception)
        force.addExclusion(first, second)
    force.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
    # ponytail: ion-ion C4 pairs are left out (ParmEd includes them); fine for one ion per
    # model (5FQL), add the ion-ion group when a system has several.
    others = [i for i in range(len(classes)) if i not in ions]
    force.addInteractionGroup(list(indices), others)
    system.addForce(force)
    return len(indices) * len(others)
