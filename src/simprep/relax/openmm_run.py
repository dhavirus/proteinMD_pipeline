"""Restrained local minimization with OpenMM (I/O edge; TASK-006, ADR-0006).

OpenMM is imported here only, lazily, so simprep runs without it until relaxation is
used (``pip install simprep[relax]``). Hydrogens are added for the calculation only
(seeded, on the Reference platform, so runs are byte-identical) and discarded: residues
whose heavy atoms moved come back without hydrogens.

With torsion restraints (loop modelling, TASK-007), a bonded-only pre-stage runs first:
the same system without non-bonded terms, so a placement with inverted chirality or cis
peptides can reach the intended basin without clashes in the way. The torsion restraints
stay on in the final minimization, or clash forces flip the centres back. With
``nonbonded: sterics_only`` the final minimization has Lennard-Jones between mobile atoms
and everything else instead of the full non-bonded terms (vacuum electrostatics stretch
a charged loop; ADR-0007).
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, replace

from simprep.relax.shell import (
    ANGSTROM_PER_NM,
    AtomKey,
    RelaxError,
    Shell,
    heavy,
    polymer_and_water,
    segments,
    untemplated,
)
from simprep.structure.model import Residue, ResidueId, Structure

COORDINATE_DECIMALS = 3  # the precision of the written files
MOVED_ANGSTROM = 0.001
RESTRAINT = "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
STERICS_ONLY = "sterics_only"
STERICS = (
    "4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); "
    "epsilon=sqrt(epsilon1*epsilon2)"
)
TORSION_RESTRAINT = "0.5*k_torsion*min(dt, 2*pi-dt)^2; dt=abs(theta-theta0); pi=3.141592653589793"


@dataclass(frozen=True)
class TorsionRestraint:
    """Harmonic restraint of the dihedral over four heavy atoms to ``target_degree``
    (pre-stage only)."""

    atoms: tuple[AtomKey, AtomKey, AtomKey, AtomKey]
    target_degree: float


@dataclass(frozen=True)
class RelaxOutcome:
    structure: Structure
    moved_residues: tuple[ResidueId, ...]
    moved_atoms: int
    max_displacement_angstrom: float
    rms_displacement_angstrom: float
    energy_kj_per_mol: tuple[float, float]
    max_mobile_force_kj_per_mol_nm: float
    left_out: tuple[dict, ...]
    openmm_version: str


def _openmm():
    try:
        import openmm
        from openmm import app, unit
    except ImportError as error:
        raise RelaxError(
            "relaxation needs OpenMM: pip install 'simprep[relax]' (openmm==8.6.1), or set "
            "the manifest's relaxation.enabled to false"
        ) from error
    return openmm, app, unit


def relax(
    structure: Structure, shell: Shell, protocol: dict, torsions: tuple[TorsionRestraint, ...] = ()
) -> RelaxOutcome:
    """Minimize ``shell`` in ``structure`` under ``protocol`` (manifest ``relaxation``, or
    the modelling protocol's ``minimization``, which adds
    ``torsion_restraint_kj_per_mol_rad2`` for ``torsions``)."""
    with_altlocs = [r.id.label() for r in structure.residues if r.altlocs]
    if with_altlocs:
        raise RelaxError(
            f"relaxation needs one conformer per residue; altlocs remain at "
            f"{', '.join(with_altlocs[:5])}{' ...' if len(with_altlocs) > 5 else ''} "
            "(decide them in the manifest; prep collapses them)"
        )
    left_out, blocking = untemplated(structure, shell, protocol["nonbonded_cutoff_nm"])
    if blocking:
        names = ", ".join(f"{r.name} {r.id.label()} ({d} A)" for r, d in blocking)
        raise RelaxError(
            f"residues without a force-field template lie within the non-bonded cutoff of "
            f"the mobile shell: {names}. Parameterize them first (TASK-006 decision 4)."
        )
    openmm, app, unit = _openmm()
    model = _model(structure, (openmm, app, unit), protocol)
    minimized = _minimize(model, shell, (protocol, torsions))
    return _outcome(structure, model, minimized, (left_out, openmm.__version__))


@dataclass(frozen=True)
class Model:
    """The OpenMM side: topology with hydrogens, positions, and our residue per residue."""

    topology: object
    positions: object
    residue_ids: list[ResidueId]
    openmm: tuple


@dataclass(frozen=True)
class Minimized:
    positions: object
    energies: tuple
    forces: object
    mobile: list[int]


def _model(structure: Structure, modules: tuple, protocol: dict) -> Model:
    openmm, app, unit = modules
    topology, positions, residue_ids = app.Topology(), [], []
    polymer, waters = polymer_and_water(structure)
    for segment in [*segments(polymer), waters]:
        chain = topology.addChain()
        for residue in segment:
            _add_residue(topology, chain, residue, (app, positions))
            residue_ids.append(residue.id)
    topology.createStandardBonds()
    quantity = unit.Quantity([openmm.Vec3(*p) for p in positions], unit.angstrom)
    topology.createDisulfideBonds(quantity)
    modeller = app.Modeller(topology, quantity)
    random.seed(protocol["seed"])  # OpenMM places new hydrogens at random positions
    reference = openmm.Platform.getPlatformByName(protocol["platform"])
    modeller.addHydrogens(pH=protocol["hydrogen_ph"], platform=reference)
    return Model(modeller.topology, modeller.positions, residue_ids, modules)


def _add_residue(topology, chain, residue: Residue, sink: tuple) -> None:
    app, positions = sink
    new = topology.addResidue(residue.name, chain, id=f"{residue.id.seq_num}{residue.id.ins_code}")
    for atom in heavy(residue):
        topology.addAtom(atom.name, app.Element.getBySymbol(atom.element.title()), new)
        positions.append(atom.position)


def _minimize(model: Model, shell: Shell, settings: tuple):
    protocol, torsions = settings
    openmm, app, unit = model.openmm
    forcefield = app.ForceField(*protocol["force_field_files"])
    system = forcefield.createSystem(
        model.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=protocol["nonbonded_cutoff_nm"] * unit.nanometer,
        constraints=None,
        rigidWater=False,  # rigid waters among fixed (massless) atoms stop the minimizer
        ignoreExternalBonds=True,  # chains are split at gaps
        residueTemplates=_disulfide_templates(model.topology),
    )
    keys = _atom_keys(model)
    mobile = [a.index for a in model.topology.atoms() if _is_mobile(a, keys, shell)]
    if protocol.get("nonbonded") == STERICS_ONLY:
        _sterics_only(system, mobile, (openmm, unit, model.positions, protocol))
    _restrain(system, keys, shell, (openmm, model, protocol))
    if torsions:
        system.addForce(_torsion_force(openmm, keys, torsions, protocol))
    for atom in model.topology.atoms():
        if atom.index not in mobile:
            system.setParticleMass(atom.index, 0.0)
    platform = openmm.Platform.getPlatformByName(protocol["platform"])
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
    context.setPositions(model.positions)
    initial = context.getState(getEnergy=True).getPotentialEnergy()
    if torsions:
        context.setPositions(_bonded_prestage(system, model, protocol))
    openmm.LocalEnergyMinimizer.minimize(
        context, protocol["tolerance_kj_per_mol_nm"], protocol["max_iterations"]
    )
    state = context.getState(getEnergy=True, getPositions=True, getForces=True)
    _check_margin(model.positions, state.getPositions(), (mobile, protocol, unit))
    return Minimized(
        state.getPositions(), (initial, state.getPotentialEnergy()), state.getForces(), mobile
    )


def _sterics_only(system, mobile: list[int], context: tuple) -> None:
    """Replace the NonbondedForce by Lennard-Jones between mobile atoms and every atom
    within the cutoff plus ``sterics_margin_nm`` of their start: no charges, no 1-2/1-3/1-4
    pairs, no fixed-fixed pairs (loop modelling; exact while no atom moves farther than
    the margin, which _check_margin enforces)."""
    openmm, unit, positions, protocol = context
    index = next(
        i for i, f in enumerate(system.getForces()) if isinstance(f, openmm.NonbondedForce)
    )
    nonbonded = system.getForce(index)
    sterics = openmm.CustomNonbondedForce(STERICS)
    for parameter in ("sigma", "epsilon"):
        sterics.addPerParticleParameter(parameter)
    for particle in range(nonbonded.getNumParticles()):
        _, sigma, epsilon = nonbonded.getParticleParameters(particle)
        sterics.addParticle(
            [sigma.value_in_unit(unit.nanometer), epsilon.value_in_unit(unit.kilojoule_per_mole)]
        )
    for exception in range(nonbonded.getNumExceptions()):
        first, second, *_ = nonbonded.getExceptionParameters(exception)
        sterics.addExclusion(first, second)
    sterics.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
    sterics.setCutoffDistance(protocol["sterics_cutoff_nm"])
    sterics.setUseSwitchingFunction(True)
    sterics.setSwitchingDistance(protocol["sterics_switch_nm"])
    reach_nm = protocol["sterics_cutoff_nm"] + protocol["sterics_margin_nm"]
    sterics.addInteractionGroup(mobile, _within(positions, mobile, (reach_nm, unit)))
    system.removeForce(index)
    system.addForce(sterics)


def _within(positions, mobile: list[int], reach: tuple) -> list[int]:
    """Indices of atoms within ``reach`` (nm, unit module) of any mobile atom."""
    import numpy

    reach_nm, unit = reach
    xyz = numpy.array(positions.value_in_unit(unit.nanometer))
    nearest = numpy.full(len(xyz), numpy.inf)
    for index in mobile:
        nearest = numpy.minimum(nearest, numpy.linalg.norm(xyz - xyz[index], axis=1))
    return [int(i) for i in numpy.flatnonzero(nearest <= reach_nm)]


def _check_margin(start, final, context: tuple) -> None:
    """Raise if a mobile atom moved farther than the sterics margin (the partner list
    would then miss pairs within the cutoff)."""
    mobile, protocol, unit = context
    if protocol.get("nonbonded") != STERICS_ONLY:
        return
    before, after = (p.value_in_unit(unit.nanometer) for p in (start, final))
    moved_nm = max(math.dist(before[i], after[i]) for i in mobile)
    if moved_nm > protocol["sterics_margin_nm"]:
        raise RelaxError(
            f"a mobile atom moved {moved_nm:.2f} nm, more than sterics_margin_nm "
            f"({protocol['sterics_margin_nm']}); raise it in the protocol and run again"
        )


def _bonded_prestage(system, model: Model, protocol: dict):
    """Minimize a copy of ``system`` without its non-bonded terms; return the positions."""
    openmm = model.openmm[0]
    staged = openmm.XmlSerializer.clone(system)
    for index in reversed(range(staged.getNumForces())):
        if isinstance(staged.getForce(index), openmm.NonbondedForce | openmm.CustomNonbondedForce):
            staged.removeForce(index)
    platform = openmm.Platform.getPlatformByName(protocol["platform"])
    context = openmm.Context(staged, openmm.VerletIntegrator(0.001), platform)
    context.setPositions(model.positions)
    openmm.LocalEnergyMinimizer.minimize(
        context, protocol["tolerance_kj_per_mol_nm"], protocol["max_iterations"]
    )
    return context.getState(getPositions=True).getPositions()


def _torsion_force(openmm, keys: dict, torsions: tuple, protocol: dict):
    index = {(rid, name): i for i, (rid, name, is_heavy) in keys.items() if is_heavy}
    force = openmm.CustomTorsionForce(TORSION_RESTRAINT)
    force.addGlobalParameter("k_torsion", protocol["torsion_restraint_kj_per_mol_rad2"])
    force.addPerTorsionParameter("theta0")
    for torsion in torsions:
        force.addTorsion(
            *(index[key] for key in torsion.atoms), [math.radians(torsion.target_degree)]
        )
    return force


def _disulfide_templates(topology) -> dict:
    """CYX for disulfide cysteines: with external bonds ignored, CYM would match too."""
    return {
        atom.residue: "CYX"
        for bond in topology.bonds()
        if bond[0].name == "SG" and bond[1].name == "SG"
        for atom in bond
    }


def _atom_keys(model: Model) -> dict:
    """OpenMM atom index -> (ResidueId, heavy-atom name); a hydrogen maps to its heavy atom."""
    parents = {}
    for first, second in model.topology.bonds():
        for hydrogen, other in ((first, second), (second, first)):
            if hydrogen.element.symbol == "H":
                parents[hydrogen.index] = other
    keys = {}
    for atom in model.topology.atoms():
        owner = parents.get(atom.index, atom)
        keys[atom.index] = (model.residue_ids[atom.residue.index], owner.name, atom is owner)
    return keys


def _is_mobile(atom, keys: dict, shell: Shell) -> bool:
    rid, name, _ = keys[atom.index]
    return (rid, name) in shell.mobile


def _restrain(system, keys: dict, shell: Shell, context: tuple) -> None:
    """Harmonic restraint to the starting position on mobile heavy atoms, except the
    mutated residues' side chains."""
    openmm, model, protocol = context
    _, _, unit = model.openmm
    start = model.positions.value_in_unit(unit.nanometer)
    force = openmm.CustomExternalForce(RESTRAINT)
    force.addGlobalParameter("k", protocol["restraint_kj_per_mol_nm2"])
    for parameter in ("x0", "y0", "z0"):
        force.addPerParticleParameter(parameter)
    for index, (rid, name, is_heavy) in sorted(keys.items()):
        if is_heavy and (rid, name) in shell.mobile and (rid, name) not in shell.unrestrained:
            force.addParticle(index, list(start[index]))
    system.addForce(force)


def _outcome(
    structure: Structure, model: Model, minimized: Minimized, extra: tuple
) -> RelaxOutcome:
    left_out, version = extra
    _, _, unit = model.openmm
    new = _heavy_positions(model, minimized.positions, unit)
    if not all(math.isfinite(c) for position in new.values() for c in position):
        raise RelaxError("minimization produced non-finite coordinates; nothing was changed")
    residues, moved, displacements = [], [], []
    for residue in structure.residues:
        updated, shifts = _updated(residue, new)
        residues.append(updated)
        displacements += shifts
        if any(shift >= MOVED_ANGSTROM for shift in shifts):
            moved.append(residue.id)
    moving = [d for d in displacements if d >= MOVED_ANGSTROM]
    initial, final = (e.value_in_unit(unit.kilojoule_per_mole) for e in minimized.energies)
    forces = minimized.forces.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
    return RelaxOutcome(
        structure=replace(structure, residues=tuple(residues)),
        moved_residues=tuple(moved),
        moved_atoms=len(moving),
        max_displacement_angstrom=max(displacements, default=0.0),
        rms_displacement_angstrom=math.sqrt(sum(d * d for d in moving) / max(1, len(moving))),
        energy_kj_per_mol=(initial, final),
        max_mobile_force_kj_per_mol_nm=max(
            (math.hypot(*forces[i]) for i in minimized.mobile), default=0.0
        ),
        left_out=tuple(
            {**r.id.to_dict(), "res_name": r.name, "distance_angstrom": d} for r, d in left_out
        ),
        openmm_version=version,
    )


def _heavy_positions(model: Model, positions, unit) -> dict:
    """(ResidueId, atom name) -> relaxed position in A, rounded like the written files."""
    keys = _atom_keys(model)
    values = positions.value_in_unit(unit.nanometer)
    return {
        (rid, name): tuple(round(c * ANGSTROM_PER_NM, COORDINATE_DECIMALS) for c in values[i])
        for i, (rid, name, is_heavy) in keys.items()
        if is_heavy
    }


def _updated(residue: Residue, new: dict) -> tuple[Residue, list[float]]:
    """The residue with relaxed heavy atoms; without hydrogens if any heavy atom moved."""
    shifts, atoms = [], []
    for atom in heavy(residue):
        position = new.get((residue.id, atom.name), atom.position)
        shifts.append(math.dist(position, atom.position))
        atoms.append(replace(atom, position=position))
    if not any(shift >= MOVED_ANGSTROM for shift in shifts):
        return residue, shifts
    return replace(residue, atoms=tuple(atoms)), shifts
