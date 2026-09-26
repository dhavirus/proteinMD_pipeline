"""Small, dependency-free geometry helpers (pure Python keeps detectors Pyodide-friendly)."""

from __future__ import annotations

import math
from collections.abc import Iterable

Point = tuple[float, float, float]


def distance(a: Point, b: Point) -> float:
    return math.dist(a, b)


def angle_degree(a: Point, vertex: Point, b: Point) -> float:
    """Angle a-vertex-b in degrees."""
    u = [a[i] - vertex[i] for i in range(3)]
    v = [b[i] - vertex[i] for i in range(3)]
    norm = math.hypot(*u) * math.hypot(*v)
    cosine = sum(u[i] * v[i] for i in range(3)) / norm
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def min_distance(points_a: Iterable[Point], points_b: Iterable[Point]) -> float | None:
    """Smallest distance between any point of ``points_a`` and any of ``points_b``."""
    list_b = list(points_b)
    distances = [distance(a, b) for a in points_a for b in list_b]
    return min(distances) if distances else None


def _sub(a: Point, b: Point) -> Point:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Point, b: Point) -> Point:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _unit(a: Point) -> Point:
    norm = math.hypot(*a)
    return (a[0] / norm, a[1] / norm, a[2] / norm)


def dihedral_degree(a: Point, b: Point, c: Point, d: Point) -> float:
    """Dihedral a-b-c-d in degrees, in (-180, 180]."""
    b1, b2, b3 = _sub(b, a), _sub(c, b), _sub(d, c)
    n1, n2 = _cross(b1, b2), _cross(b2, b3)
    x = sum(n1[i] * n2[i] for i in range(3))
    y = sum(c_i * u_i for c_i, u_i in zip(_cross(n1, n2), _unit(b2), strict=True))
    return math.degrees(math.atan2(y, x))


def place_atom(refs: tuple[Point, Point, Point], internal: tuple[float, float, float]) -> Point:
    """Position of d from a, b, c and (bond c-d, angle b-c-d, dihedral a-b-c-d) (NeRF).

    Angles in degrees. Natural extension reference frame: Parsons et al. 2005,
    J Comput Chem 26:1063-8, doi:10.1002/jcc.20237."""
    a, b, c = refs
    bond_angstrom, angle, torsion = internal
    theta, phi = math.radians(angle), math.radians(torsion)
    local = (
        -bond_angstrom * math.cos(theta),
        bond_angstrom * math.sin(theta) * math.cos(phi),
        bond_angstrom * math.sin(theta) * math.sin(phi),
    )
    bc = _unit(_sub(c, b))
    n = _unit(_cross(_sub(b, a), bc))
    m = _cross(n, bc)
    return tuple(c[i] + bc[i] * local[0] + m[i] * local[1] + n[i] * local[2] for i in range(3))
