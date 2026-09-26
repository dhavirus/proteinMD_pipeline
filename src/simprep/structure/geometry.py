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
