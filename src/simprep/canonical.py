"""Canonical JSON (RFC 8785, JCS) so Python and the browser hash documents identically.

JCS serializes with sorted object keys, no whitespace, and numbers written the way
ECMAScript writes them (``1.0`` becomes ``1``, ``1e-7`` stays exponential). Python's
``json.dumps`` differs on numbers, so numbers are formatted here explicitly.
"""

from __future__ import annotations

import hashlib
import json
import math

# ECMAScript Number::toString switches to exponent notation outside this decimal range.
ES_MAX_FIXED_EXPONENT = 21
ES_MIN_FIXED_EXPONENT = -6


def canonical_json(value: object) -> str:
    """RFC 8785 serialization of a JSON value (dict, list, str, int, float, bool, None)."""
    if value is None or isinstance(value, bool):
        return json.dumps(value)
    if isinstance(value, int | float):
        return es_number(value)
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=False)
    if isinstance(value, list | tuple):
        return "[" + ",".join(canonical_json(item) for item in value) + "]"
    if isinstance(value, dict):
        items = sorted(value.items(), key=lambda item: _utf16_key(item[0]))
        return "{" + ",".join(f"{canonical_json(k)}:{canonical_json(v)}" for k, v in items) + "}"
    raise TypeError(f"not a JSON value: {type(value).__name__}")


def sha256_canonical(value: object) -> str:
    """SHA-256 (hex) of the UTF-8 canonical serialization."""
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


def _utf16_key(key: str) -> list[int]:
    """JCS orders keys by UTF-16 code units."""
    encoded = key.encode("utf-16-be")
    return [int.from_bytes(encoded[i : i + 2], "big") for i in range(0, len(encoded), 2)]


def es_number(number: int | float) -> str:
    """Format like ECMAScript ``Number.prototype.toString`` (shortest round-trip digits)."""
    if isinstance(number, int):
        return str(number)
    if not math.isfinite(number):
        raise ValueError(f"{number} is not representable in JSON")
    if number == 0:
        return "0"
    if number < 0:
        return "-" + es_number(-number)
    digits, point = _shortest_digits(number)
    return _es_layout(digits, point)


def _shortest_digits(number: float) -> tuple[str, int]:
    """Digits d1..dk and n such that number = 0.d1..dk x 10^n (repr is shortest round-trip)."""
    mantissa, _, exponent = repr(number).partition("e")
    whole, _, fraction = mantissa.partition(".")
    all_digits = (whole + fraction).lstrip("0")
    leading_zeros = len(whole + fraction) - len(all_digits)
    point = len(whole) - leading_zeros + (int(exponent) if exponent else 0)
    return all_digits.rstrip("0"), point


def _es_layout(digits: str, point: int) -> str:
    """ECMAScript Number::toString layout for positive numbers (ECMA-262, 6.1.6.1.20)."""
    k = len(digits)
    if k <= point <= ES_MAX_FIXED_EXPONENT:
        return digits + "0" * (point - k)
    if 0 < point <= ES_MAX_FIXED_EXPONENT:
        return digits[:point] + "." + digits[point:]
    if ES_MIN_FIXED_EXPONENT < point <= 0:
        return "0." + "0" * (-point) + digits
    exponent = point - 1
    sign = "+" if exponent >= 0 else "-"
    mantissa = digits if k == 1 else digits[0] + "." + digits[1:]
    return f"{mantissa}e{sign}{abs(exponent)}"
