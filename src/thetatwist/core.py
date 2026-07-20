"""Exact calculation of the Section 12 theta twist."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence


Character = tuple[Fraction, Fraction, Fraction, Fraction, Fraction]
SCHEMA_VERSION = "thetatwist-certificate-1"


def _character(values: Iterable[int | Fraction]) -> Character:
    result = tuple(Fraction(value) for value in values)
    if len(result) != 5:
        raise ValueError("a theta character has five coefficients")
    return result  # type: ignore[return-value]


def _text(value: int | Fraction) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def character_text(character: Character) -> list[str]:
    return [_text(value) for value in character]


def product(left: Character, right: Character) -> Character:
    result = [Fraction(0)] * 5
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            if i + j <= 4:
                result[i + j] += x * y
    return _character(result)


def exponential(t: int | Fraction) -> Character:
    t = Fraction(t)
    return _character(t**degree / math.factorial(degree) for degree in range(5))


def twist(character: Character, t: int | Fraction) -> Character:
    return product(character, exponential(t))


def dual(character: Character) -> Character:
    return _character((-1) ** degree * value for degree, value in enumerate(character))


def normalized_character(character: Character) -> Character:
    if character[0] == 0:
        raise ValueError("normalization requires nonzero rank")
    return twist(character, -character[1] / character[0])


def self_euler(character: Character) -> int:
    value = 24 * product(dual(character), character)[4]
    if value.denominator != 1:
        raise ArithmeticError("self-Euler pairing is not integral")
    return value.numerator


def presentation_character(n: int, r: int, t: int) -> Character:
    """Character of the cyclic presentation with final twist ``O(t Theta)``."""

    if n < 6:
        raise ValueError("the displayed cyclic incidence formula assumes n >= 6")
    if r < 0:
        raise ValueError("r must be nonnegative")
    n_q = Fraction(n)
    t_q = Fraction(t)
    return _character(
        (
            1,
            t_q,
            t_q**2 / 2 - n_q,
            t_q**3 / 6 - n_q * t_q + 2 * n_q,
            t_q**4 / 24
            - n_q * t_q**2 / 2
            + 2 * n_q * t_q
            + n_q**2 / 2
            - 55 * n_q / 12
            - Fraction(r, 24),
        )
    )


def printed_parameters(d: int) -> tuple[int, int]:
    if d < 3 or d % 2 == 0:
        raise ValueError("Section 12 assumes odd d >= 3")
    return (d + 9) // 2, (d + 9) * (2 * d - 1)


def printed_character(d: int = 3) -> Character:
    n, r = printed_parameters(d)
    return presentation_character(n, r, 1)


def corrected_character(d: int = 3) -> Character:
    n, r = printed_parameters(d)
    return presentation_character(n, r, 3)


def secant_character(d: int, beta: int | Fraction) -> Character:
    if d <= 0:
        raise ValueError("d must be positive")
    beta = Fraction(beta)
    return _character((1, beta, -Fraction(d, 2), -Fraction(d, 6) * beta, Fraction(d * d, 24)))


def secant_membership(character: Character, d: int) -> dict[str, object]:
    alpha = character[0]
    beta = character[1]
    expected = _character(
        (
            alpha,
            beta,
            -Fraction(d, 2) * alpha,
            -Fraction(d, 6) * beta,
            Fraction(d * d, 24) * alpha,
        )
    )
    residual = _character(
        value - target for value, target in zip(character, expected, strict=True)
    )
    return {
        "member": all(value == 0 for value in residual),
        "alpha_coefficient": _text(alpha),
        "beta_coefficient": _text(beta),
        "residual": character_text(residual),
    }


def unique_repair(d: int = 3) -> dict[str, object]:
    n_printed, r_printed = printed_parameters(d)
    t = 3
    n_required = Fraction(d + t * t, 2)
    r_required = (
        t**4
        - 12 * n_required * t**2
        + 48 * n_required * t
        + 12 * n_required**2
        - 110 * n_required
        - d**2
    )
    corrected = corrected_character(d)
    expected = secant_character(d, 3)
    return {
        "degree2_constraint": "n=(d+t^2)/2",
        "degree3_factorization": "(t-3)(t^2+d)=0",
        "unique_positive_d_real_integral_twist": t,
        "printed_n": n_printed,
        "required_n": int(n_required),
        "printed_r": r_printed,
        "required_r": int(r_required),
        "counts_close": n_required == n_printed and r_required == r_printed,
        "proposed_repair": "O_X(Theta) -> O_X(3 Theta)",
        "same_nondual_operation": "tensor by O_X(2 Theta)",
        "corrected": character_text(corrected),
        "alpha_plus_3beta": character_text(expected),
        "coefficients_close": corrected == expected,
    }


def natural_repair_sweep(d: int = 3) -> dict[str, object]:
    """Calculate the results for the specified transformations."""

    printed = printed_character(d)
    candidates = (
        ("nondual_plus_2", "tensor by O_X(2 Theta)", twist(printed, 2)),
        ("nondual_minus_4", "tensor by O_X(-4 Theta)", twist(printed, -4)),
        ("dual_minus_2", "dualize, then tensor by O_X(-2 Theta)", twist(dual(printed), -2)),
        ("dual_plus_4", "dualize, then tensor by O_X(4 Theta)", twist(dual(printed), 4)),
    )
    return {
        "explicit_scope": (
            "theta twists, dualization, shifts/cone reversal, translation, Pic0 "
            "tensor, Todd/Mukai convention, multiplication pullback, parameter "
            "relabeling, component count, and point-normalization count"
        ),
        "integral_theta_twist_reduction": (
            "degree two gives k=2 and k=-4. Only k=2 agrees in degrees three and four"
        ),
        "candidates": [
            {
                "id": identifier,
                "operation": operation,
                "coefficients": character_text(character),
                "secant_membership": secant_membership(character, d),
                "self_euler": self_euler(character),
            }
            for identifier, operation, character in candidates
        ],
        "cohomological_no_ops": {
            "translation": "identity on rational cohomology",
            "Pic0_tensor": "Chern character one",
            "Todd_Mukai": "td(X)=1 for an abelian variety",
            "differential_sign": "a chain isomorphism gives the same K-class",
        },
        "cone_reversal": "negates rank and therefore cannot give the rank-one target",
        "point_normalization": "changes degree four only",
        "parameter_relabel_no_go": (
            "degrees two and three have no permitted common d-relabel solution"
        ),
        "multiplication_pullback_no_go": (
            "after c1 normalization, degree two requires "
            "(d+9)m^4=d+1. This equation has no integer solution when d>=3"
        ),
        "exact_alpha_plus_beta_autoequivalence_no_go": {
            "printed_euler": self_euler(printed),
            "target_euler": self_euler(secant_character(d, 1)),
            "reason": "derived autoequivalences keep the Euler pairing unchanged",
        },
    }


def _matrix_rank(matrix: Sequence[Sequence[Fraction]]) -> int:
    rows = [list(row) for row in matrix]
    pivot_row = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(pivot_row, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        divisor = rows[pivot_row][column]
        rows[pivot_row] = [value / divisor for value in rows[pivot_row]]
        for i in range(len(rows)):
            if i == pivot_row or rows[i][column] == 0:
                continue
            factor = rows[i][column]
            rows[i] = [
                value - factor * pivot_value
                for value, pivot_value in zip(rows[i], rows[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def _primitive_kernel_vector(matrix: tuple[tuple[Fraction, ...], ...]) -> tuple[int, int, int]:
    for first, second in ((0, 1), (0, 2), (1, 2)):
        a = matrix[first]
        b = matrix[second]
        vector = (
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        )
        if not any(vector):
            continue
        if not all(sum(row[i] * vector[i] for i in range(3)) == 0 for row in matrix):
            continue
        common = math.lcm(*(value.denominator for value in vector))
        integers = [value.numerator * (common // value.denominator) for value in vector]
        divisor = math.gcd(*[abs(value) for value in integers])
        integers = [value // divisor for value in integers]
        if next(value for value in reversed(integers) if value) < 0:
            integers = [-value for value in integers]
        return integers[0], integers[1], integers[2]
    raise ArithmeticError("expected a one-dimensional kernel")


def hodge_contraction(character: Character) -> dict[str, object]:
    c0, c1, c2, c3, c4 = character
    block = (
        (c0, 2 * c1, -2 * c2),
        (c1, 4 * c2, -6 * c3),
        (-2 * c2, -12 * c3, 24 * c4),
    )
    rank = _matrix_rank(block)
    kernel = _primitive_kernel_vector(block)
    return {
        "block": [[_text(value) for value in row] for row in block],
        "block_rank": rank,
        "full_rank": 6 * rank,
        "full_kernel_dimension": 28 - 6 * rank,
        "mixed_kernel_vector_A_Bantisym_C": list(kernel),
    }


def certificate(d: int = 3) -> dict[str, object]:
    printed = printed_character(d)
    corrected = corrected_character(d)
    older = secant_character(d, 1)
    return {
        "schema_version": SCHEMA_VERSION,
        "d": d,
        "source": "arXiv:2509.23403v2, Section 12, pages 21-22",
        "printed": {
            "coefficients": character_text(printed),
            "secant_membership": secant_membership(printed, d),
            "self_euler": self_euler(printed),
            "normalized": character_text(normalized_character(printed)),
            "hodge_contraction": hodge_contraction(printed),
        },
        "repair_certificate": unique_repair(d),
        "declared_natural_repair_sweep": natural_repair_sweep(d),
        "corrected": {
            "coefficients": character_text(corrected),
            "secant_membership": secant_membership(corrected, d),
            "self_euler": self_euler(corrected),
            "normalized": character_text(normalized_character(corrected)),
            "hodge_contraction": hodge_contraction(corrected),
        },
        "older_alpha_plus_beta_comparator": {
            "coefficients": character_text(older),
            "self_euler": self_euler(older),
            "source": "arXiv:2502.03415v2, Example 8.2.3",
        },
        "normalized_character_preserved": (
            normalized_character(printed) == normalized_character(corrected)
        ),
        "claim_boundary": (
            "This result applies only to the Chern character and Hodge contraction. "
            "It makes no claim about semiregularity, the Atiyah map, a group action, "
            "or the Hodge conjecture."
        ),
    }


def write_results(directory: str | Path, d: int = 3) -> dict[str, Path]:
    output = Path(directory)
    output.mkdir(parents=True, exist_ok=True)
    certificate_path = output / f"d{d}_certificate.json"
    certificate_path.write_text(
        json.dumps(certificate(d), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    manifest_path = output / "manifest.json"
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "artifacts": [
            {
                "path": certificate_path.name,
                "sha256": hashlib.sha256(certificate_path.read_bytes()).hexdigest(),
                "bytes": certificate_path.stat().st_size,
            }
        ],
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {"certificate": certificate_path, "manifest": manifest_path}
