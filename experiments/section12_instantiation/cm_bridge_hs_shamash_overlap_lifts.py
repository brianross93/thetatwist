"""Construct the exact low-degree Shamash overlap comparison.

This producer works in the exact Hermite decomposition of every overlap
algebra into dual-number factors.  Polynomial/Laurent exponents remain in
the source affine coordinates; only the finite coefficient algebra is put in
its two-by-two local basis.  The standard-monomial maps are recovered by the
stored chart and overlap Hermite matrices.

For the target endpoint it constructs the graph compounds ``J1_K``,
``J2_K``, and ``J3_K`` from the rational divided-difference block.  It then
uses the diagonal six-equation change matrix ``B`` and the ordered Laurent
Koszul contraction to construct

    Delta_j = J1_K H'_{j,0} - H^e_{j,0} B_j J0,
    Gamma_j = -h_1(Delta_j),

    Xi_j = J2_K H'_{j,1} + Gamma_j d'_1
           - H^e_{j,1} B_j J1_K,
    Omega_j = h_2(Xi_j).

The exact checks are ``d2 Gamma_j=Delta_j`` and
``d3 Omega_j=Xi_j``.  These give the full raw maps

    J2 = [[J2_K, Gamma], [0, B tensor J0]],
    J3 = [[J3_K, Omega], [0, B tensor J1_K]]

through the Ext^2 window.  No transported Hartshorne--Serre row, bounded
projective replacement, endpoint map, or Yoneda value is inferred.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import multiprocessing
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_graph_koszul_j1_overlap_lifts as j1_lifts,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_iw_presentation_solver as iw_solver,
)


P = 31
VARIABLE_COUNT = 8
FACTOR_COUNT = 4
SURFACE_COUNT = 6
EDGE_COUNT = 32
ZERO_EXPONENT = (0,) * VARIABLE_COUNT

INPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_iw_hermite_border_basis"
)
OVERLAP_DIRECTORY = INPUT_DIRECTORY / "overlaps"
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_shamash_overlap_lifts"
)

WEDGES = {
    degree: tuple(itertools.combinations(range(VARIABLE_COUNT), degree))
    for degree in range(4)
}
WEDGE_INDEX = {
    degree: {wedge: index for index, wedge in enumerate(WEDGES[degree])}
    for degree in range(4)
}


class HSShamashOverlapError(RuntimeError):
    """A Laurent identity, local homotopy, or bound input changed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def array_sha256(value: np.ndarray) -> str:
    return sha256_bytes(
        np.ascontiguousarray(np.asarray(value)).tobytes()
    )


def canonical_sha256(value) -> str:
    return sha256_bytes(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode(
            "utf-8"
        )
    )


def write_csv(path: Path, rows, fieldnames=None) -> None:
    rows = list(rows)
    if not rows and fieldnames is None:
        raise ValueError("fieldnames are required for an empty CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames or list(rows[0].keys())
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json_atomic(path: Path, value) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _shift(exponent, variable, amount):
    result = list(map(int, exponent))
    result[int(variable)] += int(amount)
    return tuple(result)


def _sum_exponents(left, right):
    return tuple(int(a) + int(b) for a, b in zip(left, right))


def _add_poly_term(polynomial, exponent, matrix, scale=1) -> None:
    exponent = tuple(map(int, exponent))
    value = int(scale) * np.asarray(matrix, dtype=np.int64) % P
    if exponent in polynomial:
        value = (polynomial[exponent] + value) % P
    if np.any(value):
        polynomial[exponent] = value
    else:
        polynomial.pop(exponent, None)


def _add_poly(target, source, scale=1, shift=ZERO_EXPONENT) -> None:
    for exponent, matrix in source.items():
        _add_poly_term(
            target, _sum_exponents(exponent, shift), matrix, scale
        )


def _compose_poly(left, right):
    result = {}
    for left_exponent, left_matrix in left.items():
        for right_exponent, right_matrix in right.items():
            _add_poly_term(
                result,
                _sum_exponents(left_exponent, right_exponent),
                left_matrix @ right_matrix,
            )
    return result


def _add_block(target, output_wedge, input_wedge, polynomial, scale=1) -> None:
    if not polynomial:
        return
    key = (tuple(output_wedge), tuple(input_wedge))
    value = target.setdefault(key, {})
    _add_poly(value, polynomial, scale)
    if not value:
        target.pop(key, None)


def add_block_maps(left, right, right_scale=1):
    result = {}
    for key, polynomial in left.items():
        _add_block(result, key[0], key[1], polynomial)
    for key, polynomial in right.items():
        _add_block(
            result, key[0], key[1], polynomial, scale=right_scale
        )
    return result


def shift_block_map(blocks, exponent, scale=1):
    result = {}
    for key, polynomial in blocks.items():
        shifted = {}
        _add_poly(shifted, polynomial, scale=scale, shift=exponent)
        _add_block(result, key[0], key[1], shifted)
    return result


def compose_block_maps(left, right):
    """Compose block maps ``left after right``."""

    result = {}
    right_by_output = {}
    for (middle, input_wedge), polynomial in right.items():
        right_by_output.setdefault(middle, []).append(
            (input_wedge, polynomial)
        )
    for (output_wedge, middle), left_polynomial in left.items():
        for input_wedge, right_polynomial in right_by_output.get(middle, []):
            _add_block(
                result,
                output_wedge,
                input_wedge,
                _compose_poly(left_polynomial, right_polynomial),
            )
    return result


def block_map_sha256(blocks) -> str:
    digest = hashlib.sha256()
    for key in sorted(blocks):
        digest.update(json.dumps(key, separators=(",", ":")).encode("ascii"))
        for exponent in sorted(blocks[key]):
            digest.update(
                json.dumps(exponent, separators=(",", ":")).encode("ascii")
            )
            digest.update(
                np.ascontiguousarray(
                    np.asarray(blocks[key][exponent], dtype=np.uint8)
                ).tobytes()
            )
    return digest.hexdigest()


def block_map_term_count(blocks) -> int:
    return sum(len(polynomial) for polynomial in blocks.values())


def block_map_scalar_nonzeros(blocks) -> int:
    return sum(
        int(np.count_nonzero(matrix % P))
        for polynomial in blocks.values()
        for matrix in polynomial.values()
    )


def assert_block_maps_equal(left, right, label) -> None:
    keys = set(left) | set(right)
    for key in keys:
        left_poly = left.get(key, {})
        right_poly = right.get(key, {})
        exponents = set(left_poly) | set(right_poly)
        for exponent in exponents:
            left_matrix = left_poly.get(exponent)
            right_matrix = right_poly.get(exponent)
            if left_matrix is None:
                left_matrix = np.zeros_like(right_matrix)
            if right_matrix is None:
                right_matrix = np.zeros_like(left_matrix)
            if not np.array_equal(left_matrix % P, right_matrix % P):
                raise HSShamashOverlapError(
                    f"{label} failed at block {key}, exponent {exponent}"
                )


def matrix_power(matrix, exponent, inverse=None):
    exponent = int(exponent)
    if exponent < 0:
        if inverse is None:
            raise HSShamashOverlapError(
                "a negative coefficient-algebra power lacks an inverse"
            )
        return matrix_power(inverse, -exponent)
    result = np.eye(matrix.shape[0], dtype=np.int64)
    factor = np.asarray(matrix, dtype=np.int64) % P
    while exponent:
        if exponent & 1:
            result = result @ factor % P
        factor = factor @ factor % P
        exponent //= 2
    return result


def inverse_2x2(matrix):
    matrix = np.asarray(matrix, dtype=np.int64) % P
    determinant = (
        int(matrix[0, 0]) * int(matrix[1, 1])
        - int(matrix[0, 1]) * int(matrix[1, 0])
    ) % P
    if not determinant:
        raise HSShamashOverlapError("the localized coordinate is not invertible")
    inverse_determinant = pow(determinant, -1, P)
    inverse = inverse_determinant * np.asarray(
        [[matrix[1, 1], -matrix[0, 1]], [-matrix[1, 0], matrix[0, 0]]],
        dtype=np.int64,
    ) % P
    if not np.array_equal(matrix @ inverse % P, np.eye(2, dtype=np.int64)):
        raise HSShamashOverlapError("the local two-by-two inverse failed")
    return inverse


def local_multiplication(values, tangents):
    result = np.zeros((VARIABLE_COUNT, 2, 2), dtype=np.int64)
    for variable, (value, tangent) in enumerate(zip(values, tangents)):
        result[variable] = np.asarray(
            [[int(value), 0], [int(tangent), int(value)]], dtype=np.int64
        ) % P
    return result


def evaluate_matrix_polynomial(polynomial, multiplication, inverse_variable=None):
    result = np.zeros((2, 2), dtype=np.int64)
    inverse = None
    if inverse_variable is not None:
        inverse = inverse_2x2(multiplication[inverse_variable])
    for exponent, coefficient in polynomial.items():
        operator = np.eye(2, dtype=np.int64)
        for variable, power in enumerate(exponent):
            operator = (
                operator
                @ matrix_power(
                    multiplication[variable],
                    power,
                    inverse if variable == inverse_variable else None,
                )
            ) % P
        result = (result + int(coefficient) * operator) % P
    return result


def divided_difference_matrix(polynomial, multiplication):
    """Return the ordered H0 map ``K0 -> K1`` as polynomial blocks."""

    maximum = [
        max([int(exponent[variable]) for exponent in polynomial] or [0])
        for variable in range(VARIABLE_COUNT)
    ]
    powers = {
        (variable, power): matrix_power(
            multiplication[variable], power
        )
        for variable in range(VARIABLE_COUNT)
        for power in range(maximum[variable] + 1)
    }
    identity = np.eye(2, dtype=np.int64)
    result = {}
    for exponent, coefficient in polynomial.items():
        exponent = tuple(map(int, exponent))
        for variable, power in enumerate(exponent):
            if not power:
                continue
            prefix = identity
            for earlier in range(variable):
                prefix = (
                    prefix @ powers[(earlier, exponent[earlier])]
                ) % P
            for quotient_power in range(power):
                matrix = (
                    int(coefficient)
                    * prefix
                    @ powers[(variable, quotient_power)]
                ) % P
                output_exponent = [0] * VARIABLE_COUNT
                output_exponent[variable] = power - 1 - quotient_power
                for later in range(variable + 1, VARIABLE_COUNT):
                    output_exponent[later] = exponent[later]
                _add_block(
                    result,
                    (variable,),
                    (),
                    {tuple(output_exponent): matrix},
                )
    return result


def h1_from_h0(h0):
    """Return repository ``H1=-h wedge (-)`` from one H0 column."""

    result = {}
    for input_variable in range(VARIABLE_COUNT):
        for h_variable in range(VARIABLE_COUNT):
            if h_variable == input_variable:
                continue
            polynomial = h0.get(((h_variable,), ()), {})
            if not polynomial:
                continue
            output = tuple(sorted((h_variable, input_variable)))
            sign = -1 if h_variable < input_variable else 1
            _add_block(
                result, output, (input_variable,), polynomial, scale=sign
            )
    return result


def graph_differential(degree, multiplication):
    """Return the repository graph differential in degrees one through three."""

    if degree not in (1, 2, 3):
        raise ValueError("the emitted graph window stops in degree three")
    repository_sign = -1 if degree == 2 else 1
    identity = np.eye(2, dtype=np.int64)
    result = {}
    for source in WEDGES[degree]:
        for position, variable in enumerate(source):
            target = source[:position] + source[position + 1 :]
            sign = repository_sign * (-1 if position % 2 else 1)
            polynomial = {
                _shift(ZERO_EXPONENT, variable, 1): sign * identity % P,
                ZERO_EXPONENT: -sign * multiplication[variable] % P,
            }
            _add_block(result, target, source, polynomial)
    return result


def scalar_identity_block_map(degree, polynomial):
    """Multiplication by a scalar polynomial on ``K_degree``."""

    identity = np.eye(2, dtype=np.int64)
    result = {}
    for wedge in WEDGES[degree]:
        blocks = {}
        for exponent, coefficient in polynomial.items():
            _add_poly_term(
                blocks, exponent, int(coefficient) * identity
            )
        _add_block(result, wedge, wedge, blocks)
    return result


def verify_h_nullhomotopy(polynomial, h0, h1, differentials, label):
    """Verify ``dH+Hd=f`` in graph degrees zero and one."""

    scalar0 = scalar_identity_block_map(0, polynomial)
    assert_block_maps_equal(
        compose_block_maps(differentials[1], h0),
        scalar0,
        f"{label} d1*H0=f",
    )
    scalar1 = scalar_identity_block_map(1, polynomial)
    degree_one = add_block_maps(
        compose_block_maps(differentials[2], h1),
        compose_block_maps(h0, differentials[1]),
    )
    assert_block_maps_equal(
        degree_one,
        scalar1,
        f"{label} d2*H1+H0*d1=f",
    )


def _wedge_sort_sign(values):
    if len(set(values)) != len(values):
        return None, 0
    inversions = sum(
        int(values[left] > values[right])
        for left in range(len(values))
        for right in range(left + 1, len(values))
    )
    return tuple(sorted(values)), -1 if inversions % 2 else 1


def graph_compound(degree, endpoint, changed, chart_multiplication):
    """Return the exact source or target graph map in local jet bases."""

    identity = np.eye(2, dtype=np.int64)
    if degree == 0:
        return {((), ()): {ZERO_EXPONENT: identity}}
    u = 2 * changed
    v = u + 1
    denominator = _shift(ZERO_EXPONENT, v, -1)
    generator_terms = {}
    for variable in range(VARIABLE_COUNT):
        if endpoint == "source" or variable not in (u, v):
            generator_terms[variable] = [
                (variable, ZERO_EXPONENT, 1, ZERO_EXPONENT)
            ]
        elif variable == u:
            operator_u = _shift(ZERO_EXPONENT, u, 1)
            generator_terms[variable] = [
                (u, denominator, 1, ZERO_EXPONENT),
                (v, denominator, -1, operator_u),
            ]
        else:
            operator_v = _shift(ZERO_EXPONENT, v, 1)
            generator_terms[variable] = [
                (v, denominator, -1, operator_v)
            ]

    result = {}
    for source_wedge in WEDGES[degree]:
        for selected in itertools.product(
            *(generator_terms[variable] for variable in source_wedge)
        ):
            output, wedge_sign = _wedge_sort_sign(
                tuple(term[0] for term in selected)
            )
            if not wedge_sign:
                continue
            laurent = ZERO_EXPONENT
            operator_exponent = ZERO_EXPONENT
            scalar = wedge_sign
            for _out, exponent, term_scalar, operator in selected:
                laurent = _sum_exponents(laurent, exponent)
                operator_exponent = _sum_exponents(
                    operator_exponent, operator
                )
                scalar *= int(term_scalar)
            coefficient = evaluate_matrix_polynomial(
                {operator_exponent: 1}, chart_multiplication
            )
            _add_block(
                result,
                output,
                source_wedge,
                {laurent: coefficient},
                scale=scalar,
            )
    return result


def pull_target_exponent(exponent, changed):
    result = list(map(int, exponent))
    u = 2 * changed
    v = u + 1
    u_power = result[u]
    v_power = result[v]
    result[u] = u_power
    result[v] = -(u_power + v_power)
    return tuple(result)


def pull_target_block_map(blocks, changed):
    result = {}
    for key, polynomial in blocks.items():
        pulled = {}
        for exponent, matrix in polynomial.items():
            _add_poly_term(
                pulled, pull_target_exponent(exponent, changed), matrix
            )
        _add_block(result, key[0], key[1], pulled)
    return result


def pull_target_scalar_polynomial(polynomial, changed):
    result = {}
    for exponent, coefficient in polynomial.items():
        target = pull_target_exponent(exponent, changed)
        value = (result.get(target, 0) + int(coefficient)) % P
        if value:
            result[target] = value
        else:
            result.pop(target, None)
    return result


def shift_scalar_polynomial(polynomial, exponent):
    return {
        _sum_exponents(source, exponent): int(coefficient) % P
        for source, coefficient in polynomial.items()
        if int(coefficient) % P
    }


def ordered_laurent_divided(polynomial, variable, multiplication, inverse_variable):
    """Apply the ordered Laurent quotient ``D_variable`` termwise."""

    result = {}
    inverse = inverse_2x2(multiplication[inverse_variable])
    for exponent, coefficient_matrix in polynomial.items():
        exponent = tuple(map(int, exponent))
        power = exponent[variable]
        if not power:
            continue
        prefix = np.eye(2, dtype=np.int64)
        for earlier in range(variable):
            prefix = (
                prefix
                @ matrix_power(
                    multiplication[earlier],
                    exponent[earlier],
                    inverse if earlier == inverse_variable else None,
                )
            ) % P
        if power > 0:
            summands = [
                (power - 1 - r, r, 1) for r in range(power)
            ]
        else:
            summands = [
                (-r, power + r - 1, -1)
                for r in range(1, -power + 1)
            ]
        for polynomial_power, matrix_exponent, sign in summands:
            output_exponent = [0] * VARIABLE_COUNT
            output_exponent[variable] = polynomial_power
            for later in range(variable + 1, VARIABLE_COUNT):
                output_exponent[later] = exponent[later]
            matrix = (
                prefix
                @ matrix_power(
                    multiplication[variable],
                    matrix_exponent,
                    inverse if variable == inverse_variable else None,
                )
                @ coefficient_matrix
            ) % P
            _add_poly_term(
                result, tuple(output_exponent), matrix, scale=sign
            )
    return result


def tensor_ordered_homotopy(blocks, degree, multiplication, inverse_variable):
    """Apply ``h_q(e_I p)=sum_(i<min I) e_i wedge e_I D_i p``."""

    result = {}
    for (output_wedge, input_wedge), polynomial in blocks.items():
        if len(output_wedge) != degree:
            raise HSShamashOverlapError("a homotopy input has the wrong degree")
        for variable in range(min(output_wedge)):
            divided = ordered_laurent_divided(
                polynomial, variable, multiplication, inverse_variable
            )
            _add_block(
                result, (variable, *output_wedge), input_wedge, divided
            )
    return result


def equation_change_exponent(equation, changed):
    if equation == changed:
        power = -3
    elif equation >= 4:
        power = -1
    else:
        power = 0
    return _shift(ZERO_EXPONENT, 2 * changed + 1, power)


def equation_rows(surface, chart, equations, equation_labels):
    rows = []
    for equation, polynomial in enumerate(equations):
        for exponent, coefficient in sorted(polynomial.items()):
            rows.append(
                {
                    "surface": f"Z_{surface}",
                    "chart": chart,
                    "equation_index": equation,
                    "equation_label": equation_labels[equation],
                    "exponents": json.dumps(list(exponent), separators=(",", ":")),
                    "coefficient_mod_31": int(coefficient) % P,
                }
            )
    return rows


def surface_equations(context, surface, chart):
    bits = context["chart_bits"][chart]
    labels = list(
        map(
            int,
            json.loads(context["seeds"][surface]["surface_section_labels"]),
        )
    )
    return [
        *iw_solver.curve_polynomials(bits),
        *[
            iw_solver.theta_polynomial(context["theta"][label], bits)
            for label in labels
        ],
    ], [
        "curve_0",
        "curve_1",
        "curve_2",
        "curve_3",
        f"theta_{labels[0]}",
        f"theta_{labels[1]}",
    ]


def local_support_maps(source_values, source_tangents, target_values, target_tangents, equations_source, equations_target, changed):
    """Construct and verify all local maps on one dual-number support."""

    if len(equations_source) != SURFACE_COUNT or len(equations_target) != SURFACE_COUNT:
        raise HSShamashOverlapError("every surface chart must supply six equations")

    edge_multiplication = local_multiplication(source_values, source_tangents)
    target_multiplication = local_multiplication(target_values, target_tangents)
    inverse_variable = 2 * changed + 1
    inverse_2x2(edge_multiplication[inverse_variable])

    edge_differentials = {
        degree: graph_differential(degree, edge_multiplication)
        for degree in (1, 2, 3)
    }
    target_differentials_raw = {
        degree: graph_differential(degree, target_multiplication)
        for degree in (1, 2, 3)
    }
    target_differentials = {
        degree: pull_target_block_map(
            target_differentials_raw[degree], changed
        )
        for degree in (1, 2, 3)
    }
    graph_maps = {
        degree: graph_compound(
            degree, "target", changed, target_multiplication
        )
        for degree in range(4)
    }
    for degree in (1, 2, 3):
        left = compose_block_maps(edge_differentials[degree], graph_maps[degree])
        right = compose_block_maps(
            graph_maps[degree - 1], target_differentials[degree]
        )
        assert_block_maps_equal(left, right, f"target graph J{degree}")

    gamma_maps = []
    omega_maps = []
    h_rows = []
    nonzero_delta = []
    nonzero_xi = []
    for equation, (source_polynomial, target_polynomial) in enumerate(
        zip(equations_source, equations_target)
    ):
        if np.any(
            evaluate_matrix_polynomial(
                source_polynomial, edge_multiplication
            )
        ):
            raise HSShamashOverlapError(
                "a source equation does not annihilate its local algebra"
            )
        if np.any(
            evaluate_matrix_polynomial(
                target_polynomial, target_multiplication
            )
        ):
            raise HSShamashOverlapError(
                "a target equation does not annihilate its local algebra"
            )

        change = equation_change_exponent(equation, changed)
        pulled_equation = pull_target_scalar_polynomial(
            target_polynomial, changed
        )
        expected_equation = shift_scalar_polynomial(
            source_polynomial, change
        )
        if pulled_equation != expected_equation:
            raise HSShamashOverlapError(
                "the diagonal regular-sequence change failed"
            )

        edge_h0 = divided_difference_matrix(
            source_polynomial, edge_multiplication
        )
        target_h0_raw = divided_difference_matrix(
            target_polynomial, target_multiplication
        )
        target_h0 = pull_target_block_map(target_h0_raw, changed)
        edge_h1 = h1_from_h0(edge_h0)
        target_h1_raw = h1_from_h0(target_h0_raw)
        target_h1 = pull_target_block_map(
            target_h1_raw,
            changed,
        )
        verify_h_nullhomotopy(
            source_polynomial,
            edge_h0,
            edge_h1,
            edge_differentials,
            "edge",
        )
        verify_h_nullhomotopy(
            target_polynomial,
            target_h0_raw,
            target_h1_raw,
            target_differentials_raw,
            "target",
        )
        h_rows.append(
            {
                "equation": equation,
                "edge_H0_sha256": block_map_sha256(edge_h0),
                "edge_H1_sha256": block_map_sha256(edge_h1),
                "target_H0_pulled_sha256": block_map_sha256(target_h0),
                "target_H1_pulled_sha256": block_map_sha256(target_h1),
                "edge_H0_terms": block_map_term_count(edge_h0),
                "edge_H1_terms": block_map_term_count(edge_h1),
                "target_H0_terms": block_map_term_count(target_h0),
                "target_H1_terms": block_map_term_count(target_h1),
            }
        )

        delta = add_block_maps(
            compose_block_maps(graph_maps[1], target_h0),
            shift_block_map(
                compose_block_maps(edge_h0, graph_maps[0]), change
            ),
            right_scale=-1,
        )
        if delta:
            nonzero_delta.append(equation)
        if compose_block_maps(edge_differentials[1], delta):
            raise HSShamashOverlapError("Delta is not a d1-cycle")
        gamma = shift_block_map(
            tensor_ordered_homotopy(
                delta, 1, edge_multiplication, inverse_variable
            ),
            ZERO_EXPONENT,
            scale=-1,
        )
        assert_block_maps_equal(
            compose_block_maps(edge_differentials[2], gamma),
            delta,
            "d2 Gamma=Delta",
        )
        # This is the lower-right equation-sector block of
        # D3_edge*J3=J2*D3_target.  The Laurent entry of B commutes with the
        # graph differential, so there is no derivative-of-B term.
        assert_block_maps_equal(
            compose_block_maps(
                edge_differentials[1],
                shift_block_map(graph_maps[1], change),
            ),
            shift_block_map(
                compose_block_maps(
                    graph_maps[0], target_differentials[1]
                ),
                change,
            ),
            "equation-sector lower D3*J3=J2*D3",
        )

        xi = compose_block_maps(graph_maps[2], target_h1)
        xi = add_block_maps(
            xi,
            compose_block_maps(gamma, target_differentials[1]),
        )
        xi = add_block_maps(
            xi,
            shift_block_map(
                compose_block_maps(edge_h1, graph_maps[1]), change
            ),
            right_scale=-1,
        )
        if xi:
            nonzero_xi.append(equation)
        if compose_block_maps(edge_differentials[2], xi):
            raise HSShamashOverlapError("Xi is not a d2-cycle")
        omega = tensor_ordered_homotopy(
            xi, 2, edge_multiplication, inverse_variable
        )
        assert_block_maps_equal(
            compose_block_maps(edge_differentials[3], omega),
            xi,
            "d3 Omega=Xi",
        )
        gamma_maps.append(gamma)
        omega_maps.append(omega)

    return {
        "graph_maps": graph_maps,
        "gamma": gamma_maps,
        "omega": omega_maps,
        "h_rows": h_rows,
        "nonzero_delta": nonzero_delta,
        "nonzero_xi": nonzero_xi,
        "full_J2_component_checks": SURFACE_COUNT + 1,
        "full_J3_component_checks": 2 * SURFACE_COUNT + 1,
    }


def local_source_support_checks(values, tangents, equations):
    """Verify the strict same-coordinate source endpoint on one support.

    In Hermite jet coordinates the chart-to-edge quotient restriction is the
    identity on a retained support.  Consequently every graph compound is the
    identity, ``B=I_6``, and the correction maps must vanish.  We run this
    calculation rather than recording the zero controls by convention.
    """

    if len(equations) != SURFACE_COUNT:
        raise HSShamashOverlapError("every source chart must supply six equations")
    multiplication = local_multiplication(values, tangents)
    graph_maps = {
        degree: graph_compound(
            degree, "source", 0, multiplication
        )
        for degree in range(4)
    }
    differentials = {
        degree: graph_differential(degree, multiplication)
        for degree in (1, 2, 3)
    }
    for degree in (1, 2, 3):
        assert_block_maps_equal(
            compose_block_maps(differentials[degree], graph_maps[degree]),
            compose_block_maps(
                graph_maps[degree - 1], differentials[degree]
            ),
            f"source graph J{degree}",
        )

    for equation, polynomial in enumerate(equations):
        if np.any(evaluate_matrix_polynomial(polynomial, multiplication)):
            raise HSShamashOverlapError(
                "a source equation does not annihilate its local algebra"
            )
        h0 = divided_difference_matrix(polynomial, multiplication)
        h1 = h1_from_h0(h0)
        verify_h_nullhomotopy(
            polynomial, h0, h1, differentials, "source"
        )
        delta = add_block_maps(
            compose_block_maps(graph_maps[1], h0),
            compose_block_maps(h0, graph_maps[0]),
            right_scale=-1,
        )
        xi = add_block_maps(
            compose_block_maps(graph_maps[2], h1),
            compose_block_maps(h1, graph_maps[1]),
            right_scale=-1,
        )
        if delta or xi:
            raise HSShamashOverlapError(
                f"source equation {equation} did not give strict zero corrections"
            )
    return {
        "graph_maps": graph_maps,
        "gamma_are_zero": True,
        "omega_are_zero": True,
    }


def _support_positions(points):
    return {tuple(map(int, point)): index for index, point in enumerate(points)}


def _verify_rational_support_change(
    source_values,
    source_tangents,
    target_values,
    target_tangents,
    changed,
):
    """Verify values and first jets for ``u'=u/v, v'=1/v``."""

    source_values = np.asarray(source_values, dtype=np.int64) % P
    source_tangents = np.asarray(source_tangents, dtype=np.int64) % P
    target_values = np.asarray(target_values, dtype=np.int64) % P
    target_tangents = np.asarray(target_tangents, dtype=np.int64) % P
    expected_values = source_values.copy()
    expected_tangents = source_tangents.copy()
    u = 2 * int(changed)
    v = u + 1
    denominator = int(source_values[v])
    if not denominator:
        raise HSShamashOverlapError("an overlap contains v=0")
    inverse = pow(denominator, -1, P)
    expected_values[u] = int(source_values[u]) * inverse % P
    expected_values[v] = inverse
    expected_tangents[u] = (
        int(source_tangents[u]) * denominator
        - int(source_values[u]) * int(source_tangents[v])
    ) * inverse * inverse % P
    expected_tangents[v] = (
        -int(source_tangents[v]) * inverse * inverse
    ) % P
    if not np.array_equal(expected_values, target_values):
        raise HSShamashOverlapError("target support values violate u'=u/v,v'=1/v")
    if not np.array_equal(expected_tangents, target_tangents):
        raise HSShamashOverlapError("target support jets violate u'=u/v,v'=1/v")


def load_problem_data():
    """Load and bind the actual chart/edge quotient and Hermite data."""

    context = iw_solver.load_problem_context()
    j1_charts, j1_overlaps, input_paths = j1_lifts.load_inputs()
    charts = {}
    for key, chart in j1_charts.items():
        with np.load(chart["binary"]) as artifact:
            supports = np.asarray(
                artifact["support_point_indices"], dtype=np.int64
            )
            values = np.asarray(artifact["affine_values"], dtype=np.int64) % P
            tangents = np.asarray(
                artifact["affine_tangents"], dtype=np.int64
            ) % P
            hermite = np.asarray(
                artifact["hermite_matrix"], dtype=np.int64
            ) % P
            hermite_inverse = np.asarray(
                artifact["hermite_inverse"], dtype=np.int64
            ) % P
        if values.shape != tangents.shape or values.shape != (
            len(supports),
            VARIABLE_COUNT,
        ):
            raise HSShamashOverlapError("a chart Hermite support array changed")
        identity = np.eye(2 * len(supports), dtype=np.int64)
        if not np.array_equal(hermite @ hermite_inverse % P, identity):
            raise HSShamashOverlapError("a chart Hermite inverse changed")
        charts[key] = {
            **chart,
            "supports": tuple(tuple(map(int, point)) for point in supports),
            "positions": _support_positions(supports),
            "values": values,
            "tangents": tangents,
            "hermite": hermite,
            "hermite_inverse": hermite_inverse,
        }

    overlaps = {}
    selector_checks = 0
    rational_jet_checks = 0
    for key, overlap in j1_overlaps.items():
        surface, _edge = key
        source_chart = charts[(surface, overlap["source"])]
        target_chart = charts[(surface, overlap["target"])]
        with np.load(overlap["binary"]) as artifact:
            supports = np.asarray(
                artifact["support_point_indices"], dtype=np.int64
            )
            hermite = np.asarray(
                artifact["hermite_matrix"], dtype=np.int64
            ) % P
        support_tuples = tuple(tuple(map(int, point)) for point in supports)
        endpoint_positions = {}
        for endpoint, chart in (
            ("source", source_chart),
            ("target", target_chart),
        ):
            try:
                positions = np.asarray(
                    [chart["positions"][point] for point in support_tuples],
                    dtype=np.int64,
                )
            except KeyError as error:
                raise HSShamashOverlapError(
                    f"an overlap support is absent from its {endpoint} chart"
                ) from error
            rows = np.asarray(
                [2 * position + offset for position in positions for offset in (0, 1)],
                dtype=np.int64,
            )
            selected = chart["hermite"][rows, :]
            restriction = overlap[f"{endpoint}_map"]
            if not np.array_equal(hermite @ restriction % P, selected % P):
                raise HSShamashOverlapError(
                    f"the actual {endpoint} restriction is not its Hermite selector"
                )
            endpoint_positions[endpoint] = positions
            selector_checks += 1

        for source_position, target_position in zip(
            endpoint_positions["source"], endpoint_positions["target"]
        ):
            _verify_rational_support_change(
                source_chart["values"][source_position],
                source_chart["tangents"][source_position],
                target_chart["values"][target_position],
                target_chart["tangents"][target_position],
                overlap["changed"],
            )
            rational_jet_checks += 1
        overlaps[key] = {
            **overlap,
            "supports": support_tuples,
            "hermite": hermite,
            "source_positions": endpoint_positions["source"],
            "target_positions": endpoint_positions["target"],
        }

    local_symbols = context["local_symbols"]
    equation_input_paths = [
        Path(getattr(local_symbols, name))
        for name in (
            "GLOBAL_CERTIFICATE",
            "GLOBAL_GENERATORS",
            "ATLAS",
            "ATLAS_EDGES",
            "SEEDS",
            "SUPPORTS",
            "DESCENT",
            "THETA_COEFFICIENTS",
        )
    ]
    j1_manifest = j1_lifts.OUTPUT_DIRECTORY / "manifest.json"
    if not j1_manifest.exists():
        raise HSShamashOverlapError("the exact J1 manifest is missing")
    all_inputs = tuple(
        sorted(set((*input_paths, *equation_input_paths, j1_manifest)))
    )
    return {
        "context": context,
        "charts": charts,
        "overlaps": overlaps,
        "input_paths": all_inputs,
        "selector_checks": selector_checks,
        "rational_jet_checks": rational_jet_checks,
    }


def _serialize_block_map(blocks, degree_out, degree_in, prefix):
    """Serialize a local polynomial block map as compact signed records."""

    rows = []
    for (output_wedge, input_wedge), polynomial in sorted(blocks.items()):
        output_index = WEDGE_INDEX[degree_out][tuple(output_wedge)]
        input_index = WEDGE_INDEX[degree_in][tuple(input_wedge)]
        for exponent, matrix in sorted(polynomial.items()):
            nonzero_rows, nonzero_columns = np.nonzero(
                np.asarray(matrix, dtype=np.int64) % P
            )
            for matrix_row, matrix_column in zip(
                nonzero_rows, nonzero_columns
            ):
                rows.append(
                    [
                        *map(int, prefix),
                        int(output_index),
                        int(input_index),
                        *map(int, exponent),
                        int(matrix_row),
                        int(matrix_column),
                        int(matrix[matrix_row, matrix_column]) % P,
                    ]
                )
    return rows


def deserialize_block_map(records, degree_out, degree_in, prefix_columns):
    """Recover one selected local block map from serialized records."""

    result = {}
    records = np.asarray(records, dtype=np.int64)
    for row in records:
        offset = int(prefix_columns)
        output = WEDGES[degree_out][int(row[offset])]
        input_wedge = WEDGES[degree_in][int(row[offset + 1])]
        exponent = tuple(map(int, row[offset + 2 : offset + 10]))
        matrix_row, matrix_column, coefficient = map(
            int, row[offset + 10 : offset + 13]
        )
        key = (output, input_wedge)
        polynomial = result.setdefault(key, {})
        matrix = polynomial.setdefault(
            exponent, np.zeros((2, 2), dtype=np.int64)
        )
        matrix[matrix_row, matrix_column] = (
            matrix[matrix_row, matrix_column] + coefficient
        ) % P
    return result


def _surface_equation_payload(data):
    rows = []
    equations = {}
    for surface in range(SURFACE_COUNT):
        for chart_name in sorted(data["context"]["chart_bits"]):
            chart_equations, labels = surface_equations(
                data["context"], surface, chart_name
            )
            equations[(surface, chart_name)] = (
                chart_equations,
                labels,
            )
            rows.extend(
                equation_rows(
                    surface, chart_name, chart_equations, labels
                )
            )
    return equations, rows


_WORK_DATA = None
_WORK_EQUATIONS = None
_WORK_OUTPUT_DIRECTORY = None


def _edge_artifact_name(surface, edge):
    return f"Z_{surface}_E{edge:02d}_target_shamash.npz"


def _edge_summary_name(surface, edge):
    return f"Z_{surface}_E{edge:02d}_comparison.json"


def _build_edge_task(task):
    """Build one source/target pair; used directly or by a fork pool."""

    surface, edge = map(int, task)
    data = _WORK_DATA
    equations = _WORK_EQUATIONS
    output_directory = Path(_WORK_OUTPUT_DIRECTORY)
    overlap = data["overlaps"][(surface, edge)]
    source_chart = data["charts"][(surface, overlap["source"])]
    target_chart = data["charts"][(surface, overlap["target"])]
    source_equations, _source_labels = equations[
        (surface, overlap["source"])
    ]
    target_equations, _target_labels = equations[
        (surface, overlap["target"])
    ]

    graph_records = []
    gamma_records = []
    omega_records = []
    verification_counts = []
    verification_hashes = []
    source_zero_checks = 0
    graph_term_count = 0
    gamma_term_count = 0
    omega_term_count = 0
    gamma_scalar_nonzeros = 0
    omega_scalar_nonzeros = 0
    delta_nonzero = Counter()
    xi_nonzero = Counter()

    for support_position, (
        source_position,
        target_position,
    ) in enumerate(
        zip(overlap["source_positions"], overlap["target_positions"])
    ):
        source_values = source_chart["values"][source_position]
        source_tangents = source_chart["tangents"][source_position]
        target_values = target_chart["values"][target_position]
        target_tangents = target_chart["tangents"][target_position]

        source_result = local_source_support_checks(
            source_values, source_tangents, source_equations
        )
        if not (
            source_result["gamma_are_zero"]
            and source_result["omega_are_zero"]
        ):
            raise HSShamashOverlapError("a source zero control changed")
        source_zero_checks += SURFACE_COUNT

        target_result = local_support_maps(
            source_values,
            source_tangents,
            target_values,
            target_tangents,
            source_equations,
            target_equations,
            overlap["changed"],
        )
        for degree, graph_map in target_result["graph_maps"].items():
            graph_term_count += block_map_term_count(graph_map)
            graph_records.extend(
                _serialize_block_map(
                    graph_map,
                    degree,
                    degree,
                    (support_position, degree),
                )
            )
        for equation in range(SURFACE_COUNT):
            gamma = target_result["gamma"][equation]
            omega = target_result["omega"][equation]
            gamma_terms = block_map_term_count(gamma)
            omega_terms = block_map_term_count(omega)
            gamma_nonzeros = block_map_scalar_nonzeros(gamma)
            omega_nonzeros = block_map_scalar_nonzeros(omega)
            gamma_term_count += gamma_terms
            omega_term_count += omega_terms
            gamma_scalar_nonzeros += gamma_nonzeros
            omega_scalar_nonzeros += omega_nonzeros
            gamma_records.extend(
                _serialize_block_map(
                    gamma, 2, 0, (support_position, equation)
                )
            )
            omega_records.extend(
                _serialize_block_map(
                    omega, 3, 1, (support_position, equation)
                )
            )
            delta_flag = int(equation in target_result["nonzero_delta"])
            xi_flag = int(equation in target_result["nonzero_xi"])
            delta_nonzero[equation] += delta_flag
            xi_nonzero[equation] += xi_flag
            verification_counts.append(
                [
                    support_position,
                    equation,
                    delta_flag,
                    xi_flag,
                    gamma_terms,
                    gamma_nonzeros,
                    omega_terms,
                    omega_nonzeros,
                ]
            )
            h_row = target_result["h_rows"][equation]
            verification_hashes.append(
                [
                    h_row["edge_H0_sha256"],
                    h_row["edge_H1_sha256"],
                    h_row["target_H0_pulled_sha256"],
                    h_row["target_H1_pulled_sha256"],
                    block_map_sha256(gamma),
                    block_map_sha256(omega),
                ]
            )

    graph_array = np.asarray(graph_records, dtype=np.int16)
    gamma_array = np.asarray(gamma_records, dtype=np.int16)
    omega_array = np.asarray(omega_records, dtype=np.int16)
    if graph_array.ndim != 2 or graph_array.shape[1] != 15:
        raise HSShamashOverlapError("the graph record schema changed")
    if gamma_array.size:
        gamma_array = gamma_array.reshape((-1, 15))
    else:
        gamma_array = np.zeros((0, 15), dtype=np.int16)
    if omega_array.size:
        omega_array = omega_array.reshape((-1, 15))
    else:
        omega_array = np.zeros((0, 15), dtype=np.int16)
    verification_counts_array = np.asarray(
        verification_counts, dtype=np.int32
    )
    verification_hashes_array = np.asarray(
        verification_hashes, dtype="S64"
    )
    b_exponents = np.asarray(
        [
            equation_change_exponent(equation, overlap["changed"])
            for equation in range(SURFACE_COUNT)
        ],
        dtype=np.int8,
    )
    binary_path = output_directory / _edge_artifact_name(surface, edge)
    temporary = binary_path.with_name(f".{binary_path.name}.{os.getpid()}.tmp.npz")
    np.savez_compressed(
        temporary,
        field=np.asarray([P], dtype=np.uint8),
        surface=np.asarray([surface], dtype=np.uint8),
        edge=np.asarray([edge], dtype=np.uint8),
        changed_factor=np.asarray([overlap["changed"]], dtype=np.uint8),
        support_point_indices=np.asarray(overlap["supports"], dtype=np.uint8),
        source_chart_support_positions=np.asarray(
            overlap["source_positions"], dtype=np.uint8
        ),
        target_chart_support_positions=np.asarray(
            overlap["target_positions"], dtype=np.uint8
        ),
        B_diagonal_laurent_exponents=b_exponents,
        target_graph_J0_J1_J2_J3_terms=graph_array,
        target_Gamma_terms=gamma_array,
        target_Omega_terms=omega_array,
        verification_counts=verification_counts_array,
        verification_hashes=verification_hashes_array,
    )
    temporary.replace(binary_path)

    support_count = len(overlap["supports"])
    summary = {
        "schema": "section12-HS-Shamash-overlap-edge-v1",
        "field": "F_31",
        "surface": f"Z_{surface}",
        "surface_index": surface,
        "edge_id": edge,
        "source_chart": overlap["source"],
        "target_chart": overlap["target"],
        "changed_factor": int(overlap["changed"]),
        "support_points": support_count,
        "endpoint_maps": {
            "source": {
                "J0_J1_J2_J3": "support-selector identities",
                "B": "I_6",
                "Gamma": "strict zero",
                "Omega": "strict zero",
                "equation_support_checks": source_zero_checks,
            },
            "target": {
                "B_diagonal_exponents": b_exponents.tolist(),
                "graph_polynomial_terms": graph_term_count,
                "graph_scalar_records": int(len(graph_array)),
                "Gamma_polynomial_terms": gamma_term_count,
                "Gamma_scalar_records": int(len(gamma_array)),
                "Omega_polynomial_terms": omega_term_count,
                "Omega_scalar_records": int(len(omega_array)),
                "Delta_nonzero_supports_by_equation": {
                    str(key): int(value) for key, value in sorted(delta_nonzero.items())
                },
                "Xi_nonzero_supports_by_equation": {
                    str(key): int(value) for key, value in sorted(xi_nonzero.items())
                },
            },
        },
        "exact_checks": {
            "actual_source_and_target_J0_are_Hermite_support_selectors": True,
            "source_J1_J2_J3_graph_maps_are_exact": True,
            "source_Gamma_and_Omega_are_strict_zero": True,
            "target_J1_J2_J3_graph_chain_identities": True,
            "all_six_source_and_target_equations_annihilate_every_local_factor": True,
            "all_edge_and_target_H0_H1_nullhomotopy_identities": True,
            "pulled_target_regular_sequence_equals_edge_sequence_times_B": True,
            "all_Delta_are_d1_cycles": True,
            "all_d2_Gamma_equals_Delta": True,
            "all_Xi_are_d2_cycles": True,
            "all_d3_Omega_equals_Xi": True,
            "full_D2_J2_equals_J1_D2": True,
            "full_D3_J3_equals_J2_D3": True,
        },
        "hash_roots": {
            "verification_counts_sha256": array_sha256(verification_counts_array),
            "verification_hashes_sha256": array_sha256(
                np.frombuffer(verification_hashes_array.tobytes(), dtype=np.uint8)
            ),
            "graph_terms_sha256": array_sha256(graph_array),
            "Gamma_terms_sha256": array_sha256(gamma_array),
            "Omega_terms_sha256": array_sha256(omega_array),
        },
        "binary_artifact": {
            "path": binary_path.name,
            "bytes": binary_path.stat().st_size,
            "sha256": sha256(binary_path),
        },
    }
    summary_path = output_directory / _edge_summary_name(surface, edge)
    write_json_atomic(summary_path, summary)
    summary["summary_artifact"] = {
        "path": summary_path.name,
        "bytes": summary_path.stat().st_size,
        "sha256": sha256(summary_path),
    }
    return summary


def comparison_program():
    """Return the executable low-degree comparison contract."""

    return {
        "schema": "section12-HS-Shamash-overlap-program-v1",
        "field": "F_31",
        "coordinate_ring": (
            "F_31[x0,...,x7,v_changed^-1] in source-chart coordinates"
        ),
        "target_coordinate_change": [
            "u_target=u_source/v_source",
            "v_target=1/v_source",
        ],
        "graph_differential_convention": (
            "d_q(e_I*a)=sum_r (-1)^(q-1-r) "
            "e_(I without i_r)*(x_i_r-M_i_r)*a"
        ),
        "maps": {
            "JqK": "(Lambda^q D)(I tensor J0), q=0,1,2,3",
            "Delta_j": (
                "J1K*H_target[j,0]-sum_l H_edge[l,0]*B[l,j]*J0"
            ),
            "Gamma_j": "-ordered_tensor_homotopy_degree1(Delta_j)",
            "J2": "[[J2K,Gamma],[0,B tensor J0]]",
            "Xi_j": (
                "J2K*H_target[j,1]+Gamma_j*d1_target-"
                "sum_l H_edge[l,1]*B[l,j]*J1K"
            ),
            "Omega_j": "ordered_tensor_homotopy_degree2(Xi_j)",
            "J3": "[[J3K,Omega],[0,B tensor J1K]]",
        },
        "regular_sequence_change": {
            "equation_order": [
                "curve_0",
                "curve_1",
                "curve_2",
                "curve_3",
                "theta_0",
                "theta_1",
            ],
            "B": (
                "diagonal: changed curve v^-3, other curves 1, "
                "both theta equations v^-1"
            ),
        },
        "ordered_Laurent_contraction": {
            "implemented_by": "tensor_ordered_homotopy",
            "negative_exponents": (
                "q_n(z,A)=(z^n-A^n)/(z-A), including n<0 with A invertible"
            ),
            "signs": {"degree_1": -1, "degree_2": 1},
        },
        "edge_binary_record_columns": {
            "target_graph_J0_J1_J2_J3_terms": [
                "support_position",
                "degree",
                "output_wedge_index",
                "input_wedge_index",
                *[f"x{i}_exponent" for i in range(VARIABLE_COUNT)],
                "local_matrix_row",
                "local_matrix_column",
                "coefficient_mod_31",
            ],
            "target_Gamma_terms": [
                "support_position",
                "equation_index",
                "output_K2_wedge_index",
                "input_K0_wedge_index",
                *[f"x{i}_exponent" for i in range(VARIABLE_COUNT)],
                "local_matrix_row",
                "local_matrix_column",
                "coefficient_mod_31",
            ],
            "target_Omega_terms": [
                "support_position",
                "equation_index",
                "output_K3_wedge_index",
                "input_K1_wedge_index",
                *[f"x{i}_exponent" for i in range(VARIABLE_COUNT)],
                "local_matrix_row",
                "local_matrix_column",
                "coefficient_mod_31",
            ],
            "verification_counts": [
                "support_position",
                "equation_index",
                "Delta_nonzero",
                "Xi_nonzero",
                "Gamma_polynomial_terms",
                "Gamma_scalar_records",
                "Omega_polynomial_terms",
                "Omega_scalar_records",
            ],
        },
        "jet_to_standard_basis_reconstruction": {
            "degree_q_graph": (
                "Jq_standard=(I_C(8,q) tensor H_edge^-1) * "
                "Jq_jet * (I_C(8,q) tensor H_chart)"
            ),
            "degree_2_full": (
                "use diag(I_28 tensor H_edge,I_6 tensor H_edge) and the "
                "corresponding chart matrix"
            ),
            "degree_3_full": (
                "use diag(I_56 tensor H_edge,I_48 tensor H_edge) and the "
                "corresponding chart matrix"
            ),
            "support_selection": (
                "unsupported chart dual-number factors map to zero; retained "
                "factors map by I_2"
            ),
            "bound_inputs": (
                "each edge artifact stores chart support positions; input_hashes.csv "
                "binds every chart/edge Hermite artifact"
            ),
        },
        "interpretation": (
            "The NPZ correction records and this program reconstruct the exact "
            "supportwise/Hermite-quotient J2 and J3 maps with Laurent-polynomial "
            "coefficients.  They are not hashes or zero placeholders."
        ),
    }


def next_coherence_contract():
    """Record exactly what remains after the edgewise comparison."""

    return {
        "schema": "section12-HS-Shamash-higher-descent-next-contract-v1",
        "status": "pending-after-exact-edgewise-J0-through-J3",
        "warning": (
            "Pairwise edge maps are not yet a global descent object.  Pure graph "
            "and diagonal-B blocks are strict; Gamma and Omega can create the "
            "remaining square and cube defects."
        ),
        "totalization_relations": [
            "partial*A1+A1*partial=0",
            "A1^2+partial*A2+A2*partial=0",
            "A1*A2+A2*A1+partial*A3+A3*partial=0",
            "A1*A3+A2^2+A3*A1+partial*A4+A4*partial=0",
        ],
        "cubical_interval_counts_all_six_surfaces": {
            "A1_edges": 1296,
            "A2_squares": 1296,
            "A3_cubes": 576,
            "A4_hypercubes": 96,
        },
        "nontrivial_square_blocks": {
            "degree_2_upper_right": (
                "J2K_AB*Gamma_BD+Gamma_AB*(B_BD tensor J0_BD)-"
                "J2K_AC*Gamma_CD-Gamma_AC*(B_CD tensor J0_CD)"
            ),
            "degree_3_upper_right": (
                "J3K_AB*Omega_BD+Omega_AB*(B_BD tensor J1K_BD)-"
                "J3K_AC*Omega_CD-Omega_AC*(B_CD tensor J1K_CD)"
            ),
            "coordinate_requirement": "compare in common face Laurent coordinates",
        },
        "degree_window_warning": (
            "The P0-through-P3 carrier closes the Ext^2 comparison window.  A "
            "degree-3 square defect may require P4 unless it vanishes."
        ),
        "not_yet_claimed": [
            "square coherence",
            "cube coherence",
            "a global Hartshorne-Serre descent object",
            "transported lambda",
            "a totalized Yoneda value",
            "Question 11.4",
            "the Hodge conjecture",
        ],
    }


def write_artifacts(
    output_directory=OUTPUT_DIRECTORY,
    workers=None,
    selected_edges=None,
    reuse_existing=False,
):
    """Build and write the exact sparse comparison on the requested edges."""

    global _WORK_DATA, _WORK_EQUATIONS, _WORK_OUTPUT_DIRECTORY

    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    data = load_problem_data()
    equations, equation_output_rows = _surface_equation_payload(data)
    tasks = (
        [(surface, edge) for surface in range(SURFACE_COUNT) for edge in range(EDGE_COUNT)]
        if selected_edges is None
        else [tuple(map(int, task)) for task in selected_edges]
    )
    _WORK_DATA = data
    _WORK_EQUATIONS = equations
    _WORK_OUTPUT_DIRECTORY = str(output_directory)

    equation_path = output_directory / "surface_equations_laurent_input.csv"
    input_path = output_directory / "input_hashes.csv"
    program_path = output_directory / "comparison_program.json"
    next_path = output_directory / "higher_descent_next_contract.json"
    write_csv(equation_path, equation_output_rows)
    write_csv(
        input_path,
        [
            {
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in data["input_paths"]
        ],
    )
    write_json_atomic(program_path, comparison_program())
    write_json_atomic(next_path, next_coherence_contract())

    if workers is None:
        workers = min(8, os.cpu_count() or 1)
    workers = max(1, int(workers))
    if reuse_existing:
        summaries = []
        for surface, edge in tasks:
            summary_path = output_directory / _edge_summary_name(surface, edge)
            binary_path = output_directory / _edge_artifact_name(surface, edge)
            if not summary_path.exists() or not binary_path.exists():
                raise HSShamashOverlapError(
                    "metadata-only refresh is missing an edge artifact"
                )
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            if summary.get("schema") != "section12-HS-Shamash-overlap-edge-v1":
                raise HSShamashOverlapError(
                    "metadata-only refresh found an unknown edge schema"
                )
            if summary["binary_artifact"]["sha256"] != sha256(binary_path):
                raise HSShamashOverlapError(
                    "metadata-only refresh found a changed edge binary"
                )
            summary["summary_artifact"] = {
                "path": summary_path.name,
                "bytes": summary_path.stat().st_size,
                "sha256": sha256(summary_path),
            }
            summaries.append(summary)
    elif workers == 1 or len(tasks) == 1:
        summaries = []
        for index, task in enumerate(tasks, 1):
            summaries.append(_build_edge_task(task))
            print(
                f"Shamash edge {index}/{len(tasks)}: Z_{task[0]} E{task[1]:02d}",
                flush=True,
            )
    else:
        fork_context = multiprocessing.get_context("fork")
        with fork_context.Pool(processes=min(workers, len(tasks))) as pool:
            summaries = []
            for index, summary in enumerate(
                pool.imap_unordered(_build_edge_task, tasks, chunksize=1), 1
            ):
                summaries.append(summary)
                if index == 1 or index % 8 == 0 or index == len(tasks):
                    print(
                        f"Shamash edges complete: {index}/{len(tasks)}",
                        flush=True,
                    )
    summaries.sort(key=lambda row: (row["surface_index"], row["edge_id"]))

    registry_rows = []
    for summary in summaries:
        source = summary["endpoint_maps"]["source"]
        target = summary["endpoint_maps"]["target"]
        common = {
            "surface": summary["surface"],
            "edge_id": summary["edge_id"],
            "source_chart": summary["source_chart"],
            "target_chart": summary["target_chart"],
            "changed_factor": summary["changed_factor"],
            "support_points": summary["support_points"],
            "binary_artifact": summary["binary_artifact"]["path"],
            "binary_sha256": summary["binary_artifact"]["sha256"],
            "summary_artifact": summary["summary_artifact"]["path"],
            "summary_sha256": summary["summary_artifact"]["sha256"],
        }
        registry_rows.append(
            {
                **common,
                "endpoint": "source",
                "graph_polynomial_terms": 93 * summary["support_points"],
                "Gamma_polynomial_terms": 0,
                "Omega_polynomial_terms": 0,
                "status": "EXACT_J0_J1_J2_J3_SOURCE_ZERO_CONTROL",
                "equation_support_checks": source["equation_support_checks"],
            }
        )
        registry_rows.append(
            {
                **common,
                "endpoint": "target",
                "graph_polynomial_terms": target["graph_polynomial_terms"],
                "Gamma_polynomial_terms": target["Gamma_polynomial_terms"],
                "Omega_polynomial_terms": target["Omega_polynomial_terms"],
                "status": "EXACT_J0_J1_J2_J3_TARGET_SHAMASH_COMPARISON",
                "equation_support_checks": 6 * summary["support_points"],
            }
        )
    registry_path = output_directory / "endpoint_registry.csv"
    write_csv(registry_path, registry_rows)

    complete = selected_edges is None
    support_points = sum(int(row["support_points"]) for row in summaries)
    target_graph_terms = sum(
        row["endpoint_maps"]["target"]["graph_polynomial_terms"]
        for row in summaries
    )
    gamma_terms = sum(
        row["endpoint_maps"]["target"]["Gamma_polynomial_terms"]
        for row in summaries
    )
    omega_terms = sum(
        row["endpoint_maps"]["target"]["Omega_polynomial_terms"]
        for row in summaries
    )
    gamma_records = sum(
        row["endpoint_maps"]["target"]["Gamma_scalar_records"]
        for row in summaries
    )
    omega_records = sum(
        row["endpoint_maps"]["target"]["Omega_scalar_records"]
        for row in summaries
    )
    certificate = {
        "schema": "section12-HS-Shamash-overlap-lifts-certificate-v1",
        "field": "F_31",
        "status": (
            "all-384-endpoint-supportwise-Hermite-J0-through-J3-comparisons-exact"
            if complete
            else "selected-edge-debug-build"
        ),
        "counts": {
            "surfaces": SURFACE_COUNT,
            "edges": len(summaries),
            "endpoint_maps": 2 * len(summaries),
            "source_endpoint_maps": len(summaries),
            "target_endpoint_maps": len(summaries),
            "overlap_support_points": support_points,
            "source_equation_support_zero_checks": SURFACE_COUNT * support_points,
            "target_equation_support_checks": SURFACE_COUNT * support_points,
            "Hermite_J0_selector_checks": data["selector_checks"],
            "rational_value_and_jet_transition_checks": data["rational_jet_checks"],
            "target_graph_polynomial_terms": target_graph_terms,
            "target_Gamma_polynomial_terms": gamma_terms,
            "target_Omega_polynomial_terms": omega_terms,
            "target_Gamma_scalar_records": gamma_records,
            "target_Omega_scalar_records": omega_records,
            "full_J2_component_identities": 7 * support_points,
            "full_J3_component_identities": 13 * support_points,
            "surface_equation_Laurent_terms": len(equation_output_rows),
            "input_files_bound": len(data["input_paths"]),
        },
        "exact_checks": {
            "actual_standard_basis_J0_maps_recover_support_selectors_in_Hermite_bases": True,
            "all_source_endpoint_corrections_are_computed_strict_zero": True,
            "all_target_J1_J2_J3_graph_compounds_are_exact_chain_maps": True,
            "all_six_equations_and_diagonal_B_changes_are_exact": True,
            "all_edge_source_target_H0_H1_identities_are_exact": True,
            "all_Delta_and_Xi_cycles_are_exact": True,
            "all_Gamma_and_Omega_homotopy_equations_are_exact": True,
            "all_full_D2_edge_J2_equals_J1_D2_target_component_ledgers_close": True,
            "all_full_D3_edge_J3_equals_J2_D3_target_component_ledgers_close": True,
            "all_serialized_correction_records_are_materialized_not_hash_placeholders": True,
        },
        "artifacts": {
            "surface_equations": equation_path.name,
            "input_hashes": input_path.name,
            "comparison_program": program_path.name,
            "endpoint_registry": registry_path.name,
            "per_edge_binary_artifacts": len(summaries),
            "per_edge_verification_summaries": len(summaries),
            "higher_descent_next_contract": next_path.name,
        },
        "claim_boundary": {
            "proved": (
                "the exact supportwise/Hermite-quotient J0,J1,J2,J3 Shamash "
                "comparison with Laurent-polynomial coefficients on both endpoints "
                "of each emitted original atlas edge"
            ),
            "not_claimed": (
                "square/cube coherence, a global descent object, transported "
                "lambda, a totalized Yoneda value, Question 11.4, or Hodge"
            ),
        },
    }
    certificate_path = output_directory / "shamash_overlap_lifts_certificate.json"
    write_json_atomic(certificate_path, certificate)

    artifact_paths = [
        equation_path,
        input_path,
        program_path,
        next_path,
        registry_path,
        certificate_path,
    ]
    for summary in summaries:
        artifact_paths.extend(
            (
                output_directory / summary["binary_artifact"]["path"],
                output_directory / summary["summary_artifact"]["path"],
            )
        )
    manifest = {
        "schema": "section12-HS-Shamash-overlap-lifts-manifest-v1",
        "generator": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "files": {
            path.name: {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in artifact_paths
        },
    }
    write_json_atomic(output_directory / "manifest.json", manifest)
    return certificate, summaries


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-directory", type=Path, default=OUTPUT_DIRECTORY
    )
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--surface", type=int)
    parser.add_argument("--edge", type=int)
    args = parser.parse_args()
    if (args.surface is None) != (args.edge is None):
        parser.error("--surface and --edge must be supplied together")
    return args


def main():
    args = parse_args()
    selected = (
        None
        if args.surface is None
        else [(int(args.surface), int(args.edge))]
    )
    certificate, _summaries = write_artifacts(
        args.output_directory,
        workers=args.workers,
        selected_edges=selected,
        reuse_existing=args.metadata_only,
    )
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
