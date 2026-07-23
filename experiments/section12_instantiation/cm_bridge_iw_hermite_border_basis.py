"""Build exact Hermite border bases for the Section 12 curvilinear schemes.

This is a finite-linear-algebra alternative to constructing ``I_W`` by an
81-fold ideal intersection.  On every affine product chart the rational
support points and their invariant tangent directions define the algebra

    O(W cap U) = product_p F_31[epsilon]/(epsilon^2).

Value and directional-derivative evaluation selects a deterministic order
ideal of monomials.  The induced variable multiplication matrices certify a
border basis: they commute, are cyclic from 1, satisfy the six ambient/surface
equations, and reproduce every value-and-derivative functional exactly.

The calculation also evaluates all 162 translated normalized sections.  It
proves that they lie in ``I_W`` and, together with the existing conormal-rank
certificate, that their stalks equal ``I_W`` at every point of ``W``.  It does
not by itself exclude an additional common zero away from ``W``.  Consequently
it replaces the expensive construction of the *target* ideal ``I_W`` and gives
normal-form multiplication, but it does not replace the separate saturation
check ``J=I_W`` or a Schreyer presentation of ``J``.

No claim about Question 11.4 or the Hodge conjecture is made here.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_iw_presentation_solver as iw_solver,
)


P = 31
VARIABLE_COUNT = 8
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_iw_hermite_border_basis"
)


class HermiteBorderBasisError(RuntimeError):
    """An exact input or border-basis invariant failed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank_mod(matrix: np.ndarray) -> int:
    work = np.asarray(matrix, dtype=np.int64).copy() % P
    row = 0
    for column in range(work.shape[1]):
        candidates = np.flatnonzero(work[row:, column])
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        if pivot != row:
            work[[row, pivot]] = work[[pivot, row]]
        work[row] = work[row] * pow(int(work[row, column]), -1, P) % P
        nonzero = np.flatnonzero(work[:, column])
        nonzero = nonzero[nonzero != row]
        if len(nonzero):
            multipliers = work[nonzero, column].copy()
            work[nonzero] = (
                work[nonzero] - multipliers[:, None] * work[row]
            ) % P
        row += 1
        if row == work.shape[0]:
            break
    return row


def inverse_mod(matrix: np.ndarray) -> np.ndarray:
    source = np.asarray(matrix, dtype=np.int64) % P
    if source.shape[0] != source.shape[1]:
        raise ValueError("only square matrices can be inverted")
    size = source.shape[0]
    work = np.concatenate((source.copy(), np.eye(size, dtype=np.int64)), axis=1)
    for column in range(size):
        candidates = np.flatnonzero(work[column:, column])
        if not len(candidates):
            raise HermiteBorderBasisError("the selected Hermite basis is singular")
        pivot = column + int(candidates[0])
        if pivot != column:
            work[[column, pivot]] = work[[pivot, column]]
        work[column] = work[column] * pow(int(work[column, column]), -1, P) % P
        nonzero = np.flatnonzero(work[:, column])
        nonzero = nonzero[nonzero != column]
        if len(nonzero):
            multipliers = work[nonzero, column].copy()
            work[nonzero] = (
                work[nonzero] - multipliers[:, None] * work[column]
            ) % P
    inverse = work[:, size:]
    if not np.array_equal(source @ inverse % P, np.eye(size, dtype=np.int64)):
        raise HermiteBorderBasisError("the modular inverse failed substitution")
    return inverse


def monomials_of_degree(degree: int):
    """Yield eight-variable exponent tuples in deterministic graded lex order."""

    def visit(remaining, slots, prefix):
        if slots == 1:
            yield tuple((*prefix, remaining))
            return
        for exponent in range(remaining + 1):
            yield from visit(remaining - exponent, slots - 1, (*prefix, exponent))

    yield from visit(degree, VARIABLE_COUNT, ())


def affine_point_and_tangent(point, tangent_coefficients, bits, point_values, point_first):
    values = []
    derivatives = []
    for factor, (point_index, coefficient, bit) in enumerate(
        zip(point, tangent_coefficients, bits)
    ):
        z, x, y = map(int, point_values[point_index])
        dz, dx, dy = map(int, point_first[point_index])
        coefficient = int(coefficient) % P
        dz, dx, dy = (
            coefficient * dz % P,
            coefficient * dx % P,
            coefficient * dy % P,
        )
        denominator = y if bit else z
        if not denominator:
            raise HermiteBorderBasisError(
                f"factor {factor} lies outside the selected affine chart"
            )
        inverse = pow(denominator, -1, P)
        inverse_squared = inverse * inverse % P
        if bit:
            values.extend((x * inverse % P, z * inverse % P))
            derivatives.extend(
                (
                    (dx * y - x * dy) * inverse_squared % P,
                    (dz * y - z * dy) * inverse_squared % P,
                )
            )
        else:
            values.extend((x * inverse % P, y * inverse % P))
            derivatives.extend(
                (
                    (dx * z - x * dz) * inverse_squared % P,
                    (dy * z - y * dz) * inverse_squared % P,
                )
            )
    return tuple(values), tuple(derivatives)


def monomial_jet(exponent, values, derivatives):
    value = 1
    for coordinate, power in zip(values, exponent):
        value = value * pow(int(coordinate), int(power), P) % P
    derivative = 0
    for direction, power in enumerate(exponent):
        if not power or not derivatives[direction]:
            continue
        term = int(power) * int(derivatives[direction])
        for coordinate_index, coordinate_power in enumerate(exponent):
            reduced_power = coordinate_power - int(coordinate_index == direction)
            term = term * pow(int(values[coordinate_index]), reduced_power, P) % P
        derivative += term
    return value, derivative % P


def polynomial_jet(polynomial, values, derivatives):
    value = 0
    derivative = 0
    for exponent, coefficient in polynomial.items():
        monomial_value, monomial_derivative = monomial_jet(
            exponent, values, derivatives
        )
        value += int(coefficient) * monomial_value
        derivative += int(coefficient) * monomial_derivative
    return value % P, derivative % P


def hermite_column(exponent, affine_points):
    result = []
    for values, derivatives in affine_points:
        result.extend(monomial_jet(exponent, values, derivatives))
    return np.asarray(result, dtype=np.int64) % P


def select_standard_monomials(affine_points, maximum_degree=16):
    dimension = 2 * len(affine_points)
    echelon = {}
    basis = []
    columns = []
    tested = 0
    final_degree = None
    for degree in range(maximum_degree + 1):
        for exponent in monomials_of_degree(degree):
            tested += 1
            column = hermite_column(exponent, affine_points)
            reduced = column.copy()
            for pivot, pivot_vector in echelon.items():
                if reduced[pivot]:
                    reduced = (reduced - reduced[pivot] * pivot_vector) % P
            nonzero = np.flatnonzero(reduced)
            if not len(nonzero):
                continue
            pivot = int(nonzero[0])
            reduced = reduced * pow(int(reduced[pivot]), -1, P) % P
            echelon[pivot] = reduced
            basis.append(exponent)
            columns.append(column)
            if len(basis) == dimension:
                final_degree = degree
                break
        if final_degree is not None:
            break
    if len(basis) != dimension:
        raise HermiteBorderBasisError(
            f"Hermite rank stopped at {len(basis)} of {dimension} by degree {maximum_degree}"
        )
    basis_set = set(basis)
    for exponent in basis:
        for variable, power in enumerate(exponent):
            if power:
                divisor = list(exponent)
                divisor[variable] -= 1
                if tuple(divisor) not in basis_set:
                    raise HermiteBorderBasisError("the selected monomials are not an order ideal")
    matrix = np.stack(columns, axis=1) % P
    return basis, matrix, tested, final_degree


def border_multiplication(basis, hermite_matrix, affine_points):
    dimension = len(basis)
    inverse = inverse_mod(hermite_matrix)
    basis_index = {exponent: index for index, exponent in enumerate(basis)}
    multiplication = np.zeros(
        (VARIABLE_COUNT, dimension, dimension), dtype=np.int64
    )
    border_normal_forms = {}
    for variable in range(VARIABLE_COUNT):
        for source, exponent in enumerate(basis):
            product = list(exponent)
            product[variable] += 1
            product = tuple(product)
            if product in basis_index:
                coefficients = np.zeros(dimension, dtype=np.int64)
                coefficients[basis_index[product]] = 1
            else:
                coefficients = inverse @ hermite_column(product, affine_points) % P
                border_normal_forms.setdefault(product, coefficients)
            multiplication[variable, :, source] = coefficients
    for left in range(VARIABLE_COUNT):
        for right in range(left + 1, VARIABLE_COUNT):
            if np.any(
                (multiplication[left] @ multiplication[right]
                 - multiplication[right] @ multiplication[left])
                % P
            ):
                raise HermiteBorderBasisError("two border multiplication matrices do not commute")
    one = basis_index.get((0,) * VARIABLE_COUNT)
    if one is None:
        raise HermiteBorderBasisError("the order ideal does not contain one")
    unit = np.zeros(dimension, dtype=np.int64)
    unit[one] = 1
    for target, exponent in enumerate(basis):
        vector = unit.copy()
        for variable, power in enumerate(exponent):
            for _ in range(power):
                vector = multiplication[variable] @ vector % P
        expected = np.zeros(dimension, dtype=np.int64)
        expected[target] = 1
        if not np.array_equal(vector, expected):
            raise HermiteBorderBasisError("the border algebra is not cyclic on its order ideal")
    return multiplication, border_normal_forms


def evaluate_polynomials(polynomials, affine_points):
    """Evaluate a polynomial batch after caching every distinct monomial jet."""

    exponents = sorted({exponent for polynomial in polynomials for exponent in polynomial})
    exponent_index = {exponent: index for index, exponent in enumerate(exponents)}
    coefficients = np.zeros((len(polynomials), len(exponents)), dtype=np.int64)
    for row, polynomial in enumerate(polynomials):
        for exponent, coefficient in polynomial.items():
            coefficients[row, exponent_index[exponent]] = int(coefficient) % P
    jets = np.stack(
        [hermite_column(exponent, affine_points) for exponent in exponents], axis=0
    )
    return coefficients @ jets % P


def prepare_border_context(context):
    """Cache the group action and all 972 translated coefficient vectors once."""

    if "_border_records" in context:
        return context
    (
        _global_certificate,
        descent,
        _seeds,
        _support_rows,
        _atlas,
        edges,
        _theta,
        all_generators,
    ) = context["local_symbols"].load_inputs()
    points, point_values, point_first = context["local_symbols"].seed_emitter.curve_data()
    records, _element_index, group_checks = context["local_symbols"].group_data(
        descent, points
    )
    translated = {}
    for surface in range(6):
        translated[surface] = [
            context["local_symbols"].apply_degree2_action(
                all_generators[(surface, generator)], record
            )
            for record in records
            for generator in ("X", "Q")
        ]
    context.update(
        {
            "_border_points": points,
            "_border_point_values": point_values,
            "_border_point_first": point_first,
            "_border_records": records,
            "_border_group_checks": group_checks,
            "_border_translated_original": translated,
            "_border_edges": edges,
        }
    )
    return context


def chart_border_basis(context, surface: int, chart: str):
    prepare_border_context(context)
    bits = context["chart_bits"][chart]
    tangent = json.loads(context["seeds"][surface]["tangent_fixed_frame_coefficients"])
    points = context["_border_points"]
    point_values = context["_border_point_values"]
    point_first = context["_border_point_first"]
    supports = [
        point
        for point in context["support_by_surface"][surface]
        if context["local_symbols"].chart_membership(point, bits, point_values)
    ]
    supports.sort()
    affine_points = [
        affine_point_and_tangent(point, tangent, bits, point_values, point_first)
        for point in supports
    ]
    basis, hermite, tested, final_degree = select_standard_monomials(affine_points)
    multiplication, border = border_multiplication(basis, hermite, affine_points)

    curve = iw_solver.curve_polynomials(bits)
    labels = list(
        map(int, json.loads(context["seeds"][surface]["surface_section_labels"]))
    )
    surface_polynomials = [
        iw_solver.theta_polynomial(context["theta"][label], bits) for label in labels
    ]
    ambient_jets = evaluate_polynomials([*curve, *surface_polynomials], affine_points)
    if np.any(ambient_jets):
        raise HermiteBorderBasisError("an ambient or surface equation has a nonzero Hermite jet")

    original_generators = []
    for vector in context["_border_translated_original"][surface]:
        original_generators.append(iw_solver.degree2_polynomial(vector, bits))
    generator_jets = evaluate_polynomials(original_generators, affine_points)
    if np.any(generator_jets):
        raise HermiteBorderBasisError("a translated normalized generator does not vanish on W")

    hermite_inverse = inverse_mod(hermite)
    if not np.array_equal(hermite_inverse @ hermite % P, np.eye(len(basis), dtype=np.int64)):
        raise HermiteBorderBasisError("the two-sided Hermite inverse check failed")
    point_idempotents = hermite_inverse[:, 0::2]
    point_nilpotents = hermite_inverse[:, 1::2]
    one_index = basis.index((0,) * VARIABLE_COUNT)
    expected_one = np.zeros(len(basis), dtype=np.int64)
    expected_one[one_index] = 1
    if not np.array_equal(point_idempotents.sum(axis=1) % P, expected_one):
        raise HermiteBorderBasisError("the point idempotents do not sum to one")
    dual_coordinates = np.empty_like(hermite_inverse)
    dual_coordinates[:, 0::2] = point_idempotents
    dual_coordinates[:, 1::2] = point_nilpotents
    if not np.array_equal(hermite @ dual_coordinates % P, np.eye(len(basis), dtype=np.int64)):
        raise HermiteBorderBasisError("the pointwise dual-number coordinates changed")
    return {
        "surface": f"Z_{surface}",
        "surface_index": surface,
        "chart": chart,
        "chart_bits": bits,
        "support_point_indices": supports,
        "affine_points": affine_points,
        "basis": basis,
        "hermite_matrix": hermite,
        "hermite_inverse": hermite_inverse,
        "point_idempotents": point_idempotents,
        "point_nilpotents": point_nilpotents,
        "multiplication": multiplication,
        "border_normal_forms": border,
        "summary": {
            "support_points": len(supports),
            "quotient_length": len(basis),
            "selected_standard_monomials": len(basis),
            "highest_selected_degree": final_degree,
            "monomials_tested": tested,
            "border_monomials": len(border),
            "ambient_equation_jet_rank": rank_mod(ambient_jets),
            "translated_normalized_generators": len(original_generators),
            "translated_generator_hermite_rank": rank_mod(generator_jets),
            "commuting_variable_pairs_checked": 28,
            "cyclic_basis_vectors_checked": len(basis),
            "point_idempotents": len(supports),
            "point_nilpotents": len(supports),
        },
    }


def coordinate_multiplication_from_vector(coefficients, basis, multiplication):
    """Return multiplication by a quotient element in the standard basis."""

    dimension = len(basis)
    result = np.zeros((dimension, dimension), dtype=np.int64)
    identity = np.eye(dimension, dtype=np.int64)
    for coefficient, exponent in zip(coefficients, basis):
        if not coefficient:
            continue
        operator = identity
        for variable, power in enumerate(exponent):
            for _ in range(power):
                operator = multiplication[variable] @ operator % P
        result = (result + int(coefficient) * operator) % P
    return result


def overlap_restriction(context, surface, edge, chart_results):
    """Build exact quotient-sheaf restriction maps on one cube edge."""

    source = edge["source_chart"]
    target = edge["target_chart"]
    source_result = chart_results[(surface, source)]
    target_result = chart_results[(surface, target)]
    source_bits = context["chart_bits"][source]
    target_bits = context["chart_bits"][target]
    tangent = json.loads(context["seeds"][surface]["tangent_fixed_frame_coefficients"])
    point_values = context["_border_point_values"]
    point_first = context["_border_point_first"]
    source_support = set(source_result["support_point_indices"])
    target_support = set(target_result["support_point_indices"])
    supports = sorted(source_support & target_support)
    if not supports:
        raise HermiteBorderBasisError("a cube edge has no support points on its overlap")
    source_affine = [
        affine_point_and_tangent(point, tangent, source_bits, point_values, point_first)
        for point in supports
    ]
    target_affine = [
        affine_point_and_tangent(point, tangent, target_bits, point_values, point_first)
        for point in supports
    ]
    edge_basis, edge_hermite, tested, final_degree = select_standard_monomials(
        source_affine
    )
    edge_multiplication, edge_border = border_multiplication(
        edge_basis, edge_hermite, source_affine
    )
    edge_inverse = inverse_mod(edge_hermite)

    source_columns = np.stack(
        [hermite_column(exponent, source_affine) for exponent in source_result["basis"]],
        axis=1,
    )
    target_columns = np.stack(
        [hermite_column(exponent, target_affine) for exponent in target_result["basis"]],
        axis=1,
    )
    source_map = edge_inverse @ source_columns % P
    target_map = edge_inverse @ target_columns % P
    edge_dimension = len(edge_basis)
    if rank_mod(source_map) != edge_dimension or rank_mod(target_map) != edge_dimension:
        raise HermiteBorderBasisError("a chart-to-overlap restriction is not surjective")

    for variable in range(VARIABLE_COUNT):
        if np.any(
            (
                source_map @ source_result["multiplication"][variable]
                - edge_multiplication[variable] @ source_map
            )
            % P
        ):
            raise HermiteBorderBasisError("the source restriction is not an algebra map")

        target_coordinate_jet = hermite_column(
            tuple(int(index == variable) for index in range(VARIABLE_COUNT)),
            target_affine,
        )
        target_coordinate = edge_inverse @ target_coordinate_jet % P
        target_operator = coordinate_multiplication_from_vector(
            target_coordinate, edge_basis, edge_multiplication
        )
        if np.any(
            (
                target_map @ target_result["multiplication"][variable]
                - target_operator @ target_map
            )
            % P
        ):
            raise HermiteBorderBasisError("the target restriction is not an algebra map")
    return {
        "surface": f"Z_{surface}",
        "surface_index": surface,
        "edge_id": int(edge["edge_id"]),
        "source_chart": source,
        "target_chart": target,
        "changed_factor": int(edge["changed_factor"]),
        "support_point_indices": supports,
        "basis": edge_basis,
        "hermite_matrix": edge_hermite,
        "multiplication": edge_multiplication,
        "border_normal_forms": edge_border,
        "source_restriction": source_map,
        "target_restriction": target_map,
        "summary": {
            "overlap_support_points": len(supports),
            "overlap_quotient_length": edge_dimension,
            "highest_selected_degree": final_degree,
            "monomials_tested": tested,
            "source_quotient_length": len(source_result["basis"]),
            "target_quotient_length": len(target_result["basis"]),
            "source_restriction_rank": rank_mod(source_map),
            "target_restriction_rank": rank_mod(target_map),
            "source_kernel_dimension": len(source_result["basis"]) - edge_dimension,
            "target_kernel_dimension": len(target_result["basis"]) - edge_dimension,
            "coordinate_intertwining_checks": 16,
        },
    }


def sparse_entries(matrix):
    rows, columns = np.nonzero(np.asarray(matrix) % P)
    return np.asarray(
        [[int(row), int(column), int(matrix[row, column] % P)] for row, column in zip(rows, columns)],
        dtype=np.int64,
    )


def write_chart_artifact(result, output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    stem = f"{result['surface']}_{result['chart']}"
    path = output_directory / f"{stem}.npz"
    border_exponents = sorted(result["border_normal_forms"])
    border_coefficients = np.stack(
        [result["border_normal_forms"][exponent] for exponent in border_exponents], axis=0
    )
    np.savez_compressed(
        path,
        field=np.asarray([P], dtype=np.uint8),
        chart_bits=np.asarray(result["chart_bits"], dtype=np.uint8),
        support_point_indices=np.asarray(result["support_point_indices"], dtype=np.uint8),
        affine_values=np.asarray([row[0] for row in result["affine_points"]], dtype=np.uint8),
        affine_tangents=np.asarray([row[1] for row in result["affine_points"]], dtype=np.uint8),
        standard_monomials=np.asarray(result["basis"], dtype=np.uint8),
        hermite_matrix=np.asarray(result["hermite_matrix"], dtype=np.uint8),
        hermite_inverse=np.asarray(result["hermite_inverse"], dtype=np.uint8),
        point_idempotent_coefficients=np.asarray(result["point_idempotents"], dtype=np.uint8),
        point_nilpotent_coefficients=np.asarray(result["point_nilpotents"], dtype=np.uint8),
        multiplication_matrices=np.asarray(result["multiplication"], dtype=np.uint8),
        border_monomials=np.asarray(border_exponents, dtype=np.uint8),
        border_normal_form_coefficients=np.asarray(border_coefficients, dtype=np.uint8),
    )
    summary_path = output_directory / f"{stem}.json"
    summary_path.write_text(
        json.dumps(
            {
                "schema": "section12-IW-Hermite-border-basis-v1",
                "field": "F_31",
                "surface": result["surface"],
                "chart": result["chart"],
                "summary": result["summary"],
                "binary_artifact": {"path": path.name, "sha256": sha256(path)},
                "exact_checks": {
                    "Hermite_matrix_invertible": True,
                    "standard_monomials_form_order_ideal": True,
                    "multiplication_matrices_commute": True,
                    "order_ideal_is_cyclic_from_one": True,
                    "point_idempotents_sum_to_one_and_dual_numbers_are_interpolated": True,
                    "four_curve_and_two_surface_equations_vanish_to_first_order": True,
                    "all_162_translated_normalized_sections_vanish_to_first_order": True,
                },
                "claim_boundary": {
                    "proved": (
                        "the border relations present the exact curvilinear ideal I_W "
                        "on this chart; the 162 translated normalized sections lie in I_W"
                    ),
                    "not_proved": (
                        "the translated-section ideal J has no additional base locus; "
                        "J=I_W still needs saturation/length or explicit reverse ideal membership"
                    ),
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return path, summary_path


def write_overlap_artifact(result, output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory) / "overlaps"
    output_directory.mkdir(parents=True, exist_ok=True)
    stem = f"{result['surface']}_E{result['edge_id']:02d}"
    path = output_directory / f"{stem}.npz"
    border_exponents = sorted(result["border_normal_forms"])
    border_coefficients = np.stack(
        [result["border_normal_forms"][exponent] for exponent in border_exponents], axis=0
    )
    np.savez_compressed(
        path,
        field=np.asarray([P], dtype=np.uint8),
        support_point_indices=np.asarray(result["support_point_indices"], dtype=np.uint8),
        standard_monomials=np.asarray(result["basis"], dtype=np.uint8),
        hermite_matrix=np.asarray(result["hermite_matrix"], dtype=np.uint8),
        multiplication_matrices=np.asarray(result["multiplication"], dtype=np.uint8),
        border_monomials=np.asarray(border_exponents, dtype=np.uint8),
        border_normal_form_coefficients=np.asarray(border_coefficients, dtype=np.uint8),
        source_restriction=np.asarray(result["source_restriction"], dtype=np.uint8),
        target_restriction=np.asarray(result["target_restriction"], dtype=np.uint8),
    )
    summary_path = output_directory / f"{stem}.json"
    summary_path.write_text(
        json.dumps(
            {
                "schema": "section12-IW-Hermite-overlap-restriction-v1",
                "field": "F_31",
                "surface": result["surface"],
                "edge_id": result["edge_id"],
                "source_chart": result["source_chart"],
                "target_chart": result["target_chart"],
                "changed_factor": result["changed_factor"],
                "summary": result["summary"],
                "binary_artifact": {"path": path.name, "sha256": sha256(path)},
                "exact_checks": {
                    "both_restrictions_surjective": True,
                    "source_coordinate_multiplications_intertwine": True,
                    "target_coordinate_multiplications_intertwine_after_rational_chart_change": True,
                    "overlap_multiplication_matrices_commute": True,
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return path, summary_path


def run_all(context, output_directory=OUTPUT_DIRECTORY, include_overlaps=True):
    output_directory = Path(output_directory)
    chart_results = {}
    chart_rows = []
    for surface in range(6):
        for chart in sorted(context["chart_bits"]):
            result = chart_border_basis(context, surface, chart)
            path, summary_path = write_chart_artifact(result, output_directory)
            chart_results[(surface, chart)] = result
            chart_rows.append(
                {
                    "surface": result["surface"],
                    "chart": chart,
                    **result["summary"],
                    "artifact": path.name,
                    "artifact_sha256": sha256(path),
                    "summary_artifact": summary_path.name,
                }
            )
            print(
                f"{result['surface']} {chart}: length {len(result['basis'])}, "
                f"degree {result['summary']['highest_selected_degree']}",
                flush=True,
            )

    overlap_rows = []
    if include_overlaps:
        for surface in range(6):
            for edge in context["_border_edges"]:
                result = overlap_restriction(context, surface, edge, chart_results)
                path, summary_path = write_overlap_artifact(result, output_directory)
                overlap_rows.append(
                    {
                        "surface": result["surface"],
                        "edge_id": result["edge_id"],
                        "source_chart": result["source_chart"],
                        "target_chart": result["target_chart"],
                        **result["summary"],
                        "artifact": str(path.relative_to(output_directory)),
                        "artifact_sha256": sha256(path),
                        "summary_artifact": str(summary_path.relative_to(output_directory)),
                    }
                )
                print(
                    f"{result['surface']} E{result['edge_id']:02d}: overlap length "
                    f"{result['summary']['overlap_quotient_length']}",
                    flush=True,
                )

    manifest = output_directory / "hermite_border_basis_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "section12-IW-Hermite-border-basis-manifest-v1",
                "field": "F_31",
                "chart_count": len(chart_rows),
                "overlap_count": len(overlap_rows),
                "charts": chart_rows,
                "overlaps": overlap_rows,
                "exact_scope": {
                    "I_W_target_presentations": "complete on all 96 affine surface charts",
                    "O_W_quotient_restrictions": (
                        "complete on all 192 oriented surface-edge instances"
                        if include_overlaps
                        else "not requested in this run"
                    ),
                    "normalized_J_equality": (
                        "not inferred: zero Hermite jets prove J subset I_W, but an "
                        "extra base locus remains a separate saturation question"
                    ),
                    "free_I_W_Schreyer_matrices": "not emitted by this border-basis tranche",
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest


def verify_manifest(path=OUTPUT_DIRECTORY / "hermite_border_basis_manifest.json"):
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    charts = data.get("charts", [])
    overlaps = data.get("overlaps", [])
    chart_keys = {(row["surface"], row["chart"]) for row in charts}
    expected_charts = {
        (f"Z_{surface}", f"C{chart:04b}")
        for surface in range(6)
        for chart in range(16)
    }
    if len(charts) != 96 or chart_keys != expected_charts:
        raise HermiteBorderBasisError("the manifest does not contain the 96 chart jobs")
    overlap_keys = {(row["surface"], int(row["edge_id"])) for row in overlaps}
    expected_overlaps = {
        (f"Z_{surface}", edge) for surface in range(6) for edge in range(32)
    }
    if len(overlaps) != 192 or overlap_keys != expected_overlaps:
        raise HermiteBorderBasisError("the manifest does not contain the 192 overlap jobs")

    for row in charts:
        artifact = path.parent / row["artifact"]
        if sha256(artifact) != row["artifact_sha256"]:
            raise HermiteBorderBasisError("a chart border-basis artifact hash changed")
        with np.load(artifact) as stored:
            dimension = int(row["quotient_length"])
            support_count = int(row["support_points"])
            if stored["standard_monomials"].shape != (dimension, VARIABLE_COUNT):
                raise HermiteBorderBasisError("a chart standard-monomial shape changed")
            if stored["hermite_matrix"].shape != (dimension, dimension):
                raise HermiteBorderBasisError("a chart Hermite matrix shape changed")
            if stored["multiplication_matrices"].shape != (
                VARIABLE_COUNT,
                dimension,
                dimension,
            ):
                raise HermiteBorderBasisError("a chart multiplication shape changed")
            if stored["point_idempotent_coefficients"].shape != (
                dimension,
                support_count,
            ):
                raise HermiteBorderBasisError("a chart idempotent shape changed")

    for row in overlaps:
        artifact = path.parent / row["artifact"]
        if sha256(artifact) != row["artifact_sha256"]:
            raise HermiteBorderBasisError("an overlap artifact hash changed")
        with np.load(artifact) as stored:
            dimension = int(row["overlap_quotient_length"])
            if stored["standard_monomials"].shape != (dimension, VARIABLE_COUNT):
                raise HermiteBorderBasisError("an overlap standard-monomial shape changed")
            if stored["source_restriction"].shape != (
                dimension,
                int(row["source_quotient_length"]),
            ):
                raise HermiteBorderBasisError("a source restriction shape changed")
            if stored["target_restriction"].shape != (
                dimension,
                int(row["target_quotient_length"]),
            ):
                raise HermiteBorderBasisError("a target restriction shape changed")
    return data


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface", type=int, default=0)
    parser.add_argument("--chart", default="C0000")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--without-overlaps", action="store_true")
    parser.add_argument("--output-directory", type=Path, default=OUTPUT_DIRECTORY)
    args = parser.parse_args(argv)
    if args.verify:
        data = verify_manifest(args.output_directory / "hermite_border_basis_manifest.json")
        print(
            json.dumps(
                {"charts": len(data["charts"]), "overlaps": len(data["overlaps"])},
                sort_keys=True,
            )
        )
        return 0
    context = iw_solver.load_problem_context()
    prepare_border_context(context)
    if args.all:
        manifest = run_all(
            context,
            args.output_directory,
            include_overlaps=not args.without_overlaps,
        )
        print(manifest)
        return 0
    result = chart_border_basis(context, args.surface, args.chart)
    path, summary = write_chart_artifact(result, args.output_directory)
    print(json.dumps(result["summary"], indent=2, sort_keys=True))
    print(path)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
