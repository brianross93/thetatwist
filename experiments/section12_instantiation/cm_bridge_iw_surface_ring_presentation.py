"""Build one exact surface-ring presentation of the Hermite ideal ``I_W``.

The persisted Hermite artifact gives a finite algebra ``A_W`` together with
eight commuting coordinate multiplication matrices.  Over the affine
polynomial ring ``P`` the graph equations ``x_i I-M_i`` have an exact Koszul
resolution.  After base change to the six-equation surface ring, the first
homology is ``A_W^6``.  Six divided-difference cycles, one for each surface
equation, kill exactly those six classes.

Contracting an order-ideal spanning tree and identifying repeated border
transitions first gives ``S^2750 -> S^297 -> I_W``.  Exact unit Tietze
cancellations then remove the 127 non-corner border relations and give

    S^1887 -> S^170 -> I_W -> 0

on the representative chart ``Z_2/C0010``.  Every entry is retained as an
exact sparse polynomial over ``F_31``.  The producer checks that the first
2744 composites are zero already over ``P`` and that the last six composites
are the six defining equations, hence zero over ``S``.

This is a presentation of the independently certified Hermite ideal.  It
does not use the failed two- or three-generator translated-section ideals.
The left map is not claimed injective; this is not a length-one resolution.
"""

from __future__ import annotations

import argparse
import hashlib
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
SURFACE = 2
CHART = "C0010"
ZERO_EXPONENT = (0,) * VARIABLE_COUNT
BORDER_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_iw_hermite_border_basis"
)
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_iw_surface_ring_presentation"
)


class IWSurfaceRingPresentationError(RuntimeError):
    """An upstream invariant or exact polynomial identity failed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _add(poly, exponent, coefficient):
    exponent = tuple(map(int, exponent))
    value = (int(poly.get(exponent, 0)) + int(coefficient)) % P
    if value:
        poly[exponent] = value
    else:
        poly.pop(exponent, None)


def _scaled(poly, scalar):
    scalar = int(scalar) % P
    if not scalar:
        return {}
    return {
        exponent: int(coefficient) * scalar % P
        for exponent, coefficient in poly.items()
        if int(coefficient) * scalar % P
    }


def _multiply(left, right):
    result = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                int(a) + int(b) for a, b in zip(left_exponent, right_exponent)
            )
            _add(result, exponent, int(left_coefficient) * int(right_coefficient))
    return result


def _add_scaled(target, source, scalar=1):
    for exponent, coefficient in source.items():
        _add(target, exponent, int(scalar) * int(coefficient))


def _unit_exponent(variable):
    exponent = [0] * VARIABLE_COUNT
    exponent[int(variable)] = 1
    return tuple(exponent)


def _matrix_power(matrix, power):
    result = np.eye(matrix.shape[0], dtype=np.int64)
    base = np.asarray(matrix, dtype=np.int64) % P
    for _ in range(int(power)):
        result = result @ base % P
    return result


def _evaluate_matrix_polynomial(polynomial, multiplication):
    dimension = multiplication.shape[1]
    result = np.zeros((dimension, dimension), dtype=np.int64)
    powers = {
        (variable, power): _matrix_power(multiplication[variable], power)
        for variable in range(VARIABLE_COUNT)
        for power in range(
            1
            + max(
                [int(exponent[variable]) for exponent in polynomial]
                or [0]
            )
        )
    }
    for exponent, coefficient in polynomial.items():
        operator = np.eye(dimension, dtype=np.int64)
        for variable, power in enumerate(exponent):
            operator = operator @ powers[(variable, int(power))] % P
        result = (result + int(coefficient) * operator) % P
    return result


def load_input():
    path = BORDER_DIRECTORY / f"Z_{SURFACE}_{CHART}.npz"
    summary_path = BORDER_DIRECTORY / f"Z_{SURFACE}_{CHART}.json"
    if not path.exists() or not summary_path.exists():
        raise IWSurfaceRingPresentationError("the bound Hermite chart is absent")
    with np.load(path) as artifact:
        field = int(artifact["field"][0])
        standard = np.asarray(artifact["standard_monomials"], dtype=np.int64)
        multiplication = np.asarray(
            artifact["multiplication_matrices"], dtype=np.int64
        ) % P
        border = np.asarray(artifact["border_monomials"], dtype=np.int64)
        normal_forms = np.asarray(
            artifact["border_normal_form_coefficients"], dtype=np.int64
        ) % P
    if field != P:
        raise IWSurfaceRingPresentationError("the Hermite artifact field changed")
    if standard.shape != (98, 8) or multiplication.shape != (8, 98, 98):
        raise IWSurfaceRingPresentationError("the Hermite graph dimensions changed")
    if border.shape != (297, 8) or normal_forms.shape != (297, 98):
        raise IWSurfaceRingPresentationError("the Hermite border dimensions changed")

    context = iw_solver.load_problem_context()
    bits = context["chart_bits"][CHART]
    labels = list(
        map(int, json.loads(context["seeds"][SURFACE]["surface_section_labels"]))
    )
    surface_equations = [
        *iw_solver.curve_polynomials(bits),
        *[
            iw_solver.theta_polynomial(context["theta"][label], bits)
            for label in labels
        ],
    ]
    if len(surface_equations) != 6:
        raise IWSurfaceRingPresentationError("the surface stopped having six equations")

    regular_path = ROOT / (
        "results/section12_instantiation/cm_bridge_o8_section_ring_extension/"
        "o8_section_ring_extension_certificate.json"
    )
    regular = json.loads(regular_path.read_text(encoding="utf-8"))
    if regular.get("schema") != "section12-o8-section-ring-extension-v1":
        raise IWSurfaceRingPresentationError(
            "the bound surface regular-sequence certificate changed schema"
        )
    if regular.get("O8_surface_Koszul_contract", {}).get(
        "exactness_source"
    ) != (
        "the two certified theta equations form the regular sequence defining "
        "the smooth complete-intersection surface"
    ):
        raise IWSurfaceRingPresentationError(
            "the bound surface regular-sequence hypothesis changed"
        )
    return {
        "path": path,
        "summary_path": summary_path,
        "standard": [tuple(map(int, row)) for row in standard],
        "multiplication": multiplication,
        "border": [tuple(map(int, row)) for row in border],
        "normal_forms": normal_forms,
        "surface_equations": surface_equations,
        "surface_labels": labels,
        "regular_path": regular_path,
    }


def border_generators(data):
    generators = []
    for border_exponent, coefficients in zip(
        data["border"], data["normal_forms"]
    ):
        polynomial = {border_exponent: 1}
        for standard_exponent, coefficient in zip(data["standard"], coefficients):
            if int(coefficient) % P:
                _add(polynomial, standard_exponent, -int(coefficient))
        generators.append(polynomial)
    return generators


def transition_quotient(data):
    """Map 784 graph relations to 297 unique border relations or zero."""

    standard_index = {exponent: index for index, exponent in enumerate(data["standard"])}
    border_index = {exponent: index for index, exponent in enumerate(data["border"])}
    quotient = -np.ones((VARIABLE_COUNT, len(data["standard"])), dtype=np.int64)
    internal = 0
    border_transitions = 0
    for variable in range(VARIABLE_COUNT):
        for source, exponent in enumerate(data["standard"]):
            product = list(exponent)
            product[variable] += 1
            product = tuple(product)
            if product in standard_index:
                internal += 1
                continue
            if product not in border_index:
                raise IWSurfaceRingPresentationError(
                    "a graph transition is neither standard nor in the border"
                )
            quotient[variable, source] = border_index[product]
            border_transitions += 1
    if internal != 195 or border_transitions != 589:
        raise IWSurfaceRingPresentationError("the graph transition census changed")
    if set(map(int, quotient[quotient >= 0])) != set(range(len(data["border"]))):
        raise IWSurfaceRingPresentationError("the graph-to-border map is not surjective")

    # One divisor edge for every nonunit standard monomial is a triangular
    # spanning tree.  It proves that the graph generator contraction is a
    # genuine Tietze contraction, not a dimension count.
    standard_set = set(data["standard"])
    tree = []
    for target, exponent in enumerate(data["standard"]):
        if exponent == ZERO_EXPONENT:
            continue
        for variable, power in enumerate(exponent):
            if not power:
                continue
            parent = list(exponent)
            parent[variable] -= 1
            parent = tuple(parent)
            if parent in standard_set:
                tree.append((variable, data["standard"].index(parent), target))
                break
        else:
            raise IWSurfaceRingPresentationError("the order ideal lost a parent")
    if len(tree) != len(data["standard"]) - 1:
        raise IWSurfaceRingPresentationError("the order-ideal tree changed size")
    return quotient, tree, internal, border_transitions


def _add_matrix_term(column, row, exponent, coefficient):
    row = int(row)
    if row < 0:
        return
    polynomial = column.setdefault(row, {})
    _add(polynomial, exponent, coefficient)
    if not polynomial:
        column.pop(row, None)


def koszul_columns(data, quotient):
    """Project the 2744 graph-Koszul columns to border generators."""

    multiplication = data["multiplication"]
    dimension = len(data["standard"])
    columns = []

    def add_L(column, block_variable, operator_variable, source, scale):
        target_row = int(quotient[block_variable, source])
        if target_row >= 0:
            _add_matrix_term(
                column,
                target_row,
                _unit_exponent(operator_variable),
                scale,
            )
        matrix_column = multiplication[operator_variable, :, source]
        for target in np.flatnonzero(matrix_column):
            border_row = int(quotient[block_variable, int(target)])
            if border_row >= 0:
                _add_matrix_term(
                    column,
                    border_row,
                    ZERO_EXPONENT,
                    -int(scale) * int(matrix_column[target]),
                )

    for left in range(VARIABLE_COUNT):
        for right in range(left + 1, VARIABLE_COUNT):
            for source in range(dimension):
                column = {}
                # d(e_left wedge e_right) = L_right e_left - L_left e_right.
                add_L(column, left, right, source, 1)
                add_L(column, right, left, source, -1)
                columns.append(column)
    if len(columns) != dimension * 28:
        raise IWSurfaceRingPresentationError("the Koszul column count changed")
    return columns


def divided_difference_cycle(polynomial, multiplication, unit_index):
    """Return h with sum_i (x_i I-M_i) h_i = f(x)e_unit."""

    dimension = multiplication.shape[1]
    maximum_powers = [
        max([int(exponent[variable]) for exponent in polynomial] or [0])
        for variable in range(VARIABLE_COUNT)
    ]
    powers = {
        (variable, power): _matrix_power(multiplication[variable], power)
        for variable in range(VARIABLE_COUNT)
        for power in range(maximum_powers[variable] + 1)
    }
    unit = np.zeros(dimension, dtype=np.int64)
    unit[unit_index] = 1
    blocks = [[{} for _ in range(dimension)] for _ in range(VARIABLE_COUNT)]
    identity = np.eye(dimension, dtype=np.int64)

    for exponent, coefficient in polynomial.items():
        exponent = tuple(map(int, exponent))
        for variable, power in enumerate(exponent):
            if not power:
                continue
            left = identity
            for earlier in range(variable):
                left = left @ powers[(earlier, exponent[earlier])] % P
            for t in range(power):
                vector = left @ powers[(variable, t)] @ unit % P
                polynomial_exponent = [0] * VARIABLE_COUNT
                polynomial_exponent[variable] = power - 1 - t
                for later in range(variable + 1, VARIABLE_COUNT):
                    polynomial_exponent[later] = exponent[later]
                for target in np.flatnonzero(vector):
                    _add(
                        blocks[variable][int(target)],
                        tuple(polynomial_exponent),
                        int(coefficient) * int(vector[target]),
                    )
    return blocks


def graph_d1(cycle, multiplication):
    dimension = multiplication.shape[1]
    result = [{} for _ in range(dimension)]
    for variable, block in enumerate(cycle):
        x_exponent = _unit_exponent(variable)
        for source, polynomial in enumerate(block):
            if not polynomial:
                continue
            for exponent, coefficient in polynomial.items():
                shifted = tuple(a + b for a, b in zip(exponent, x_exponent))
                _add(result[source], shifted, coefficient)
            matrix_column = multiplication[variable, :, source]
            for target in np.flatnonzero(matrix_column):
                _add_scaled(
                    result[int(target)],
                    polynomial,
                    -int(matrix_column[target]),
                )
    return result


def surface_columns(data, quotient):
    unit_index = data["standard"].index(ZERO_EXPONENT)
    columns = []
    graph_cycles = []
    for equation in data["surface_equations"]:
        if np.any(_evaluate_matrix_polynomial(equation, data["multiplication"])):
            raise IWSurfaceRingPresentationError(
                "a surface equation does not annihilate the Hermite algebra"
            )
        cycle = divided_difference_cycle(
            equation, data["multiplication"], unit_index
        )
        image = graph_d1(cycle, data["multiplication"])
        expected = [dict() for _ in data["standard"]]
        expected[unit_index] = {
            tuple(map(int, exponent)): int(coefficient) % P
            for exponent, coefficient in equation.items()
            if int(coefficient) % P
        }
        if image != expected:
            raise IWSurfaceRingPresentationError(
                "a divided-difference cycle has the wrong graph boundary"
            )
        projected = {}
        for variable, block in enumerate(cycle):
            for source, polynomial in enumerate(block):
                border_row = int(quotient[variable, source])
                if border_row < 0:
                    continue
                for exponent, coefficient in polynomial.items():
                    _add_matrix_term(
                        projected, border_row, exponent, coefficient
                    )
        columns.append(projected)
        graph_cycles.append(cycle)
    return columns, graph_cycles


def compose(generator_polynomials, column):
    result = {}
    for row, coefficient_polynomial in column.items():
        _add_scaled(
            result,
            _multiply(generator_polynomials[int(row)], coefficient_polynomial),
        )
    return result


def verify_composition(generators, columns, surface_equations):
    koszul_count = 28 * 98
    maximum_terms = 0
    for index, column in enumerate(columns):
        composite = compose(generators, column)
        maximum_terms = max(maximum_terms, len(composite))
        if index < koszul_count:
            expected = {}
        else:
            expected = {
                tuple(map(int, exponent)): int(coefficient) % P
                for exponent, coefficient in surface_equations[index - koszul_count].items()
                if int(coefficient) % P
            }
        if composite != expected:
            raise IWSurfaceRingPresentationError(
                f"presentation composite failed at column {index}"
            )
    return maximum_terms


def pure_constant_pivots(columns):
    pivots = set()
    occurrences = 0
    for column_index, column in enumerate(columns):
        for row, polynomial in column.items():
            if set(polynomial) == {ZERO_EXPONENT}:
                pivots.add(int(row))
                occurrences += 1
    return len(pivots), occurrences


def corner_rows(border):
    border_set = set(border)
    result = []
    for row, exponent in enumerate(border):
        is_corner = True
        for variable, power in enumerate(exponent):
            if not power:
                continue
            divisor = list(exponent)
            divisor[variable] -= 1
            if tuple(divisor) in border_set:
                is_corner = False
                break
        if is_corner:
            result.append(row)
    return result


def _column_term_count(column):
    return sum(len(polynomial) for polynomial in column.values())


def cancel_noncorner_units(generators, border, columns, surface_equations):
    """Tietze-cancel every non-corner border relation through a unit pivot.

    The corner set is intrinsic to the selected order ideal.  A cancellation
    is performed only when the relevant matrix entry is a *pure* nonzero
    field constant.  Hence every step is an invertible elementary operation
    over the surface ring, not a rank heuristic.
    """

    corners = set(corner_rows(border))
    noncorners = sorted(
        set(range(len(border))) - corners,
        key=lambda row: border[row],
        reverse=True,
    )
    active_rows = set(range(len(border)))
    active_columns = set(range(len(columns)))
    row_columns = [set() for _ in border]
    for column_index, column in enumerate(columns):
        for row in column:
            row_columns[int(row)].add(column_index)

    # Boundary ledgers express each P-level composite as a combination of the
    # six surface equations.  Koszul columns start at zero; the last six start
    # at the corresponding standard basis vector.
    boundaries = [{} for _ in columns]
    for equation_index in range(len(surface_equations)):
        boundaries[2744 + equation_index] = {
            equation_index: {ZERO_EXPONENT: 1}
        }

    cancellation_rows = []
    cancellation_columns = []
    for row in noncorners:
        candidates = [
            column_index
            for column_index in row_columns[row]
            if column_index in active_columns
            and set(columns[column_index].get(row, {})) == {ZERO_EXPONENT}
        ]
        if not candidates:
            raise IWSurfaceRingPresentationError(
                f"non-corner border row {row} has no exact unit pivot"
            )
        # Prefer a zero-boundary Koszul column and then the smallest exact
        # update.  This keeps the surface-equation ledger sparse.
        pivot_column = min(
            candidates,
            key=lambda column_index: (
                bool(boundaries[column_index]),
                _column_term_count(columns[column_index]) * len(row_columns[row]),
                _column_term_count(columns[column_index]),
                len(columns[column_index]),
                column_index,
            ),
        )
        pivot_coefficient = int(
            columns[pivot_column][row][ZERO_EXPONENT]
        ) % P
        inverse = pow(pivot_coefficient, -1, P)
        for target_row in list(columns[pivot_column]):
            columns[pivot_column][target_row] = _scaled(
                columns[pivot_column][target_row], inverse
            )
        for equation_index in list(boundaries[pivot_column]):
            boundaries[pivot_column][equation_index] = _scaled(
                boundaries[pivot_column][equation_index], inverse
            )
        pivot = {
            target_row: dict(polynomial)
            for target_row, polynomial in columns[pivot_column].items()
        }
        pivot_boundary = {
            equation_index: dict(polynomial)
            for equation_index, polynomial in boundaries[pivot_column].items()
        }

        for other_column in list(row_columns[row] - {pivot_column}):
            if other_column not in active_columns:
                continue
            multiplier = columns[other_column].pop(row)
            row_columns[row].discard(other_column)
            for target_row, pivot_polynomial in pivot.items():
                if target_row == row:
                    continue
                existed = target_row in columns[other_column]
                product = _multiply(multiplier, pivot_polynomial)
                if not existed:
                    columns[other_column][target_row] = {}
                _add_scaled(columns[other_column][target_row], product, -1)
                if not columns[other_column][target_row]:
                    columns[other_column].pop(target_row, None)
                    if existed:
                        row_columns[target_row].discard(other_column)
                elif not existed:
                    row_columns[target_row].add(other_column)
            for equation_index, pivot_polynomial in pivot_boundary.items():
                product = _multiply(multiplier, pivot_polynomial)
                boundary_polynomial = boundaries[other_column].setdefault(
                    equation_index, {}
                )
                _add_scaled(boundary_polynomial, product, -1)
                if not boundary_polynomial:
                    boundaries[other_column].pop(equation_index, None)

        for target_row in list(columns[pivot_column]):
            row_columns[target_row].discard(pivot_column)
        active_rows.remove(row)
        active_columns.remove(pivot_column)
        row_columns[row].clear()
        columns[pivot_column] = {}
        boundaries[pivot_column] = {}
        cancellation_rows.append(row)
        cancellation_columns.append(pivot_column)

    if active_rows != corners:
        raise IWSurfaceRingPresentationError("the Tietze reduction did not leave the corners")
    if any(column >= 2744 for column in cancellation_columns):
        raise IWSurfaceRingPresentationError("a surface Tor column was unexpectedly cancelled")

    # Verify the transformed matrix against its exact surface-equation ledger.
    for column_index in sorted(active_columns):
        composite = compose(generators, columns[column_index])
        expected = {}
        for equation_index, coefficient_polynomial in boundaries[column_index].items():
            _add_scaled(
                expected,
                _multiply(coefficient_polynomial, surface_equations[equation_index]),
            )
        if composite != expected:
            raise IWSurfaceRingPresentationError(
                f"the reduced presentation composite failed at column {column_index}"
            )

    row_order = sorted(active_rows)
    column_order = sorted(active_columns)
    row_relabel = {old: new for new, old in enumerate(row_order)}
    reduced_columns = [
        {
            row_relabel[int(row)]: polynomial
            for row, polynomial in columns[column_index].items()
        }
        for column_index in column_order
    ]
    reduced_generators = [generators[row] for row in row_order]
    reduced_boundaries = [boundaries[column_index] for column_index in column_order]
    return {
        "generators": reduced_generators,
        "columns": reduced_columns,
        "boundaries": reduced_boundaries,
        "original_corner_rows": row_order,
        "original_remaining_columns": column_order,
        "cancelled_rows": cancellation_rows,
        "cancelled_columns": cancellation_columns,
    }


def prune_zero_and_field_proportional_columns(columns, boundaries, original_columns):
    """Remove zero columns and repeated columns up to an F_31 scalar.

    This changes neither the generated syzygy module nor exactness.  It is a
    domain-only compression; no claim of a minimal first-syzygy module is
    made.
    """

    signatures = {}
    kept_columns = []
    kept_boundaries = []
    kept_original = []
    zero_count = 0
    proportional_count = 0
    zero_original_columns = []
    proportional_ledger = []
    for column, boundary, original in zip(columns, boundaries, original_columns):
        terms = [
            (int(row), tuple(map(int, exponent)), int(coefficient) % P)
            for row, polynomial in sorted(column.items())
            for exponent, coefficient in sorted(polynomial.items())
        ]
        if not terms:
            zero_count += 1
            zero_original_columns.append(int(original))
            continue
        inverse = pow(terms[0][2], -1, P)
        signature = tuple(
            (row, exponent, coefficient * inverse % P)
            for row, exponent, coefficient in terms
        )
        if signature in signatures:
            proportional_count += 1
            representative_index, representative_first = signatures[signature]
            scalar = terms[0][2] * pow(int(representative_first), -1, P) % P
            proportional_ledger.append(
                (
                    int(original),
                    int(kept_original[representative_index]),
                    int(scalar),
                )
            )
            continue
        signatures[signature] = (len(kept_columns), terms[0][2])
        kept_columns.append(column)
        kept_boundaries.append(boundary)
        kept_original.append(original)
    return {
        "columns": kept_columns,
        "boundaries": kept_boundaries,
        "original_columns": kept_original,
        "zero_columns_removed": zero_count,
        "field_proportional_columns_removed": proportional_count,
        "zero_original_columns": zero_original_columns,
        "proportional_ledger": proportional_ledger,
    }


def _flatten_generators(generators):
    rows = []
    coefficients = []
    exponents = []
    for row, polynomial in enumerate(generators):
        for exponent, coefficient in sorted(polynomial.items()):
            rows.append(row)
            coefficients.append(int(coefficient) % P)
            exponents.append(exponent)
    return (
        np.asarray(rows, dtype=np.uint16),
        np.asarray(coefficients, dtype=np.uint8),
        np.asarray(exponents, dtype=np.uint8),
    )


def _flatten_matrix(columns):
    rows = []
    column_indices = []
    coefficients = []
    exponents = []
    for column_index, column in enumerate(columns):
        for row, polynomial in sorted(column.items()):
            for exponent, coefficient in sorted(polynomial.items()):
                rows.append(int(row))
                column_indices.append(column_index)
                coefficients.append(int(coefficient) % P)
                exponents.append(exponent)
    return (
        np.asarray(rows, dtype=np.uint16),
        np.asarray(column_indices, dtype=np.uint16),
        np.asarray(coefficients, dtype=np.uint8),
        np.asarray(exponents, dtype=np.uint8),
    )


def build_presentation():
    data = load_input()
    generators = border_generators(data)
    quotient, tree, internal, border_transitions = transition_quotient(data)
    koszul = koszul_columns(data, quotient)
    surface, _graph_cycles = surface_columns(data, quotient)
    columns = [*koszul, *surface]
    maximum_composite_terms = verify_composition(
        generators, columns, data["surface_equations"]
    )
    pivot_rows, pivot_occurrences = pure_constant_pivots(columns)
    reduced = cancel_noncorner_units(
        generators,
        data["border"],
        columns,
        data["surface_equations"],
    )
    pruned = prune_zero_and_field_proportional_columns(
        reduced["columns"],
        reduced["boundaries"],
        reduced["original_remaining_columns"],
    )
    return {
        "data": data,
        "generators": reduced["generators"],
        "quotient": quotient,
        "tree": tree,
        "columns": pruned["columns"],
        "boundaries": pruned["boundaries"],
        "original_corner_rows": reduced["original_corner_rows"],
        "original_remaining_columns": pruned["original_columns"],
        "cancelled_rows": reduced["cancelled_rows"],
        "cancelled_columns": reduced["cancelled_columns"],
        "removed_zero_original_columns": pruned["zero_original_columns"],
        "removed_proportional_original_columns": pruned["proportional_ledger"],
        "counts": {
            "standard_monomials": 98,
            "original_border_generators": 297,
            "corner_border_generators": len(reduced["generators"]),
            "graph_C0_rank": 98,
            "graph_C1_rank": 784,
            "graph_Koszul_C2_rank": 2744,
            "surface_Tor_cycles": 6,
            "pre_Tietze_relation_columns": 2750,
            "Tietze_unit_cancellations": len(reduced["cancelled_rows"]),
            "post_Tietze_relation_columns": len(reduced["columns"]),
            "zero_relation_columns_removed": pruned["zero_columns_removed"],
            "field_proportional_relation_columns_removed": pruned[
                "field_proportional_columns_removed"
            ],
            "presentation_relation_columns": len(pruned["columns"]),
            "order_ideal_tree_relations_contracted": len(tree),
            "graph_internal_transitions": internal,
            "graph_border_transitions": border_transitions,
            "repeated_border_transitions_merged": border_transitions - 297,
            "non_tree_internal_zero_transitions_removed": internal - len(tree),
            "pure_constant_pivot_rows_available_for_further_Tietze_reduction": pivot_rows,
            "pure_constant_pivot_occurrences": pivot_occurrences,
            "maximum_terms_in_verified_composite": maximum_composite_terms,
        },
    }


def write_artifacts(result, output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    data = result["data"]
    generator_row, generator_coefficient, generator_exponent = _flatten_generators(
        result["generators"]
    )
    (
        syzygy_row,
        syzygy_column,
        syzygy_coefficient,
        syzygy_exponent,
    ) = _flatten_matrix(result["columns"])
    surface_term_equation = []
    surface_term_coefficient = []
    surface_term_exponent = []
    for equation_index, equation in enumerate(data["surface_equations"]):
        for exponent, coefficient in sorted(equation.items()):
            surface_term_equation.append(equation_index)
            surface_term_coefficient.append(int(coefficient) % P)
            surface_term_exponent.append(exponent)

    boundary_term_column = []
    boundary_term_equation = []
    boundary_term_coefficient = []
    boundary_term_exponent = []
    for column_index, boundary in enumerate(result["boundaries"]):
        for equation_index, polynomial in sorted(boundary.items()):
            for exponent, coefficient in sorted(polynomial.items()):
                boundary_term_column.append(column_index)
                boundary_term_equation.append(equation_index)
                boundary_term_coefficient.append(int(coefficient) % P)
                boundary_term_exponent.append(exponent)

    artifact_path = output_directory / "Z_2_C0010_surface_ring_presentation.npz"
    np.savez_compressed(
        artifact_path,
        field=np.asarray([P], dtype=np.uint8),
        matrix_shape=np.asarray([170, 1887], dtype=np.uint16),
        standard_monomials=np.asarray(data["standard"], dtype=np.uint8),
        border_monomials=np.asarray(data["border"], dtype=np.uint8),
        graph_to_border=np.asarray(result["quotient"], dtype=np.int16),
        order_ideal_tree=np.asarray(result["tree"], dtype=np.uint16),
        original_corner_rows=np.asarray(result["original_corner_rows"], dtype=np.uint16),
        original_remaining_columns=np.asarray(
            result["original_remaining_columns"], dtype=np.uint16
        ),
        cancelled_noncorner_rows=np.asarray(result["cancelled_rows"], dtype=np.uint16),
        cancelled_unit_pivot_columns=np.asarray(
            result["cancelled_columns"], dtype=np.uint16
        ),
        removed_zero_original_columns=np.asarray(
            result["removed_zero_original_columns"], dtype=np.uint16
        ),
        removed_proportional_original_columns=np.asarray(
            result["removed_proportional_original_columns"], dtype=np.uint16
        ),
        generator_term_row=generator_row,
        generator_term_coefficient=generator_coefficient,
        generator_term_exponent=generator_exponent,
        syzygy_term_row=syzygy_row,
        syzygy_term_column=syzygy_column,
        syzygy_term_coefficient=syzygy_coefficient,
        syzygy_term_exponent=syzygy_exponent,
        surface_term_equation=np.asarray(surface_term_equation, dtype=np.uint8),
        surface_term_coefficient=np.asarray(surface_term_coefficient, dtype=np.uint8),
        surface_term_exponent=np.asarray(surface_term_exponent, dtype=np.uint8),
        boundary_term_column=np.asarray(boundary_term_column, dtype=np.uint16),
        boundary_term_equation=np.asarray(boundary_term_equation, dtype=np.uint8),
        boundary_term_coefficient=np.asarray(
            boundary_term_coefficient, dtype=np.uint8
        ),
        boundary_term_exponent=np.asarray(boundary_term_exponent, dtype=np.uint8),
    )

    regular_path = data["regular_path"]
    certificate = {
        "schema": "section12-IW-surface-ring-first-presentation-v1",
        "field": "F_31",
        "surface": "Z_2",
        "chart": CHART,
        "surface_ring": (
            "S=F_31[u0,v0,...,u3,v3]/(four affine elliptic-curve equations,"
            " two theta equations defining Z_2)"
        ),
        "exact_presentation": "S^1887 -> S^170 -> I_W -> 0",
        "counts": result["counts"],
        "construction": {
            "polynomial_graph_complex": "Koszul complex of x_i I-M_i on P tensor_F31 A_W",
            "graph_exactness": (
                "x_i tensor 1-1 tensor M_i are successively monic, hence a regular sequence"
            ),
            "base_change_H1": (
                "Tor_1^P(S,A_W)=A_W^6 because the six defining equations form a P-regular sequence and annihilate A_W"
            ),
            "six_extra_columns": (
                "explicit divided-difference homotopies whose graph boundaries are the six defining equations times the unit"
            ),
            "contraction": (
                "97 triangular order-ideal tree relations are contracted; 589 border transitions are merged onto 297 distinct Hermite border generators; 127 exact pure-unit Tietze pivots leave the 170 border corners; zero and F_31-proportional domain columns are pruned"
            ),
        },
        "exact_checks": {
            "all_eight_coordinate_multiplications_bound_from_Hermite_artifact": True,
            "all_six_surface_equations_annihilate_A_W": True,
            "graph_to_border_map_surjective": True,
            "order_ideal_tree_has_97_edges": True,
            "2744_Koszul_composites_zero_over_polynomial_ring": True,
            "six_divided_difference_composites_equal_surface_equations": True,
            "127_noncorner_relations_removed_by_explicit_unit_Tietze_pivots": True,
            "reduced_2623_columns_compose_to_the_stored_surface_equation_ledger": True,
            "zero_and_field_proportional_domain_columns_pruned_without_changing_image": True,
            "therefore_d_squared_zero_over_surface_ring": True,
            "first_presentation_exact_by_regular_sequence_change_of_rings_and_Tietze_contraction": True,
        },
        "upstream": {
            "Hermite_binary": {
                "path": str(data["path"].relative_to(ROOT)),
                "sha256": sha256(data["path"]),
            },
            "Hermite_summary": {
                "path": str(data["summary_path"].relative_to(ROOT)),
                "sha256": sha256(data["summary_path"]),
            },
            "surface_regular_sequence_certificate": {
                "path": str(regular_path.relative_to(ROOT)),
                "sha256": sha256(regular_path),
            },
        },
        "binary_artifact": {
            "path": artifact_path.name,
            "sha256": sha256(artifact_path),
            "encoding": (
                "COO polynomial terms: each matrix entry is a sparse F_31 polynomial in eight affine variables"
            ),
            "syzygy_column_order": (
                "the retained original column indices are stored explicitly; before 127 cancellations they were 28 variable pairs times 98 standard monomials followed by six surface divided-difference columns"
            ),
        },
        "claim_boundary": {
            "proved": (
                "one exact finite 170-generator surface-ring first presentation of the actual Hermite ideal I_W on Z_2/C0010"
            ),
            "not_claimed": (
                "injectivity of the S^1887 left map, minimality of the 170 corner generators, a Hilbert-Burch-minimal matrix, an open-overlap localization, or the final Hartshorne-Serre cross-Hom endpoint"
            ),
            "next_exact_reduction": (
                "construct the syzygy-of-syzygy/Ext cocycle from this first presentation, or contract it to a projective/Hilbert-Burch model before forming the Hartshorne-Serre extension"
            ),
        },
    }
    certificate_path = output_directory / "surface_ring_presentation_certificate.json"
    certificate_path.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return certificate_path, artifact_path


def reproduce(output_directory=OUTPUT_DIRECTORY):
    result = build_presentation()
    paths = write_artifacts(result, output_directory)
    return result, paths


def verify_artifacts(output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    certificate_path = output_directory / "surface_ring_presentation_certificate.json"
    artifact_path = output_directory / "Z_2_C0010_surface_ring_presentation.npz"
    certificate = json.loads(certificate_path.read_text(encoding="utf-8"))
    if certificate["schema"] != "section12-IW-surface-ring-first-presentation-v1":
        raise IWSurfaceRingPresentationError("the certificate schema changed")
    if certificate["binary_artifact"]["sha256"] != sha256(artifact_path):
        raise IWSurfaceRingPresentationError("the presentation artifact hash changed")
    with np.load(artifact_path) as artifact:
        if list(map(int, artifact["matrix_shape"])) != [170, 1887]:
            raise IWSurfaceRingPresentationError("the persisted matrix shape changed")
        if artifact["syzygy_term_exponent"].shape[1] != VARIABLE_COUNT:
            raise IWSurfaceRingPresentationError("the persisted exponent width changed")
        if int(artifact["syzygy_term_column"].max()) != 1886:
            raise IWSurfaceRingPresentationError("the last reduced column is absent")
    return certificate


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--output-directory", type=Path, default=OUTPUT_DIRECTORY)
    args = parser.parse_args(argv)
    if args.verify:
        certificate = verify_artifacts(args.output_directory)
        print(json.dumps(certificate["counts"], sort_keys=True))
        return 0
    result, paths = reproduce(args.output_directory)
    print(f"wrote {paths[0]}")
    print(f"wrote {paths[1]}")
    print(json.dumps(result["counts"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
