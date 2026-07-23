"""Transport the Hartshorne--Serre Ext^2 row across one cube edge.

The chart-local cochain is

    Lambda_C = [lambda_C, 0^6] : P2_C -> O(4H)_C,

where ``lambda_C=-tau_C H_5 ... H_0`` in the repository orientation.
For an oriented cube edge ``source -> target`` this producer forms

    Delta = Lambda_edge J2 - h^4 phi(Lambda_target).

It then constructs one degree-one cochain ``a_target`` and requires the
single identity ``Delta=a_target D2_target``.  In particular, the six
surface-equation columns are checks on the same cochain; they are never
repaired independently.

This is an edge pilot.  It does not assert face/cube coherence or global
Hartshorne--Serre descent.
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
    cm_bridge_hs_border_to_local_frames as hs_frames,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_extension_from_surface_presentation as hs_extension,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_shamash_overlap_lifts as shamash,
)


P = 31
VARIABLE_COUNT = 8
ZERO_EXPONENT = (0,) * VARIABLE_COUNT
DEFAULT_SURFACE = 0
DEFAULT_EDGE = 2
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_lambda_overlap_transport"
)


class HSLambdaTransportError(RuntimeError):
    """An orientation, cochain, or exact comparison identity failed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def chart_parity(chart: str) -> int:
    if not chart.startswith("C") or len(chart) != 5:
        raise HSLambdaTransportError("a chart label changed")
    return -1 if sum(map(int, chart[1:])) % 2 else 1


def _add_matrix_term(polynomial, exponent, matrix, scale=1):
    exponent = tuple(map(int, exponent))
    value = int(scale) * np.asarray(matrix, dtype=np.int64) % P
    if exponent in polynomial:
        value = (polynomial[exponent] + value) % P
    if np.any(value):
        polynomial[exponent] = value
    else:
        polynomial.pop(exponent, None)


def _add_row_map(target, source, scale=1):
    for wedge, polynomial in source.items():
        value = target.setdefault(tuple(wedge), {})
        for exponent, matrix in polynomial.items():
            _add_matrix_term(value, exponent, matrix, scale)
        if not value:
            target.pop(tuple(wedge), None)


def _compose_row_with_blocks(row, blocks):
    """Compose a 1-by-2 block row with a 2-by-2 block map."""

    result = {}
    for (output_wedge, input_wedge), block_polynomial in blocks.items():
        left = row.get(tuple(output_wedge), {})
        if not left:
            continue
        target = result.setdefault(tuple(input_wedge), {})
        for left_exponent, left_matrix in left.items():
            for right_exponent, right_matrix in block_polynomial.items():
                exponent = tuple(
                    int(a) + int(b)
                    for a, b in zip(left_exponent, right_exponent)
                )
                _add_matrix_term(
                    target, exponent, np.asarray(left_matrix) @ right_matrix
                )
        if not target:
            result.pop(tuple(input_wedge), None)
    return result


def _transform_row_exponents(row, changed):
    """Apply the involution between target and source Laurent coordinates."""

    result = {}
    for wedge, polynomial in row.items():
        target = result.setdefault(tuple(wedge), {})
        for exponent, matrix in polynomial.items():
            _add_matrix_term(
                target, shamash.pull_target_exponent(exponent, changed), matrix
            )
        if not target:
            result.pop(tuple(wedge), None)
    return result


def _shift_row(row, exponent, scale=1):
    result = {}
    for wedge, polynomial in row.items():
        target = result.setdefault(tuple(wedge), {})
        for source_exponent, matrix in polynomial.items():
            shifted = tuple(
                int(a) + int(b) for a, b in zip(source_exponent, exponent)
            )
            _add_matrix_term(target, shifted, matrix, scale)
        if not target:
            result.pop(tuple(wedge), None)
    return result


def _assert_rows_equal(left, right, label):
    wedges = set(left) | set(right)
    for wedge in wedges:
        left_poly = left.get(wedge, {})
        right_poly = right.get(wedge, {})
        for exponent in set(left_poly) | set(right_poly):
            left_matrix = left_poly.get(exponent)
            right_matrix = right_poly.get(exponent)
            if left_matrix is None:
                left_matrix = np.zeros_like(right_matrix)
            if right_matrix is None:
                right_matrix = np.zeros_like(left_matrix)
            if not np.array_equal(left_matrix % P, right_matrix % P):
                raise HSLambdaTransportError(
                    f"{label} failed at wedge {wedge}, exponent {exponent}"
                )


def _row_term_count(row):
    return sum(len(polynomial) for polynomial in row.values())


def _row_scalar_count(row):
    return sum(
        int(np.count_nonzero(matrix % P))
        for polynomial in row.values()
        for matrix in polynomial.values()
    )


def _row_hash(row):
    digest = hashlib.sha256()
    for wedge in sorted(row):
        digest.update(json.dumps(wedge, separators=(",", ":")).encode("ascii"))
        for exponent in sorted(row[wedge]):
            digest.update(
                json.dumps(exponent, separators=(",", ":")).encode("ascii")
            )
            digest.update(
                np.ascontiguousarray(row[wedge][exponent] % P).tobytes()
            )
    return digest.hexdigest()


def _lambda_row(equations, values, tangents, residue_row):
    """Build ``-tau H5...H0`` as a row on graph K2."""

    _multiplication, operators = hs_extension.build_local_h_operators(
        equations, values, tangents
    )
    cochain = hs_extension._local_change_of_rings_cochain(
        operators, residue_row, range(6)
    )
    pairs = tuple(itertools.combinations(range(VARIABLE_COUNT), 2))
    result = {}
    for column, polynomial in cochain.items():
        pair_number, local_source = divmod(int(column), 2)
        wedge = pairs[pair_number]
        target = result.setdefault(wedge, {})
        for exponent, coefficient in polynomial.items():
            matrix = np.zeros((1, 2), dtype=np.int64)
            matrix[0, local_source] = int(coefficient) % P
            _add_matrix_term(target, exponent, matrix)
    return result


def _evaluate_row_polynomial(polynomial, multiplication, inverse_variable):
    """Evaluate a Laurent-polynomial functional on one Hermite factor."""

    result = np.zeros((1, 2), dtype=np.int64)
    inverse = shamash.inverse_2x2(multiplication[inverse_variable])
    for exponent, row_matrix in polynomial.items():
        operator = np.eye(2, dtype=np.int64)
        for variable, power in enumerate(exponent):
            operator = (
                operator
                @ shamash.matrix_power(
                    multiplication[variable],
                    power,
                    inverse if variable == inverse_variable else None,
                )
            ) % P
        result = (result + np.asarray(row_matrix) @ operator) % P
    return result


def _evaluate_block_polynomial(polynomial, multiplication, inverse_variable):
    """Evaluate a Laurent-polynomial A-linear endomorphism block."""

    result = np.zeros((2, 2), dtype=np.int64)
    inverse = shamash.inverse_2x2(multiplication[inverse_variable])
    for exponent, coefficient_matrix in polynomial.items():
        operator = np.eye(2, dtype=np.int64)
        for variable, power in enumerate(exponent):
            operator = (
                operator
                @ shamash.matrix_power(
                    multiplication[variable],
                    power,
                    inverse if variable == inverse_variable else None,
                )
            ) % P
        result = (result + operator @ np.asarray(coefficient_matrix)) % P
    return result


def _rref_solve(system, right_hand_side):
    """Return the canonical free-variables-zero solution over F_31."""

    system = np.asarray(system, dtype=np.int64) % P
    right_hand_side = np.asarray(right_hand_side, dtype=np.int64).reshape(-1) % P
    augmented = np.column_stack((system, right_hand_side)) % P
    row = 0
    pivots = []
    for column in range(system.shape[1]):
        candidates = np.flatnonzero(augmented[row:, column] % P)
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        if pivot != row:
            augmented[[row, pivot]] = augmented[[pivot, row]]
        augmented[row] = (
            augmented[row] * pow(int(augmented[row, column]), -1, P)
        ) % P
        for other in range(augmented.shape[0]):
            if other == row:
                continue
            scale = int(augmented[other, column]) % P
            if scale:
                augmented[other] = (augmented[other] - scale * augmented[row]) % P
        pivots.append(column)
        row += 1
        if row == augmented.shape[0]:
            break
    for remaining in range(row, augmented.shape[0]):
        if not np.any(augmented[remaining, :-1] % P) and augmented[remaining, -1] % P:
            raise HSLambdaTransportError("the 68-by-16 Hermite system is inconsistent")
    solution = np.zeros(system.shape[1], dtype=np.int64)
    for pivot_row, pivot_column in enumerate(pivots):
        solution[pivot_column] = augmented[pivot_row, -1]
    if np.any((system @ solution - right_hand_side) % P):
        raise HSLambdaTransportError("the canonical RREF solution did not verify")
    free = [column for column in range(system.shape[1]) if column not in pivots]
    return solution, {
        "rank": len(pivots),
        "nullity": system.shape[1] - len(pivots),
        "pivot_columns": pivots,
        "free_columns": free,
        "system_sha256": hashlib.sha256(
            np.ascontiguousarray(system).tobytes()
        ).hexdigest(),
        "rhs_sha256": hashlib.sha256(
            np.ascontiguousarray(right_hand_side).tobytes()
        ).hexdigest(),
        "rref_sha256": hashlib.sha256(
            np.ascontiguousarray(augmented).tobytes()
        ).hexdigest(),
    }


def _solve_hermite_correction(
    target_d2,
    target_h0,
    delta_target,
    delta_equations_target,
    target_multiplication,
    inverse_variable,
):
    """Solve one A-linear row equation across all 34 target columns."""

    pairs = tuple(itertools.combinations(range(VARIABLE_COUNT), 2))
    columns = [
        (f"K2_{index}", wedge, delta_target.get(wedge, {}), target_d2)
        for index, wedge in enumerate(pairs)
    ]
    columns.extend(
        (
            f"H0_{equation}",
            (),
            delta_equations_target[equation].get((), {}),
            target_h0[equation],
        )
        for equation in range(6)
    )
    system = np.zeros((2 * len(columns), 2 * VARIABLE_COUNT), dtype=np.int64)
    right_hand_side = np.zeros(2 * len(columns), dtype=np.int64)
    column_labels = []
    for column, (label, input_wedge, delta_polynomial, blocks) in enumerate(columns):
        column_labels.append(label)
        right_hand_side[2 * column : 2 * column + 2] = _evaluate_row_polynomial(
            delta_polynomial, target_multiplication, inverse_variable
        )[0]
        for variable in range(VARIABLE_COUNT):
            polynomial = blocks.get(((variable,), input_wedge), {})
            if not polynomial:
                continue
            block = _evaluate_block_polynomial(
                polynomial, target_multiplication, inverse_variable
            )
            for source_coordinate in range(2):
                for target_coordinate in range(2):
                    system[
                        2 * column + target_coordinate,
                        2 * variable + source_coordinate,
                    ] = block[source_coordinate, target_coordinate] % P
    solution, ledger = _rref_solve(system, right_hand_side)
    a_target = {}
    for variable in range(VARIABLE_COUNT):
        row = solution[2 * variable : 2 * variable + 2].reshape(1, 2)
        if np.any(row % P):
            a_target[(variable,)] = {ZERO_EXPONENT: row % P}
    ledger.update(
        {
            "base_algebra": "A=F_31[epsilon]/(epsilon^2)",
            "matrix_shape": [int(value) for value in system.shape],
            "A_unknowns": VARIABLE_COUNT,
            "F31_unknowns": 2 * VARIABLE_COUNT,
            "A_target_columns": len(columns),
            "F31_scalar_equations": 2 * len(columns),
            "column_labels": column_labels,
            "canonical_free_variables_zero": True,
            "solution_sha256": hashlib.sha256(
                np.ascontiguousarray(solution).tobytes()
            ).hexdigest(),
        }
    )
    return a_target, ledger


def _verify_hermite_residual(
    residual_k2,
    residual_equations,
    target_multiplication,
    inverse_variable,
):
    pairs = tuple(itertools.combinations(range(VARIABLE_COUNT), 2))
    labels = []
    polynomial_nonzero = []
    for index, wedge in enumerate(pairs):
        polynomial = residual_k2.get(wedge, {})
        labels.append(f"K2_{index}")
        if polynomial:
            polynomial_nonzero.append(f"K2_{index}")
        if np.any(
            _evaluate_row_polynomial(
                polynomial, target_multiplication, inverse_variable
            )
        ):
            raise HSLambdaTransportError(f"the K2_{index} Hermite residual is nonzero")
    for equation in range(6):
        polynomial = residual_equations[equation].get((), {})
        labels.append(f"H0_{equation}")
        if polynomial:
            polynomial_nonzero.append(f"H0_{equation}")
        if np.any(
            _evaluate_row_polynomial(
                polynomial, target_multiplication, inverse_variable
            )
        ):
            raise HSLambdaTransportError(
                f"the H0_{equation} Hermite residual is nonzero"
            )
    return {
        "verified_A_columns": len(labels),
        "verified_F31_scalar_coordinates": 2 * len(labels),
        "column_labels": labels,
        "polynomially_nonzero_before_Hermite_reduction": polynomial_nonzero,
        "all_Hermite_residuals_zero": True,
    }


def _boundary_delta(context, surface, chart, values, tangents):
    _labels, delta, _inverse = hs_frames.boundary_delta_on_chart(
        context, surface, chart, values, tangents
    )
    epsilon = chart_parity(chart)
    residue = [epsilon * int(delta[1]) % P, epsilon * int(delta[0]) % P]
    return delta, residue, epsilon


def _verify_orientation_weight(
    source_delta,
    target_delta,
    source_residue,
    target_residue,
    source_multiplication,
    changed,
):
    v = 2 * changed + 1
    mv = source_multiplication[v]
    mt = shamash.inverse_2x2(mv)
    source_delta_row = np.asarray(
        [source_delta[1], source_delta[0]], dtype=np.int64
    ) % P
    target_delta_row = np.asarray(
        [target_delta[1], target_delta[0]], dtype=np.int64
    ) % P
    if not np.array_equal(
        target_delta_row,
        source_delta_row @ shamash.matrix_power(mt, 2) % P,
    ):
        raise HSLambdaTransportError("phi(delta_target) is not t^2 delta_source")
    j8 = -shamash.matrix_power(mt, 3) % P
    if not np.array_equal(
        np.asarray(source_residue) @ j8 % P,
        np.asarray(target_residue) @ mt % P,
    ):
        raise HSLambdaTransportError("the normalized tau/J8 weight is not +t")
    return {
        "phi_delta_target": "t^2 delta_source",
        "normalized_tau_source_J8": "t phi(tau_target)",
        "det_B": "t^5",
        "O4H_factor": "t^-4=v_source^4",
    }


def build_support(data, equations, surface, edge_id, support_position):
    overlap = data["overlaps"][(surface, edge_id)]
    source_chart = data["charts"][(surface, overlap["source"])]
    target_chart = data["charts"][(surface, overlap["target"])]
    source_position = int(overlap["source_positions"][support_position])
    target_position = int(overlap["target_positions"][support_position])
    changed = int(overlap["changed"])
    v = 2 * changed + 1

    source_values = source_chart["values"][source_position]
    source_tangents = source_chart["tangents"][source_position]
    target_values = target_chart["values"][target_position]
    target_tangents = target_chart["tangents"][target_position]
    source_equations = equations[(surface, overlap["source"])]
    target_equations = equations[(surface, overlap["target"])]

    source_delta, source_residue, epsilon_source = _boundary_delta(
        data["context"], surface, overlap["source"], source_values, source_tangents
    )
    target_delta, target_residue, epsilon_target = _boundary_delta(
        data["context"], surface, overlap["target"], target_values, target_tangents
    )
    if epsilon_target != -epsilon_source:
        raise HSLambdaTransportError("an oriented edge did not flip chart parity")

    source_multiplication = shamash.local_multiplication(
        source_values, source_tangents
    )
    target_multiplication = shamash.local_multiplication(
        target_values, target_tangents
    )
    orientation = _verify_orientation_weight(
        source_delta,
        target_delta,
        source_residue,
        target_residue,
        source_multiplication,
        changed,
    )

    comparison = shamash.local_support_maps(
        source_values,
        source_tangents,
        target_values,
        target_tangents,
        source_equations,
        target_equations,
        changed,
    )
    lambda_source = _lambda_row(
        source_equations, source_values, source_tangents, source_residue
    )
    lambda_target = _lambda_row(
        target_equations, target_values, target_tangents, target_residue
    )

    source_composed = _compose_row_with_blocks(
        lambda_source, comparison["graph_maps"][2]
    )
    target_pulled = _transform_row_exponents(lambda_target, changed)
    h4_exponent = list(ZERO_EXPONENT)
    h4_exponent[v] = 4
    target_weighted = _shift_row(target_pulled, tuple(h4_exponent))
    delta_k2 = dict(source_composed)
    _add_row_map(delta_k2, target_weighted, -1)

    delta_equations = [
        _compose_row_with_blocks(lambda_source, gamma)
        for gamma in comparison["gamma"]
    ]

    # The coordinate exponent map is an involution.  The upstream J2/J3
    # comparison is supportwise over the two-dimensional Hermite algebra, so
    # solve the complete A-linear equation there.  This is one 68-by-16
    # F_31 system: all 28 K2 columns and all six H0 columns constrain one row.
    delta_target = _transform_row_exponents(delta_k2, changed)
    target_d2 = shamash.graph_differential(2, target_multiplication)
    target_h0 = []
    delta_equations_target = [
        _transform_row_exponents(row, changed) for row in delta_equations
    ]
    for equation, target_polynomial in enumerate(target_equations):
        target_h0.append(
            shamash.divided_difference_matrix(
                target_polynomial, target_multiplication
            )
        )
    a_target, linear_system = _solve_hermite_correction(
        target_d2,
        target_h0,
        delta_target,
        delta_equations_target,
        target_multiplication,
        v,
    )
    residual_k2 = _compose_row_with_blocks(a_target, target_d2)
    _add_row_map(residual_k2, delta_target, -1)
    residual_equations = []
    for equation in range(6):
        residual = _compose_row_with_blocks(a_target, target_h0[equation])
        _add_row_map(residual, delta_equations_target[equation], -1)
        residual_equations.append(residual)
    hermite_check = _verify_hermite_residual(
        residual_k2,
        residual_equations,
        target_multiplication,
        v,
    )
    a_source_coordinates = _transform_row_exponents(a_target, changed)
    return {
        "support_position": int(support_position),
        "source_position": source_position,
        "target_position": target_position,
        "epsilon_source": epsilon_source,
        "epsilon_target": epsilon_target,
        "orientation": orientation,
        "lambda_source": lambda_source,
        "lambda_target": lambda_target,
        "delta_k2": delta_k2,
        "delta_equations": delta_equations,
        "a_target": a_target,
        "a_source_coordinates": a_source_coordinates,
        "linear_system": linear_system,
        "hermite_residual_ledger": hermite_check,
        "checks": {
            "single_a_solves_K2_block_in_Hermite_quotient": True,
            "same_a_solves_all_six_H0_blocks_in_Hermite_quotient": True,
            "lambda_Gamma_blocks_included": True,
        },
    }


def load_equations(data):
    result = {}
    for surface, chart in data["charts"]:
        result[(surface, chart)] = shamash.surface_equations(
            data["context"], surface, chart
        )[0]
    return result


def build_edge(surface=DEFAULT_SURFACE, edge_id=DEFAULT_EDGE, support=None):
    data = shamash.load_problem_data()
    key = (int(surface), int(edge_id))
    if key not in data["overlaps"]:
        raise HSLambdaTransportError(f"edge {key} is absent")
    equations = load_equations(data)
    overlap = data["overlaps"][key]
    positions = (
        [int(support)]
        if support is not None
        else list(range(len(overlap["supports"])))
    )
    if any(position not in range(len(overlap["supports"])) for position in positions):
        raise HSLambdaTransportError("a support position is outside the edge")
    rows = [
        build_support(data, equations, key[0], key[1], position)
        for position in positions
    ]
    return {
        "field": P,
        "surface": key[0],
        "edge_id": key[1],
        "source_chart": overlap["source"],
        "target_chart": overlap["target"],
        "changed_factor": int(overlap["changed"]),
        "edge_support_count": len(overlap["supports"]),
        "computed_support_count": len(rows),
        "supports": rows,
        "input_paths": data["input_paths"],
    }


def _serialize_row(rows, name):
    wedge_indices = {
        1: {(variable,): variable for variable in range(VARIABLE_COUNT)},
        2: {
            wedge: index
            for index, wedge in enumerate(
                itertools.combinations(range(VARIABLE_COUNT), 2)
            )
        },
    }
    records = []
    for support_record in rows:
        row = support_record[name]
        for wedge, polynomial in sorted(row.items()):
            for exponent, matrix in sorted(polynomial.items()):
                for local_column in range(matrix.shape[1]):
                    coefficient = int(matrix[0, local_column]) % P
                    if coefficient:
                        records.append(
                            [
                                support_record["support_position"],
                                len(wedge),
                                wedge_indices.get(len(wedge), {}).get(
                                    tuple(wedge), -1
                                ),
                                *map(int, exponent),
                                local_column,
                                coefficient,
                            ]
                        )
    return np.asarray(records, dtype=np.int16)


def write_artifacts(result, output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    stem = f"Z_{result['surface']}_E{result['edge_id']:02d}"
    payload = {
        "field": np.asarray([P], dtype=np.uint8),
        "support_positions": np.asarray(
            [row["support_position"] for row in result["supports"]], dtype=np.uint16
        ),
        "a_target_source_coordinate_terms": _serialize_row(
            result["supports"], "a_source_coordinates"
        ),
        "delta_K2_target_terms": _serialize_row(
            [
                {
                    "support_position": row["support_position"],
                    "value": _transform_row_exponents(
                        row["delta_k2"], result["changed_factor"]
                    ),
                }
                for row in result["supports"]
            ],
            "value",
        ),
    }
    for equation in range(6):
        records = []
        for row in result["supports"]:
            temporary = {"support_position": row["support_position"], "value": row["delta_equations"][equation]}
            records.append(temporary)
        payload[f"delta_equation_{equation}_terms"] = _serialize_row(
            records, "value"
        )
    binary = output_directory / f"{stem}_lambda_transport.npz"
    np.savez_compressed(binary, **payload)

    support_summaries = []
    for row in result["supports"]:
        support_summaries.append(
            {
                "support_position": row["support_position"],
                "source_position": row["source_position"],
                "target_position": row["target_position"],
                "epsilon_source": row["epsilon_source"],
                "epsilon_target": row["epsilon_target"],
                "delta_K2_terms": _row_term_count(row["delta_k2"]),
                "delta_equation_terms": [
                    _row_term_count(value) for value in row["delta_equations"]
                ],
                "a_terms": _row_term_count(row["a_source_coordinates"]),
                "a_scalar_nonzeros": _row_scalar_count(row["a_source_coordinates"]),
                "a_sha256": _row_hash(row["a_source_coordinates"]),
                "linear_system": row["linear_system"],
                "hermite_residual_ledger": row["hermite_residual_ledger"],
                **row["checks"],
            }
        )
    certificate = {
        "schema": "section12-HS-lambda-edge-transport-pilot-v1",
        "status": "exact-single-a-supportwise-Hermite-counit-transport",
        "field": "F_31",
        "surface": f"Z_{result['surface']}",
        "edge_id": result["edge_id"],
        "source_chart": result["source_chart"],
        "target_chart": result["target_chart"],
        "changed_factor": result["changed_factor"],
        "counts": {
            "edge_support_count": result["edge_support_count"],
            "computed_support_count": result["computed_support_count"],
        },
        "orientation_normalization": {
            "epsilon_C": "(-1)^popcount(chart_bits)",
            "raw_tau_relation": "tau_s_raw J8=-t phi(tau_t_raw)",
            "normalized_tau_relation": "tau_s J8=t phi(tau_t)",
            "B": "diag(t^3,1,1,1,t,t)",
            "det_B": "t^5",
            "h": "v_source=t^-1",
            "O4H_factor": "t/det(B)=t^-4=h^4",
        },
        "typed_equation": {
            "Delta": "Lambda_edge J2-h^4 phi(Lambda_target)",
            "Delta_K2": "lambda_edge J2K-h^4 phi(lambda_target)",
            "Delta_j": "lambda_edge Gamma_j",
            "single_correction": "a_target:K1_target->O(4H)_edge",
            "identity": "Delta=a_target D2_target",
            "induced_pushout_map": "[[h^4,-a_target],[0,J1]]",
        },
        "binary_term_schema": [
            "support_position",
            "wedge_degree",
            "wedge_index",
            "exponent_x0",
            "exponent_x1",
            "exponent_x2",
            "exponent_x3",
            "exponent_x4",
            "exponent_x5",
            "exponent_x6",
            "exponent_x7",
            "local_column",
            "coefficient_mod_31",
        ],
        "support_ledger": support_summaries,
        "binary_artifact": {
            "path": binary.name,
            "bytes": binary.stat().st_size,
            "sha256": sha256(binary),
        },
        "upstream": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in result["input_paths"]
            if path.exists()
        },
        "degree_window": (
            "P0-through-P3 suffices for edgewise Ext2 closure. Face/cube "
            "coherence, especially a degree-3 square defect, may require P4 "
            "unless that defect is strictly zero."
        ),
        "claim_boundary": {
            "proved": (
                "the displayed edge/support Hermite factors admit one exact "
                "degree-one cochain transporting the chart-local HS counit with "
                "O(4H) weight; one 68-by-16 system enforces every K2 and H0 column"
            ),
            "not_claimed": (
                "a polynomial or whole-surface quotient lift, all 192 edges, "
                "face/cube coherence, K-translated descent, a global HS object, "
                "Yoneda/T, Question 11.4, characteristic zero, or Hodge"
            ),
        },
    }
    certificate_path = output_directory / f"{stem}_lambda_transport_certificate.json"
    certificate_path.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest = {
        "schema": "section12-HS-lambda-edge-transport-pilot-manifest-v1",
        "producer": str(Path(__file__).resolve().relative_to(ROOT)),
        "producer_sha256": sha256(Path(__file__).resolve()),
        "artifacts": {
            certificate_path.name: sha256(certificate_path),
            binary.name: sha256(binary),
        },
    }
    manifest_path = output_directory / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return certificate_path, binary, manifest_path


def verify_artifacts(output_directory=OUTPUT_DIRECTORY, surface=DEFAULT_SURFACE, edge_id=DEFAULT_EDGE):
    output_directory = Path(output_directory)
    stem = f"Z_{int(surface)}_E{int(edge_id):02d}"
    certificate_path = output_directory / f"{stem}_lambda_transport_certificate.json"
    manifest_path = output_directory / "manifest.json"
    certificate = json.loads(certificate_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if certificate.get("schema") != "section12-HS-lambda-edge-transport-pilot-v1":
        raise HSLambdaTransportError("the certificate schema changed")
    if not all(
        row["single_a_solves_K2_block_in_Hermite_quotient"]
        and row["same_a_solves_all_six_H0_blocks_in_Hermite_quotient"]
        and row["hermite_residual_ledger"]["all_Hermite_residuals_zero"]
        for row in certificate["support_ledger"]
    ):
        raise HSLambdaTransportError("a persisted single-a check is false")
    for name, digest in manifest["artifacts"].items():
        if sha256(output_directory / name) != digest:
            raise HSLambdaTransportError(f"artifact hash mismatch: {name}")
    return certificate


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--surface", type=int, default=DEFAULT_SURFACE)
    parser.add_argument("--edge", type=int, default=DEFAULT_EDGE)
    parser.add_argument("--support", type=int)
    parser.add_argument("--verify", action="store_true")
    arguments = parser.parse_args()
    if arguments.verify:
        certificate = verify_artifacts(
            surface=arguments.surface, edge_id=arguments.edge
        )
        print(json.dumps(certificate["counts"], sort_keys=True))
        return
    result = build_edge(arguments.surface, arguments.edge, arguments.support)
    certificate, binary, manifest = write_artifacts(result)
    print(certificate)
    print(binary)
    print(manifest)


if __name__ == "__main__":
    main()
