"""Close the six-sector Shamash comparison for one nonidentity K generator.

The comparison uses one functorial law for all six equations.  On each
nontrivial factor, source affine coordinates are changed linearly to ``(a,b)``
with ``b`` equal to the pulled target frame.  Every target coordinate then has
the form ``N(a,b)/b``.  The comparison therefore lives in the exact Laurent
ring obtained by inverting the two nonconstant frame variables.

The producer constructs the graph compounds, the diagonal regular-sequence
change B, and the common corrections

    Delta_j = J1 H'_j - H_j B_j J0,
    Gamma_j = -h(Delta_j).

It checks ``d2 Gamma_j=Delta_j`` and the complete degree-two Shamash identity.
It also constructs the degree-three Omega inputs and checks their identities.
Only after J2 closes does it test the Hartshorne--Serre counit.  The single
row ``a`` is solved against all 28 graph columns and all six equation columns
at once on every retained dual-number factor.  If that system is inconsistent,
the emitted obstruction is an exact left-nullspace witness.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import os
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_all_chart_presentations as all_charts,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_lambda_overlap_transport as lambda_transport,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_rank2_k_chain_comparison as k_chain,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_rank2_k_generator_overlap_j1 as k_j1,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_shamash_overlap_lifts as shamash,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_iw_hermite_border_basis as hermite,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_iw_presentation_solver as iw_solver,
)


P = 31
VARIABLE_COUNT = 8
EQUATION_COUNT = 6
SURFACE = 2
CHART = "C0010"
GENERATOR_EXPONENTS = (1, 0, 0, 0)
ZERO = (0,) * VARIABLE_COUNT
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_rank2_k_generator_shamash"
)


class HSKGeneratorShamashError(RuntimeError):
    """A frame, Shamash, or counit identity changed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def array_sha256(value) -> str:
    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(value)).tobytes()
    ).hexdigest()


def pair_ledger_sha256(value) -> str:
    """Hash dual-number pair ledgers in one platform-independent encoding."""

    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(value, dtype="<i8")).tobytes()
    ).hexdigest()


def write_json_atomic(path: Path, value) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def write_csv(path: Path, rows) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError("cannot write an empty CSV")
    fieldnames = []
    for row in rows:
        for name in row:
            if name not in fieldnames:
                fieldnames.append(name)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def add_scalar_term(polynomial, exponent, coefficient):
    exponent = tuple(map(int, exponent))
    value = (int(polynomial.get(exponent, 0)) + int(coefficient)) % P
    if value:
        polynomial[exponent] = value
    else:
        polynomial.pop(exponent, None)


def scalar_add(left, right, scale=1):
    result = dict(left)
    for exponent, coefficient in right.items():
        add_scalar_term(result, exponent, int(scale) * int(coefficient))
    return result


def scalar_multiply(left, right):
    result = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                int(a) + int(b)
                for a, b in zip(left_exponent, right_exponent)
            )
            add_scalar_term(
                result,
                exponent,
                int(left_coefficient) * int(right_coefficient),
            )
    return result


def scalar_power(polynomial, exponent):
    exponent = int(exponent)
    if exponent < 0:
        raise ValueError("scalar polynomial powers must be nonnegative")
    result = {ZERO: 1}
    factor = dict(polynomial)
    while exponent:
        if exponent & 1:
            result = scalar_multiply(result, factor)
        factor = scalar_multiply(factor, factor)
        exponent //= 2
    return result


def linear_polynomial(form, factor):
    result = {}
    add_scalar_term(result, ZERO, form[0])
    for local in range(2):
        exponent = list(ZERO)
        exponent[2 * factor + local] = 1
        add_scalar_term(result, exponent, form[local + 1])
    return result


def substitute_scalar_polynomial(polynomial, substitutions):
    result = {}
    for exponent, coefficient in polynomial.items():
        term = {ZERO: int(coefficient) % P}
        for variable, power in enumerate(exponent):
            if power:
                term = scalar_multiply(
                    term, scalar_power(substitutions[variable], int(power))
                )
        result = scalar_add(result, term)
    return result


def substitute_block_map(blocks, substitutions):
    result = {}
    for (output_wedge, input_wedge), polynomial in blocks.items():
        target = {}
        for exponent, matrix in polynomial.items():
            scalar = substitute_scalar_polynomial({exponent: 1}, substitutions)
            for target_exponent, coefficient in scalar.items():
                shamash._add_poly_term(
                    target,
                    target_exponent,
                    int(coefficient) * np.asarray(matrix, dtype=np.int64),
                )
        shamash._add_block(result, output_wedge, input_wedge, target)
    return result


def substitute_row(row, substitutions):
    result = {}
    for wedge, polynomial in row.items():
        target = result.setdefault(tuple(wedge), {})
        for exponent, matrix in polynomial.items():
            scalar = substitute_scalar_polynomial({exponent: 1}, substitutions)
            for target_exponent, coefficient in scalar.items():
                lambda_transport._add_matrix_term(
                    target,
                    target_exponent,
                    int(coefficient) * np.asarray(matrix, dtype=np.int64),
                )
        if not target:
            result.pop(tuple(wedge), None)
    return result


def determinant_mod(matrix):
    matrix = np.asarray(matrix, dtype=np.int64) % P
    if matrix.shape != (3, 3):
        raise ValueError("the affine-frame determinant expects a 3-by-3 matrix")
    a, b, c = map(int, matrix[0])
    d, e, f = map(int, matrix[1])
    g, h, i = map(int, matrix[2])
    return (
        a * (e * i - f * h)
        - b * (d * i - f * g)
        + c * (d * h - e * g)
    ) % P


def adapted_factor(frame, numerator_u, numerator_v, factor):
    """Choose (a,b) with b=L and return exact affine substitutions."""

    frame = tuple(map(int, frame))
    numerator_u = tuple(map(int, numerator_u))
    numerator_v = tuple(map(int, numerator_v))
    identity_frame = frame == (1, 0, 0)
    if identity_frame:
        transform = np.eye(3, dtype=np.int64)
    else:
        candidates = (numerator_u, numerator_v, (0, 1, 0), (0, 0, 1))
        transform = None
        for candidate in candidates:
            trial = np.asarray(((1, 0, 0), candidate, frame), dtype=np.int64) % P
            try:
                hermite.inverse_mod(trial)
            except Exception:
                continue
            transform = trial
            break
        if transform is None:
            raise HSKGeneratorShamashError("no affine coordinate complements the frame")
    inverse = hermite.inverse_mod(transform)
    if not np.array_equal(transform @ inverse % P, np.eye(3, dtype=np.int64)):
        raise HSKGeneratorShamashError("an adapted affine frame is not invertible")
    u_adapted = tuple(map(int, np.asarray(numerator_u) @ inverse % P))
    v_adapted = tuple(map(int, np.asarray(numerator_v) @ inverse % P))
    frame_adapted = tuple(map(int, np.asarray(frame) @ inverse % P))
    expected_frame = (1, 0, 0) if identity_frame else (0, 0, 1)
    if frame_adapted != expected_frame:
        raise HSKGeneratorShamashError("the adapted denominator is not the b coordinate")

    original_u = tuple(map(int, inverse[1]))
    original_v = tuple(map(int, inverse[2]))
    original_substitutions = (
        linear_polynomial(original_u, factor),
        linear_polynomial(original_v, factor),
    )

    def rational_coordinate(form):
        if identity_frame:
            return linear_polynomial(form, factor)
        result = {}
        for index, coefficient in enumerate(form):
            if not coefficient:
                continue
            exponent = list(ZERO)
            exponent[2 * factor + 1] = -1
            if index == 1:
                exponent[2 * factor] += 1
            elif index == 2:
                exponent[2 * factor + 1] += 1
            add_scalar_term(result, exponent, coefficient)
        return result

    return {
        "factor": factor,
        "identity": identity_frame,
        "transform": transform,
        "inverse": inverse,
        "frame_original": frame,
        "numerator_u_original": numerator_u,
        "numerator_v_original": numerator_v,
        "frame_adapted": frame_adapted,
        "numerator_u_adapted": u_adapted,
        "numerator_v_adapted": v_adapted,
        "original_substitutions": original_substitutions,
        "target_substitutions": (
            rational_coordinate(u_adapted),
            rational_coordinate(v_adapted),
        ),
        "inverse_variable": None if identity_frame else 2 * factor + 1,
    }


def adapt_jet(values, tangents, factor_data):
    adapted_values = []
    adapted_tangents = []
    for factor, data in enumerate(factor_data):
        u, v = map(int, values[2 * factor : 2 * factor + 2])
        du, dv = map(int, tangents[2 * factor : 2 * factor + 2])
        vector = np.asarray((1, u, v), dtype=np.int64)
        derivative = np.asarray((0, du, dv), dtype=np.int64)
        image = data["transform"] @ vector % P
        image_derivative = data["transform"] @ derivative % P
        if int(image[0]) != 1:
            raise HSKGeneratorShamashError("an adapted affine constant changed")
        adapted_values.extend((int(image[1]), int(image[2])))
        adapted_tangents.extend((int(image_derivative[1]), int(image_derivative[2])))
    return tuple(adapted_values), tuple(adapted_tangents)


def graph_compounds_general(target_multiplication, factor_data):
    identity = np.eye(2, dtype=np.int64)
    generator_terms = {}
    for target_variable in range(VARIABLE_COUNT):
        factor, local = divmod(target_variable, 2)
        data = factor_data[factor]
        numerator = data[
            "numerator_u_adapted" if local == 0 else "numerator_v_adapted"
        ]
        terms = []
        if data["identity"]:
            for source_local in range(2):
                coefficient = int(numerator[source_local + 1]) * identity % P
                if np.any(coefficient):
                    terms.append((2 * factor + source_local, ZERO, coefficient))
        else:
            exponent = list(ZERO)
            exponent[2 * factor + 1] = -1
            exponent = tuple(exponent)
            coefficient_a = int(numerator[1]) * identity % P
            coefficient_b = (
                int(numerator[2]) * identity
                - target_multiplication[target_variable]
            ) % P
            if np.any(coefficient_a):
                terms.append((2 * factor, exponent, coefficient_a))
            if np.any(coefficient_b):
                terms.append((2 * factor + 1, exponent, coefficient_b))
        if not terms:
            raise HSKGeneratorShamashError("a target graph generator mapped to zero")
        generator_terms[target_variable] = terms

    maps = {0: {((), ()): {ZERO: identity}}}
    for degree in (1, 2, 3, 8):
        result = {}
        source_wedges = (
            shamash.WEDGES[degree]
            if degree in shamash.WEDGES
            else (tuple(range(VARIABLE_COUNT)),)
        )
        for source_wedge in source_wedges:
            for selected in itertools.product(
                *(generator_terms[variable] for variable in source_wedge)
            ):
                output, sign = shamash._wedge_sort_sign(
                    tuple(term[0] for term in selected)
                )
                if not sign:
                    continue
                exponent = ZERO
                coefficient = identity
                for _output, term_exponent, matrix in selected:
                    exponent = shamash._sum_exponents(exponent, term_exponent)
                    coefficient = coefficient @ matrix % P
                shamash._add_block(
                    result,
                    output,
                    source_wedge,
                    {exponent: coefficient},
                    scale=sign,
                )
        maps[degree] = result
    return maps


def matrix_power_multi(multiplication, variable, power, inverses):
    return shamash.matrix_power(
        multiplication[variable], int(power), inverses.get(variable)
    )


def ordered_laurent_divided_multi(polynomial, variable, multiplication, inverses):
    result = {}
    for exponent, coefficient_matrix in polynomial.items():
        exponent = tuple(map(int, exponent))
        power = exponent[variable]
        if not power:
            continue
        prefix = np.eye(2, dtype=np.int64)
        for earlier in range(variable):
            prefix = (
                prefix
                @ matrix_power_multi(
                    multiplication, earlier, exponent[earlier], inverses
                )
            ) % P
        if power > 0:
            summands = [(power - 1 - r, r, 1) for r in range(power)]
        else:
            summands = [
                (-r, power + r - 1, -1) for r in range(1, -power + 1)
            ]
        for polynomial_power, matrix_exponent, sign in summands:
            output_exponent = [0] * VARIABLE_COUNT
            output_exponent[variable] = polynomial_power
            for later in range(variable + 1, VARIABLE_COUNT):
                output_exponent[later] = exponent[later]
            matrix = (
                prefix
                @ matrix_power_multi(
                    multiplication, variable, matrix_exponent, inverses
                )
                @ coefficient_matrix
            ) % P
            shamash._add_poly_term(
                result, tuple(output_exponent), matrix, scale=sign
            )
    return result


def tensor_ordered_homotopy_multi(blocks, degree, multiplication, inverses):
    result = {}
    for (output_wedge, input_wedge), polynomial in blocks.items():
        if len(output_wedge) != degree:
            raise HSKGeneratorShamashError("a homotopy input has the wrong degree")
        for variable in range(min(output_wedge)):
            divided = ordered_laurent_divided_multi(
                polynomial, variable, multiplication, inverses
            )
            shamash._add_block(
                result, (variable, *output_wedge), input_wedge, divided
            )
    return result


def evaluate_row_multi(polynomial, multiplication, inverses):
    result = np.zeros((1, 2), dtype=np.int64)
    for exponent, row_matrix in polynomial.items():
        operator = np.eye(2, dtype=np.int64)
        for variable, power in enumerate(exponent):
            operator = (
                operator
                @ matrix_power_multi(multiplication, variable, power, inverses)
            ) % P
        result = (result + np.asarray(row_matrix) @ operator) % P
    return result


def evaluate_block_multi(polynomial, multiplication, inverses):
    result = np.zeros((2, 2), dtype=np.int64)
    for exponent, coefficient_matrix in polynomial.items():
        operator = np.eye(2, dtype=np.int64)
        for variable, power in enumerate(exponent):
            operator = (
                operator
                @ matrix_power_multi(multiplication, variable, power, inverses)
            ) % P
        result = (result + operator @ np.asarray(coefficient_matrix)) % P
    return result


def dual_matrix(pair):
    value, tangent = map(lambda item: int(item) % P, pair)
    return np.asarray([[value, 0], [tangent, value]], dtype=np.int64)


def dual_multiply(left, right):
    left_value, left_tangent = map(lambda item: int(item) % P, left)
    right_value, right_tangent = map(lambda item: int(item) % P, right)
    return (
        left_value * right_value % P,
        (left_value * right_tangent + left_tangent * right_value) % P,
    )


def dual_inverse(pair):
    value, tangent = map(lambda item: int(item) % P, pair)
    if not value:
        raise HSKGeneratorShamashError(
            "a residue/counit transport factor is not a unit"
        )
    inverse_value = pow(value, -1, P)
    return (
        inverse_value,
        (-tangent * inverse_value * inverse_value) % P,
    )


def dual_divide(numerator, denominator):
    return dual_multiply(numerator, dual_inverse(denominator))


def scalar_laurent_pair(scalar, exponent, values, tangents):
    value, tangent = hermite.monomial_jet(
        tuple(map(int, exponent)), values, tangents
    )
    return int(scalar) * int(value) % P, int(scalar) * int(tangent) % P


def frobenius_multiplier(target_residue, transported_residue):
    """Recover the unique local multiplication element between residues."""

    delta1, delta0 = map(lambda item: int(item) % P, target_residue)
    left0, left1 = map(lambda item: int(item) % P, transported_residue)
    if not delta0:
        raise HSKGeneratorShamashError(
            "the target residue is not a Frobenius generator"
        )
    value = left1 * pow(delta0, -1, P) % P
    tangent = (left0 - delta1 * value) * pow(delta0, -1, P) % P
    pair = value, tangent
    replay = (
        np.asarray(target_residue, dtype=np.int64).reshape(1, 2)
        @ dual_matrix(pair)
    ) % P
    if not np.array_equal(
        replay,
        np.asarray(transported_residue, dtype=np.int64).reshape(1, 2)
        % P,
    ):
        raise HSKGeneratorShamashError(
            "the local Frobenius multiplier failed replay"
        )
    return pair


def evaluated_lambda_row(h0_rows, residue_row, multiplication, inverses):
    """Evaluate ``-tau H5...H0`` directly in the local coefficient algebra.

    The formal source equations in the adapted Laurent frame can be large.
    Expanding their sixfold change-of-rings cochain is unnecessary: the
    counit system is imposed only after passage to the retained dual-number
    factor.  Evaluation is a ring homomorphism, so the same cochain is obtained
    by first evaluating the eight H0 blocks and then performing the exterior
    pullbacks as 1-by-2 by 2-by-2 matrix products.
    """

    full_wedge = tuple(range(VARIABLE_COUNT))
    cochain = {
        full_wedge: np.asarray(residue_row, dtype=np.int64).reshape(1, 2) % P
    }
    for step, h0 in enumerate(h0_rows):
        input_degree = 7 - step
        degree_sign = -1 if input_degree % 2 else 1
        evaluated_h = {
            variable: evaluate_block_multi(
                h0.get(((variable,), ()), {}), multiplication, inverses
            )
            for variable in range(VARIABLE_COUNT)
        }
        next_cochain = {}
        for wedge in itertools.combinations(range(VARIABLE_COUNT), input_degree):
            wedge_set = set(wedge)
            value = np.zeros((1, 2), dtype=np.int64)
            for variable in range(VARIABLE_COUNT):
                if variable in wedge_set:
                    continue
                output_wedge = tuple(sorted((*wedge, variable)))
                outer = cochain.get(output_wedge)
                if outer is None:
                    continue
                wedge_sign = -1 if sum(
                    1 for item in wedge if item < variable
                ) % 2 else 1
                value = (
                    value
                    + degree_sign
                    * wedge_sign
                    * (outer @ evaluated_h[variable])
                ) % P
            if np.any(value):
                next_cochain[tuple(wedge)] = value
        cochain = next_cochain

    # Reversing six odd exterior operators contributes (-1)^15=-1.
    result = {}
    for wedge, matrix in cochain.items():
        matrix = -np.asarray(matrix, dtype=np.int64) % P
        if np.any(matrix):
            result[tuple(wedge)] = {ZERO: matrix}
    return result


def verify_evaluated_lambda_against_canonical(
    equations,
    values,
    tangents,
    residue_row,
    evaluated_row,
    multiplication,
):
    """Cross-check the direct recurrence against the existing canonical code."""

    formal = lambda_transport._lambda_row(
        equations, values, tangents, residue_row
    )
    wedges = set(formal) | set(evaluated_row)
    for wedge in wedges:
        left = evaluate_row_multi(
            formal.get(wedge, {}), multiplication, {}
        )
        right = evaluate_row_multi(
            evaluated_row.get(wedge, {}), multiplication, {}
        )
        if np.any((left - right) % P):
            raise HSKGeneratorShamashError(
                "the evaluated lambda recurrence disagrees with the canonical code"
            )


def rref_solve_or_obstruction(system, right_hand_side):
    """Solve Ax=b or return an exact left-nullspace obstruction witness."""

    system = np.asarray(system, dtype=np.int64) % P
    right_hand_side = np.asarray(right_hand_side, dtype=np.int64).reshape(-1) % P
    row_count, column_count = system.shape
    augmented = np.column_stack(
        (
            system,
            right_hand_side,
            np.eye(row_count, dtype=np.int64),
        )
    ) % P
    row = 0
    pivots = []
    for column in range(column_count):
        candidates = np.flatnonzero(augmented[row:, column] % P)
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        if pivot != row:
            augmented[[row, pivot]] = augmented[[pivot, row]]
        augmented[row] = (
            augmented[row] * pow(int(augmented[row, column]), -1, P)
        ) % P
        for other in range(row_count):
            if other == row:
                continue
            scale = int(augmented[other, column]) % P
            if scale:
                augmented[other] = (
                    augmented[other] - scale * augmented[row]
                ) % P
        pivots.append(column)
        row += 1
        if row == row_count:
            break

    for remaining in range(row, row_count):
        if not np.any(augmented[remaining, :column_count] % P) and int(
            augmented[remaining, column_count]
        ) % P:
            witness = augmented[remaining, column_count + 1 :] % P
            obstruction = int(witness @ right_hand_side % P)
            if np.any(witness @ system % P) or not obstruction:
                raise HSKGeneratorShamashError(
                    "an inconsistent row did not give a left-nullspace witness"
                )
            return None, {
                "status": "OBSTRUCTED",
                "rank": len(pivots),
                "matrix_shape": [row_count, column_count],
                "left_null_witness": list(map(int, witness)),
                "left_null_times_rhs_mod_31": obstruction,
                "system_sha256": array_sha256(system),
                "rhs_sha256": array_sha256(right_hand_side),
            }

    solution = np.zeros(column_count, dtype=np.int64)
    for pivot_row, pivot_column in enumerate(pivots):
        solution[pivot_column] = augmented[pivot_row, column_count]
    if np.any((system @ solution - right_hand_side) % P):
        raise HSKGeneratorShamashError("the RREF solution failed substitution")
    return solution, {
        "status": "SOLVED",
        "rank": len(pivots),
        "nullity": column_count - len(pivots),
        "matrix_shape": [row_count, column_count],
        "pivot_columns": pivots,
        "free_columns": [
            column for column in range(column_count) if column not in pivots
        ],
        "canonical_free_variables_zero": True,
        "system_sha256": array_sha256(system),
        "rhs_sha256": array_sha256(right_hand_side),
        "solution_sha256": array_sha256(solution),
    }


def solve_single_counit_row(
    target_d2,
    target_h0,
    delta_k2,
    delta_equations,
    multiplication,
    inverses,
):
    """Use all 34 A-columns to solve one A-linear row on target K1."""

    pairs = tuple(itertools.combinations(range(VARIABLE_COUNT), 2))
    columns = [
        (f"K2_{index}", wedge, delta_k2.get(wedge, {}), target_d2)
        for index, wedge in enumerate(pairs)
    ]
    columns.extend(
        (
            f"H0_{equation}",
            (),
            delta_equations[equation].get((), {}),
            target_h0[equation],
        )
        for equation in range(EQUATION_COUNT)
    )
    system = np.zeros((2 * len(columns), 2 * VARIABLE_COUNT), dtype=np.int64)
    rhs = np.zeros(2 * len(columns), dtype=np.int64)
    labels = []
    for column, (label, input_wedge, delta_polynomial, blocks) in enumerate(columns):
        labels.append(label)
        rhs[2 * column : 2 * column + 2] = evaluate_row_multi(
            delta_polynomial, multiplication, inverses
        )[0]
        for variable in range(VARIABLE_COUNT):
            polynomial = blocks.get(((variable,), input_wedge), {})
            if not polynomial:
                continue
            block = evaluate_block_multi(polynomial, multiplication, inverses)
            for source_coordinate in range(2):
                for target_coordinate in range(2):
                    system[
                        2 * column + target_coordinate,
                        2 * variable + source_coordinate,
                    ] = block[source_coordinate, target_coordinate] % P
    solution, ledger = rref_solve_or_obstruction(system, rhs)
    ledger["column_labels"] = labels
    ledger["one_row_constrains_all_six_equation_sectors"] = True
    ledger["system_matrix"] = system
    ledger["rhs_vector"] = rhs
    if solution is None:
        return None, ledger
    ledger["canonical_solution"] = solution
    a_row = {}
    for variable in range(VARIABLE_COUNT):
        row = solution[2 * variable : 2 * variable + 2].reshape(1, 2)
        if np.any(row):
            a_row[(variable,)] = {ZERO: row}
    return a_row, ledger


def verify_counit_residual(
    a_row,
    target_d2,
    target_h0,
    delta_k2,
    delta_equations,
    multiplication,
    inverses,
):
    residual_k2 = lambda_transport._compose_row_with_blocks(a_row, target_d2)
    lambda_transport._add_row_map(residual_k2, delta_k2, -1)
    residual_equations = []
    polynomially_nonzero = []
    for index, wedge in enumerate(itertools.combinations(range(VARIABLE_COUNT), 2)):
        polynomial = residual_k2.get(wedge, {})
        if polynomial:
            polynomially_nonzero.append(f"K2_{index}")
        if np.any(evaluate_row_multi(polynomial, multiplication, inverses)):
            raise HSKGeneratorShamashError("a K2 counit residual is nonzero")
    for equation in range(EQUATION_COUNT):
        residual = lambda_transport._compose_row_with_blocks(
            a_row, target_h0[equation]
        )
        lambda_transport._add_row_map(residual, delta_equations[equation], -1)
        residual_equations.append(residual)
        polynomial = residual.get((), {})
        if polynomial:
            polynomially_nonzero.append(f"H0_{equation}")
        if np.any(evaluate_row_multi(polynomial, multiplication, inverses)):
            raise HSKGeneratorShamashError("an equation-sector counit residual is nonzero")
    return {
        "all_34_A_columns_zero_in_Hermite_quotient": True,
        "verified_F31_scalar_coordinates": 68,
        "polynomially_nonzero_before_Hermite_reduction": polynomially_nonzero,
    }


def verify_counit_class_invariance(
    lambda_source,
    edge_d3,
    target_d2,
    delta_k2,
    edge_multiplication,
    inverses,
):
    """Certify that the K2 residual class is independent of J2 homotopy.

    After restriction to the local coefficient algebra, every target D2 block
    is zero.  Hence any row of the form ``a D2`` dies.  Any two K2 comparison
    lifts differ by a D3 boundary because the localized graph equations remain
    a regular Koszul sequence.  The second check below certifies directly that
    lambda annihilates every such D3 boundary in the same quotient.
    """

    target_d2_nonzero = []
    for (output_wedge, input_wedge), polynomial in target_d2.items():
        value = evaluate_block_multi(
            polynomial, edge_multiplication, inverses
        )
        if np.any(value):
            target_d2_nonzero.append([list(output_wedge), list(input_wedge)])
    if target_d2_nonzero:
        raise HSKGeneratorShamashError(
            "target D2 did not vanish in the local coefficient quotient"
        )

    lambda_d3 = lambda_transport._compose_row_with_blocks(
        lambda_source, edge_d3
    )
    lambda_d3_nonzero = []
    for wedge in itertools.combinations(range(VARIABLE_COUNT), 3):
        value = evaluate_row_multi(
            lambda_d3.get(wedge, {}), edge_multiplication, inverses
        )
        if np.any(value):
            lambda_d3_nonzero.append(list(wedge))
    if lambda_d3_nonzero:
        raise HSKGeneratorShamashError(
            "lambda did not annihilate D3 in the local coefficient quotient"
        )

    residual_coordinates = []
    for pair_index, wedge in enumerate(
        itertools.combinations(range(VARIABLE_COUNT), 2)
    ):
        value = evaluate_row_multi(
            delta_k2.get(wedge, {}), edge_multiplication, inverses
        )[0]
        for local_coordinate, coefficient in enumerate(value):
            coefficient = int(coefficient) % P
            if coefficient:
                residual_coordinates.append(
                    {
                        "pair_index": pair_index,
                        "wedge": list(wedge),
                        "local_coordinate": local_coordinate,
                        "coefficient_mod_31": coefficient,
                    }
                )
    return {
        "target_D2_zero_in_local_coefficient_quotient": True,
        "target_D2_blocks_checked": len(target_d2),
        "lambda_source_D3_zero_in_local_coefficient_quotient": True,
        "source_D3_columns_checked": 56,
        "localized_graph_sequence_regular": True,
        "alternative_K2_lifts_differ_by_D3_boundaries": True,
        "K2_counit_class_independent_of_comparison_homotopy": True,
        "nonzero_obstruction_class_present": bool(residual_coordinates),
        "nonzero_K2_residual_coordinates": residual_coordinates,
    }


def equation_change_data(
    factor_data,
    record,
    source_equations,
    target_equations,
    target_substitutions,
):
    """Compute the one diagonal B law and verify all six equation identities."""

    b_entries = []
    for equation in range(EQUATION_COUNT):
        exponent = list(ZERO)
        scalar = 1
        if equation < 4:
            factor = equation
            scalar = int(record["curve_scalars"][factor]) % P
            inverse_variable = factor_data[factor]["inverse_variable"]
            if inverse_variable is not None:
                exponent[inverse_variable] = -3
        else:
            label = (2, 3)[equation - 4]
            scalar = int(record["theta_characters"][label]) % P
            for data in factor_data:
                if data["inverse_variable"] is not None:
                    exponent[data["inverse_variable"]] -= 1
        b_polynomial = {tuple(exponent): scalar}
        pulled_target = substitute_scalar_polynomial(
            target_equations[equation], target_substitutions
        )
        expected = scalar_multiply(source_equations[equation], b_polynomial)
        if pulled_target != expected:
            raise HSKGeneratorShamashError(
                f"equation {equation} did not obey the common diagonal B law"
            )
        b_entries.append(
            {
                "equation": equation,
                "exponent": tuple(exponent),
                "scalar": scalar,
                "polynomial": b_polynomial,
                "pulled_terms": len(pulled_target),
            }
        )
    return b_entries


def local_shamash_support(
    adapted_values,
    adapted_tangents,
    target_values,
    target_tangents,
    source_equations,
    target_equations,
    target_substitutions,
    factor_data,
    b_entries,
):
    """Construct J0 through J3 and the six common Gamma/Omega corrections."""

    edge_multiplication = shamash.local_multiplication(
        adapted_values, adapted_tangents
    )
    target_multiplication = shamash.local_multiplication(
        target_values, target_tangents
    )
    inverses = {
        data["inverse_variable"]: shamash.inverse_2x2(
            edge_multiplication[data["inverse_variable"]]
        )
        for data in factor_data
        if data["inverse_variable"] is not None
    }
    edge_differentials = {
        degree: shamash.graph_differential(degree, edge_multiplication)
        for degree in (1, 2, 3)
    }
    target_raw_differentials = {
        degree: shamash.graph_differential(degree, target_multiplication)
        for degree in (1, 2, 3)
    }
    target_differentials = {
        degree: substitute_block_map(
            target_raw_differentials[degree], target_substitutions
        )
        for degree in (1, 2, 3)
    }
    graph_maps = graph_compounds_general(target_multiplication, factor_data)
    for degree in (1, 2, 3):
        shamash.assert_block_maps_equal(
            shamash.compose_block_maps(edge_differentials[degree], graph_maps[degree]),
            shamash.compose_block_maps(
                graph_maps[degree - 1], target_differentials[degree]
            ),
            f"general K graph J{degree}",
        )

    gamma = []
    omega = []
    edge_h0_rows = []
    edge_h1_rows = []
    target_h0_rows = []
    target_h1_rows = []
    target_h0_raw_rows = []
    target_h1_raw_rows = []
    delta_nonzero = []
    xi_nonzero = []
    for equation in range(EQUATION_COUNT):
        if np.any(
            shamash.evaluate_matrix_polynomial(
                source_equations[equation], edge_multiplication
            )
        ):
            raise HSKGeneratorShamashError("an adapted source equation is nonzero")
        if np.any(
            shamash.evaluate_matrix_polynomial(
                target_equations[equation], target_multiplication
            )
        ):
            raise HSKGeneratorShamashError("a target equation is nonzero")
        edge_h0 = shamash.divided_difference_matrix(
            source_equations[equation], edge_multiplication
        )
        edge_h1 = shamash.h1_from_h0(edge_h0)
        target_h0_raw = shamash.divided_difference_matrix(
            target_equations[equation], target_multiplication
        )
        target_h1_raw = shamash.h1_from_h0(target_h0_raw)
        target_h0 = substitute_block_map(target_h0_raw, target_substitutions)
        target_h1 = substitute_block_map(target_h1_raw, target_substitutions)
        shamash.verify_h_nullhomotopy(
            source_equations[equation],
            edge_h0,
            edge_h1,
            edge_differentials,
            "adapted edge",
        )
        shamash.verify_h_nullhomotopy(
            target_equations[equation],
            target_h0_raw,
            target_h1_raw,
            target_raw_differentials,
            "canonical target",
        )
        b = b_entries[equation]
        shifted_edge_h0 = shamash.shift_block_map(
            shamash.compose_block_maps(edge_h0, graph_maps[0]),
            b["exponent"],
            scale=b["scalar"],
        )
        delta = shamash.add_block_maps(
            shamash.compose_block_maps(graph_maps[1], target_h0),
            shifted_edge_h0,
            right_scale=-1,
        )
        if delta:
            delta_nonzero.append(equation)
        if shamash.compose_block_maps(edge_differentials[1], delta):
            raise HSKGeneratorShamashError("a Delta is not a d1 cycle")
        gamma_j = shamash.shift_block_map(
            tensor_ordered_homotopy_multi(
                delta, 1, edge_multiplication, inverses
            ),
            ZERO,
            scale=-1,
        )
        shamash.assert_block_maps_equal(
            shamash.compose_block_maps(edge_differentials[2], gamma_j),
            delta,
            "d2 Gamma=Delta",
        )
        # This is the upper-right equation column of the complete J2 identity.
        left_equation = shamash.add_block_maps(
            shamash.compose_block_maps(edge_differentials[2], gamma_j),
            shamash.shift_block_map(
                edge_h0, b["exponent"], scale=b["scalar"]
            ),
        )
        right_equation = shamash.compose_block_maps(graph_maps[1], target_h0)
        shamash.assert_block_maps_equal(
            left_equation,
            right_equation,
            "D2_edge J2=J1 D2_target equation column",
        )

        xi = shamash.compose_block_maps(graph_maps[2], target_h1)
        xi = shamash.add_block_maps(
            xi, shamash.compose_block_maps(gamma_j, target_differentials[1])
        )
        xi = shamash.add_block_maps(
            xi,
            shamash.shift_block_map(
                shamash.compose_block_maps(edge_h1, graph_maps[1]),
                b["exponent"],
                scale=b["scalar"],
            ),
            right_scale=-1,
        )
        if xi:
            xi_nonzero.append(equation)
        if shamash.compose_block_maps(edge_differentials[2], xi):
            raise HSKGeneratorShamashError("an Xi is not a d2 cycle")
        omega_j = tensor_ordered_homotopy_multi(
            xi, 2, edge_multiplication, inverses
        )
        shamash.assert_block_maps_equal(
            shamash.compose_block_maps(edge_differentials[3], omega_j),
            xi,
            "d3 Omega=Xi",
        )
        gamma.append(gamma_j)
        omega.append(omega_j)
        edge_h0_rows.append(edge_h0)
        edge_h1_rows.append(edge_h1)
        target_h0_rows.append(target_h0)
        target_h1_rows.append(target_h1)
        target_h0_raw_rows.append(target_h0_raw)
        target_h1_raw_rows.append(target_h1_raw)

    # The graph K2 column is the seventh and only other component of J2.
    shamash.assert_block_maps_equal(
        shamash.compose_block_maps(edge_differentials[2], graph_maps[2]),
        shamash.compose_block_maps(graph_maps[1], target_differentials[2]),
        "D2_edge J2=J1 D2_target graph column",
    )
    return {
        "edge_multiplication": edge_multiplication,
        "target_multiplication": target_multiplication,
        "inverses": inverses,
        "edge_differentials": edge_differentials,
        "target_differentials": target_differentials,
        "graph_maps": graph_maps,
        "gamma": gamma,
        "omega": omega,
        "edge_h0": edge_h0_rows,
        "edge_h1": edge_h1_rows,
        "target_h0": target_h0_rows,
        "target_h1": target_h1_rows,
        "target_h0_raw": target_h0_raw_rows,
        "target_h1_raw": target_h1_raw_rows,
        "delta_nonzero": delta_nonzero,
        "xi_nonzero": xi_nonzero,
        "full_J2_component_identities": 7,
        "full_J3_bound_inputs": 13,
    }


def serialize_block_maps(support_results, key, degree_out, degree_in, equation=None):
    records = []
    for support_position, result in enumerate(support_results):
        value = result[key] if equation is None else result[key][equation]
        prefix = (support_position,) if equation is None else (support_position, equation)
        records.extend(
            shamash._serialize_block_map(
                value, degree_out, degree_in, prefix
            )
        )
    width = len(prefix) + 13
    if records:
        return np.asarray(records, dtype=np.int16).reshape((-1, width))
    return np.zeros((0, width), dtype=np.int16)


def serialize_row(support_rows, key):
    wedge_indices = {
        0: {(): 0},
        1: {(variable,): variable for variable in range(VARIABLE_COUNT)},
        2: {
            wedge: index
            for index, wedge in enumerate(
                itertools.combinations(range(VARIABLE_COUNT), 2)
            )
        },
    }
    records = []
    for support_position, record in enumerate(support_rows):
        row = record.get(key) or {}
        for wedge, polynomial in sorted(row.items()):
            for exponent, matrix in sorted(polynomial.items()):
                for local_column in range(matrix.shape[1]):
                    coefficient = int(matrix[0, local_column]) % P
                    if coefficient:
                        records.append(
                            [
                                support_position,
                                len(wedge),
                                wedge_indices[len(wedge)][tuple(wedge)],
                                *map(int, exponent),
                                local_column,
                                coefficient,
                            ]
                        )
    if records:
        return np.asarray(records, dtype=np.int16).reshape((-1, 13))
    return np.zeros((0, 13), dtype=np.int16)


def serialize_graph_degree(support_results, degree):
    records = []
    for support_position, result in enumerate(support_results):
        records.extend(
            shamash._serialize_block_map(
                result["graph_maps"][degree],
                degree,
                degree,
                (support_position,),
            )
        )
    if records:
        return np.asarray(records, dtype=np.int16).reshape((-1, 14))
    return np.zeros((0, 14), dtype=np.int16)


def serialize_equation_block_maps(
    support_results, key, degree_out, degree_in
):
    records = []
    for support_position, result in enumerate(support_results):
        for equation, value in enumerate(result[key]):
            records.extend(
                shamash._serialize_block_map(
                    value,
                    degree_out,
                    degree_in,
                    (support_position, equation),
                )
            )
    if records:
        return np.asarray(records, dtype=np.int16).reshape((-1, 15))
    return np.zeros((0, 15), dtype=np.int16)


def build():
    graph_input = k_j1.build_data()
    context = iw_solver.load_problem_context()
    canonical_equations, equation_labels = shamash.surface_equations(
        context, SURFACE, CHART
    )
    descent, _extension, _stable, _elements, _refined, _edge, theta, labels = (
        k_chain.load_inputs()
    )
    records = k_chain.group_records(descent, theta)
    record = next(
        item for item in records if tuple(item["exponents"]) == GENERATOR_EXPONENTS
    )
    if labels != [2, 3]:
        raise HSKGeneratorShamashError("the Z_2 theta labels changed")
    for left, right in zip(record["matrices"], graph_input["record"]["factor_matrices"]):
        if np.any((np.asarray(left) - np.asarray(right)) % P):
            raise HSKGeneratorShamashError("the two K matrix ledgers disagree")

    bits = tuple(map(int, context["chart_bits"][CHART]))
    factor_data = []
    frame_rows = []
    original_substitutions = [None] * VARIABLE_COUNT
    target_substitutions = [None] * VARIABLE_COUNT
    for factor, matrix in enumerate(record["matrices"]):
        frame, numerator_u, numerator_v = k_j1.affine_forms(matrix, bits[factor])
        data = adapted_factor(frame, numerator_u, numerator_v, factor)
        factor_data.append(data)
        for local in range(2):
            original_substitutions[2 * factor + local] = data[
                "original_substitutions"
            ][local]
            target_substitutions[2 * factor + local] = data[
                "target_substitutions"
            ][local]
        frame_rows.append(
            {
                "factor": factor,
                "original_frame_L": json.dumps(list(data["frame_original"])),
                "adapted_transform_rows": json.dumps(
                    np.asarray(data["transform"], dtype=int).tolist()
                ),
                "adapted_inverse_rows": json.dumps(
                    np.asarray(data["inverse"], dtype=int).tolist()
                ),
                "target_u_in_1_a_b": json.dumps(
                    list(data["numerator_u_adapted"])
                ),
                "target_v_in_1_a_b": json.dumps(
                    list(data["numerator_v_adapted"])
                ),
                "inverse_variable": (
                    "" if data["inverse_variable"] is None else data["inverse_variable"]
                ),
            }
        )

    source_equations = [
        substitute_scalar_polynomial(polynomial, original_substitutions)
        for polynomial in canonical_equations
    ]
    adapted_orientation_scalar = 1
    for data in factor_data:
        adapted_orientation_scalar = (
            adapted_orientation_scalar * determinant_mod(data["transform"])
        ) % P
    if not adapted_orientation_scalar:
        raise HSKGeneratorShamashError(
            "the adapted-coordinate top orientation is not a unit"
        )
    b_entries = equation_change_data(
        factor_data,
        record,
        source_equations,
        canonical_equations,
        target_substitutions,
    )
    equation_rows = [
        {
            "equation": entry["equation"],
            "label": equation_labels[entry["equation"]],
            "B_scalar_mod_31": entry["scalar"],
            "B_Laurent_exponents_adapted_coordinates": json.dumps(
                list(entry["exponent"])
            ),
            "pulled_target_terms": entry["pulled_terms"],
            "identity": "pull(target_equation)=B_j*source_equation",
        }
        for entry in b_entries
    ]

    orientation = all_charts.orientation_map(SURFACE)
    support_results = []
    counit_rows = []
    support_summary = []
    residue_transport_pairs = []
    orientation_coboundary_pairs = []
    raw_counit_weight_pairs = []
    target_lambda_normalization_pairs = []
    weight_exponent = list(ZERO)
    for data in factor_data:
        if data["inverse_variable"] is not None:
            weight_exponent[data["inverse_variable"]] += 4
    weight_exponent = tuple(weight_exponent)
    weight_scalar = pow(int(record["H_normalization"]), 4, P)
    b_determinant_exponent = tuple(
        sum(int(entry["exponent"][variable]) for entry in b_entries)
        for variable in range(VARIABLE_COUNT)
    )
    b_determinant_scalar = 1
    for entry in b_entries:
        b_determinant_scalar = (
            b_determinant_scalar * int(entry["scalar"])
        ) % P
    if not b_determinant_scalar or not weight_scalar:
        raise HSKGeneratorShamashError(
            "a determinant or normalized O(4H) frame is not a unit"
        )

    for support_position, point_row in enumerate(graph_input["point_rows"]):
        source_support = tuple(json.loads(point_row["source_support"]))
        target_support = tuple(json.loads(point_row["target_support"]))
        source_values = tuple(json.loads(point_row["source_affine_values"]))
        source_tangents = tuple(json.loads(point_row["source_affine_tangent"]))
        target_values = tuple(json.loads(point_row["target_affine_values"]))
        target_tangents = tuple(json.loads(point_row["target_affine_tangent"]))
        adapted_values, adapted_tangents = adapt_jet(
            source_values, source_tangents, factor_data
        )
        comparison = local_shamash_support(
            adapted_values,
            adapted_tangents,
            target_values,
            target_tangents,
            source_equations,
            canonical_equations,
            target_substitutions,
            factor_data,
            b_entries,
        )
        support_results.append(comparison)

        if source_support not in orientation or target_support not in orientation:
            raise HSKGeneratorShamashError("a K support lacks a residue orientation")
        source_delta0, source_delta1 = orientation[source_support]
        target_delta0, target_delta1 = orientation[target_support]
        source_residue = [
            adapted_orientation_scalar * source_delta1 % P,
            adapted_orientation_scalar * source_delta0 % P,
        ]
        target_residue = [target_delta1, target_delta0]
        top_wedge = tuple(range(VARIABLE_COUNT))
        top_graph = comparison["graph_maps"][8][(top_wedge, top_wedge)]
        top_graph_evaluated = evaluate_block_multi(
            top_graph,
            comparison["edge_multiplication"],
            comparison["inverses"],
        )
        top_graph_pair = (
            int(top_graph_evaluated[0, 0]) % P,
            int(top_graph_evaluated[1, 0]) % P,
        )
        if not np.array_equal(top_graph_evaluated, dual_matrix(top_graph_pair)):
            raise HSKGeneratorShamashError(
                "the top graph is not one local multiplication operator"
            )
        transported_residue = (
            np.asarray(source_residue, dtype=np.int64).reshape(1, 2)
            @ top_graph_evaluated
        )[0] % P
        residue_transport_pair = frobenius_multiplier(
            target_residue, transported_residue
        )
        source_orientation_pair = (
            int(source_delta0) % P,
            int(source_delta1) % P,
        )
        target_orientation_pair = (
            int(target_delta0) % P,
            int(target_delta1) % P,
        )
        orientation_coboundary_pair = dual_divide(
            source_orientation_pair, target_orientation_pair
        )
        formula_pair = dual_multiply(
            dual_multiply(
                (adapted_orientation_scalar, 0),
                orientation_coboundary_pair,
            ),
            top_graph_pair,
        )
        if formula_pair != residue_transport_pair:
            raise HSKGeneratorShamashError(
                "the orientation-coboundary residue law changed"
            )
        b_determinant_pair = (1, 0)
        for entry in b_entries:
            b_determinant_pair = dual_multiply(
                b_determinant_pair,
                scalar_laurent_pair(
                    entry["scalar"],
                    entry["exponent"],
                    adapted_values,
                    adapted_tangents,
                ),
            )
        local_weight_pair = scalar_laurent_pair(
            weight_scalar,
            weight_exponent,
            adapted_values,
            adapted_tangents,
        )
        raw_counit_weight_pair = dual_divide(
            residue_transport_pair, b_determinant_pair
        )
        target_lambda_normalization_pair = dual_divide(
            raw_counit_weight_pair, local_weight_pair
        )
        dual_inverse(residue_transport_pair)
        dual_inverse(raw_counit_weight_pair)
        dual_inverse(target_lambda_normalization_pair)
        if dual_multiply(
            dual_multiply(
                b_determinant_pair, local_weight_pair
            ),
            target_lambda_normalization_pair,
        ) != residue_transport_pair:
            raise HSKGeneratorShamashError(
                "the supportwise counit normalization factorization changed"
            )
        residue_transport_pairs.append(residue_transport_pair)
        orientation_coboundary_pairs.append(
            orientation_coboundary_pair
        )
        raw_counit_weight_pairs.append(raw_counit_weight_pair)
        target_lambda_normalization_pairs.append(
            target_lambda_normalization_pair
        )
        lambda_source = evaluated_lambda_row(
            comparison["edge_h0"],
            source_residue,
            comparison["edge_multiplication"],
            comparison["inverses"],
        )
        lambda_target_raw = evaluated_lambda_row(
            comparison["target_h0_raw"],
            target_residue,
            comparison["target_multiplication"],
            {},
        )
        verify_evaluated_lambda_against_canonical(
            canonical_equations,
            target_values,
            target_tangents,
            target_residue,
            lambda_target_raw,
            comparison["target_multiplication"],
        )
        normalized_target_residue = (
            np.asarray(target_residue, dtype=np.int64).reshape(1, 2)
            @ dual_matrix(target_lambda_normalization_pair)
        )[0] % P
        lambda_target_normalized = evaluated_lambda_row(
            comparison["target_h0_raw"],
            normalized_target_residue,
            comparison["target_multiplication"],
            {},
        )
        verify_evaluated_lambda_against_canonical(
            canonical_equations,
            target_values,
            target_tangents,
            normalized_target_residue,
            lambda_target_normalized,
            comparison["target_multiplication"],
        )
        source_composed = lambda_transport._compose_row_with_blocks(
            lambda_source, comparison["graph_maps"][2]
        )
        target_weighted = lambda_transport._shift_row(
            lambda_target_normalized, weight_exponent, scale=weight_scalar
        )
        delta_k2 = dict(source_composed)
        lambda_transport._add_row_map(delta_k2, target_weighted, -1)
        delta_equations = [
            lambda_transport._compose_row_with_blocks(lambda_source, gamma)
            for gamma in comparison["gamma"]
        ]
        invariance_ledger = verify_counit_class_invariance(
            lambda_source,
            comparison["edge_differentials"][3],
            comparison["target_differentials"][2],
            delta_k2,
            comparison["edge_multiplication"],
            comparison["inverses"],
        )
        a_row, solve_ledger = solve_single_counit_row(
            comparison["target_differentials"][2],
            comparison["target_h0"],
            delta_k2,
            delta_equations,
            comparison["edge_multiplication"],
            comparison["inverses"],
        )
        residual_ledger = None
        if a_row is not None:
            residual_ledger = verify_counit_residual(
                a_row,
                comparison["target_differentials"][2],
                comparison["target_h0"],
                delta_k2,
                delta_equations,
                comparison["edge_multiplication"],
                comparison["inverses"],
            )
        counit_rows.append(
            {
                "a": a_row,
                "delta_k2": delta_k2,
                "delta_equations": delta_equations,
                "solve": solve_ledger,
                "residual": residual_ledger,
                "invariance": invariance_ledger,
            }
        )
        support_summary.append(
            {
                "support_position": support_position,
                "source_support": list(source_support),
                "target_support": list(target_support),
                "residue_transport_pair": list(residue_transport_pair),
                "orientation_coboundary_pair": list(
                    orientation_coboundary_pair
                ),
                "raw_counit_weight_pair": list(
                    raw_counit_weight_pair
                ),
                "target_lambda_normalization_pair": list(
                    target_lambda_normalization_pair
                ),
                "Gamma_terms": sum(
                    shamash.block_map_term_count(value)
                    for value in comparison["gamma"]
                ),
                "Omega_terms": sum(
                    shamash.block_map_term_count(value)
                    for value in comparison["omega"]
                ),
                "Delta_nonzero_equations": comparison["delta_nonzero"],
                "Xi_nonzero_equations": comparison["xi_nonzero"],
                "counit_status": solve_ledger["status"],
                "counit_rank": solve_ledger["rank"],
                "nonzero_K2_residual_coordinates": invariance_ledger[
                    "nonzero_K2_residual_coordinates"
                ],
                "K2_counit_class_homotopy_invariant": invariance_ledger[
                    "K2_counit_class_independent_of_comparison_homotopy"
                ],
                "counit_a_terms": (
                    0 if a_row is None else lambda_transport._row_term_count(a_row)
                ),
                "counit_a_sha256": None if a_row is None else lambda_transport._row_hash(a_row),
                "counit_obstruction": (
                    None
                    if solve_ledger["status"] == "SOLVED"
                    else {
                        "left_null_witness": solve_ledger["left_null_witness"],
                        "left_null_times_rhs_mod_31": solve_ledger[
                            "left_null_times_rhs_mod_31"
                        ],
                    }
                ),
                "all_34_counit_columns_zero": bool(
                    residual_ledger
                    and residual_ledger[
                        "all_34_A_columns_zero_in_Hermite_quotient"
                    ]
                ),
            }
        )

    solved = sum(row["solve"]["status"] == "SOLVED" for row in counit_rows)
    normalization_is_constant = (
        len(set(target_lambda_normalization_pairs)) == 1
        and target_lambda_normalization_pairs[0][1] == 0
    )
    return {
        "factor_data": factor_data,
        "frame_rows": frame_rows,
        "source_equations": source_equations,
        "target_substitutions": target_substitutions,
        "b_entries": b_entries,
        "equation_rows": equation_rows,
        "support_results": support_results,
        "counit_rows": counit_rows,
        "support_summary": support_summary,
        "record": record,
        "weight_exponent": weight_exponent,
        "weight_scalar": weight_scalar,
        "residue_transport_pairs": residue_transport_pairs,
        "orientation_coboundary_pairs": orientation_coboundary_pairs,
        "raw_counit_weight_pairs": raw_counit_weight_pairs,
        "target_lambda_normalization_pairs": (
            target_lambda_normalization_pairs
        ),
        "target_lambda_normalization_scalar": (
            target_lambda_normalization_pairs[0][0]
            if normalization_is_constant
            else None
        ),
        "normalization": {
            "B_determinant_scalar": b_determinant_scalar,
            "B_determinant_exponent": b_determinant_exponent,
            "normalized_weight_scalar": weight_scalar,
            "normalized_weight_exponent": weight_exponent,
            "target_lambda_normalization_is_constant": (
                normalization_is_constant
            ),
            "target_lambda_normalization_constant_pair": (
                list(target_lambda_normalization_pairs[0])
                if normalization_is_constant
                else None
            ),
            "target_lambda_normalization_pairs_sha256": pair_ledger_sha256(
                target_lambda_normalization_pairs
            ),
            "residue_transport_pairs_sha256": pair_ledger_sha256(
                residue_transport_pairs
            ),
            "raw_counit_weight_pairs_sha256": pair_ledger_sha256(
                raw_counit_weight_pairs
            ),
            "orientation_coboundary_pairs_sha256": pair_ledger_sha256(
                orientation_coboundary_pairs
            ),
            "pair_ledger_hash_encoding": "signed-little-endian-int64-C-order",
            "orientation_coboundary_distinct_pairs": len(
                set(orientation_coboundary_pairs)
            ),
            "top_J8_residue_transport_checks": len(support_results),
            "adapted_source_orientation_scalar": (
                adapted_orientation_scalar
            ),
        },
        "counts": {
            "supports": len(support_results),
            "equations": EQUATION_COUNT,
            "J2_component_identities": 7 * len(support_results),
            "D3_compatibility_inputs": 13 * len(support_results),
            "Gamma_polynomial_terms": sum(
                shamash.block_map_term_count(value)
                for result in support_results
                for value in result["gamma"]
            ),
            "Omega_polynomial_terms": sum(
                shamash.block_map_term_count(value)
                for result in support_results
                for value in result["omega"]
            ),
            "counit_systems_solved": solved,
            "counit_systems_obstructed": len(counit_rows) - solved,
        },
    }


def write_artifacts(output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    result = build()

    maps_path = output_directory / "generator_shamash_maps.npz"
    systems_path = output_directory / "generator_counit_systems.npz"
    frames_path = output_directory / "generator_adapted_frames.csv"
    equations_path = output_directory / "generator_equation_change.csv"
    supports_path = output_directory / "generator_support_summary.csv"
    certificate_path = output_directory / "generator_shamash_certificate.json"

    np.savez_compressed(
        maps_path,
        field=np.asarray([P], dtype=np.uint8),
        generator_exponents=np.asarray(GENERATOR_EXPONENTS, dtype=np.uint8),
        graph_J1_records=serialize_graph_degree(result["support_results"], 1),
        graph_J2_records=serialize_graph_degree(result["support_results"], 2),
        graph_J3_records=serialize_graph_degree(result["support_results"], 3),
        Gamma_records=serialize_equation_block_maps(
            result["support_results"], "gamma", 2, 0
        ),
        Omega_records=serialize_equation_block_maps(
            result["support_results"], "omega", 3, 1
        ),
        counit_a_records=serialize_row(result["counit_rows"], "a"),
        counit_delta_K2_records=serialize_row(
            result["counit_rows"], "delta_k2"
        ),
        normalized_weight_scalar=np.asarray(
            [result["weight_scalar"]], dtype=np.uint8
        ),
        normalized_weight_exponent=np.asarray(
            result["weight_exponent"], dtype=np.int8
        ),
        target_lambda_normalization_pairs=np.asarray(
            result["target_lambda_normalization_pairs"], dtype=np.uint8
        ),
        residue_transport_pairs=np.asarray(
            result["residue_transport_pairs"], dtype=np.uint8
        ),
        raw_counit_weight_pairs=np.asarray(
            result["raw_counit_weight_pairs"], dtype=np.uint8
        ),
        orientation_coboundary_pairs=np.asarray(
            result["orientation_coboundary_pairs"], dtype=np.uint8
        ),
        adapted_source_orientation_scalar=np.asarray(
            [result["normalization"]["adapted_source_orientation_scalar"]],
            dtype=np.uint8,
        ),
    )
    systems = np.stack(
        [row["solve"]["system_matrix"] for row in result["counit_rows"]]
    )
    right_hand_sides = np.stack(
        [row["solve"]["rhs_vector"] for row in result["counit_rows"]]
    )
    solutions = np.stack(
        [row["solve"]["canonical_solution"] for row in result["counit_rows"]]
    )
    np.savez_compressed(
        systems_path,
        field=np.asarray([P], dtype=np.uint8),
        systems=np.asarray(systems, dtype=np.uint8),
        right_hand_sides=np.asarray(right_hand_sides, dtype=np.uint8),
        canonical_solutions=np.asarray(solutions, dtype=np.uint8),
        ranks=np.asarray(
            [row["solve"]["rank"] for row in result["counit_rows"]],
            dtype=np.uint8,
        ),
        nullities=np.asarray(
            [row["solve"]["nullity"] for row in result["counit_rows"]],
            dtype=np.uint8,
        ),
    )
    write_csv(frames_path, result["frame_rows"])
    write_csv(equations_path, result["equation_rows"])
    support_rows = []
    for row in result["support_summary"]:
        output = {}
        for key, value in row.items():
            output[key] = (
                json.dumps(value, sort_keys=True)
                if isinstance(value, (list, dict))
                else value
            )
        support_rows.append(output)
    write_csv(supports_path, support_rows)

    all_counit_residuals_zero = all(
        bool(row["residual"])
        and row["residual"][
            "all_34_A_columns_zero_in_Hermite_quotient"
        ]
        for row in result["counit_rows"]
    )
    all_homotopy_checks = all(
        row["invariance"][
            "K2_counit_class_independent_of_comparison_homotopy"
        ]
        for row in result["counit_rows"]
    )
    certificate = {
        "schema": "section12-HS-rank2-K-generator-Shamash-v1",
        "field": "F_31",
        "status": (
            "supportwise-localized-full-J2-exact;"
            "D3-compatibility-inputs-exact;counit-a-solved;"
            "tree-contracted-global-matrix-open"
        ),
        "surface": "Z_2",
        "chart": CHART,
        "K_generator_exponents": list(GENERATOR_EXPONENTS),
        "counts": result["counts"],
        "normalization": {
            key: list(value) if isinstance(value, tuple) else value
            for key, value in result["normalization"].items()
        },
        "exact_checks": {
            "one_common_diagonal_B_law_for_all_six_equations": True,
            "Gamma_not_six_fitted_repairs": True,
            "full_supportwise_J2_component_identities": result["counts"][
                "J2_component_identities"
            ],
            "D3_compatibility_inputs": result["counts"][
                "D3_compatibility_inputs"
            ],
            "top_J8_residue_transport_checks": result["normalization"][
                "top_J8_residue_transport_checks"
            ],
            "one_orientation_coboundary_law_for_all_supports": True,
            "raw_counit_weight_r_over_det_B": True,
            "target_normalization_q_over_O4H": True,
            "adapted_coordinate_orientation_determinant_exact": True,
            "canonical_target_lambda_recurrence_cross_checks": result["counts"][
                "supports"
            ],
            "single_counit_row_systems_solved": result["counts"][
                "counit_systems_solved"
            ],
            "all_34_counit_columns_zero_in_local_quotient": (
                all_counit_residuals_zero
            ),
            "K2_counit_class_homotopy_invariance_checks": (
                result["counts"]["supports"] if all_homotopy_checks else 0
            ),
            "nonzero_counit_obstruction_witnesses": result["counts"][
                "counit_systems_obstructed"
            ],
        },
        "counit_system": {
            "shape_per_support": [68, 16],
            "rank_per_support": sorted(
                {row["solve"]["rank"] for row in result["counit_rows"]}
            ),
            "nullity_per_support": sorted(
                {row["solve"]["nullity"] for row in result["counit_rows"]}
            ),
            "one_row_used_for_all_28_K2_and_six_equation_columns": True,
            "canonical_free_variables_set_to_zero": True,
            "canonical_solution_a_is_zero_on_all_supports": all(
                not row["a"] for row in result["counit_rows"]
            ),
        },
        "claim_boundary": {
            "proved": (
                f"on all {result['counts']['supports']} retained dual-number "
                f"overlap factors for generator {GENERATOR_EXPONENTS}, the "
                "exact Laurent graph comparison, one "
                "six-equation B law, six Gamma homotopies, full supportwise J2 "
                "identity, Omega/D3 compatibility inputs, and one common counit "
                "row a"
            ),
            "not_claimed": (
                "the rectangular tree-contracted J2 matrix, a "
                "square global K linearization, all K generators or coherence, "
                "the global cross-Hom/Yoneda operator, Question 11.4, or the "
                "Hodge conjecture"
            ),
            "next_exact_construction": (
                "contract the certified supportwise maps through the source and "
                "target border-basis trees and verify the rectangular global "
                "D2 J2 = J1 D2 identity without changing this local law"
            ),
        },
        "files": {
            "maps": maps_path.name,
            "counit_systems": systems_path.name,
            "adapted_frames": frames_path.name,
            "equation_change": equations_path.name,
            "support_summary": supports_path.name,
        },
        "inputs": {
            "K_J1_producer": {
                "path": str(Path(k_j1.__file__).resolve().relative_to(ROOT)),
                "sha256": sha256(Path(k_j1.__file__).resolve()),
            },
            "K_chain_producer": {
                "path": str(Path(k_chain.__file__).resolve().relative_to(ROOT)),
                "sha256": sha256(Path(k_chain.__file__).resolve()),
            },
            "ordinary_Shamash_producer": {
                "path": str(Path(shamash.__file__).resolve().relative_to(ROOT)),
                "sha256": sha256(Path(shamash.__file__).resolve()),
            },
        },
    }
    write_json_atomic(certificate_path, certificate)

    files = {}
    for path in (
        maps_path,
        systems_path,
        frames_path,
        equations_path,
        supports_path,
        certificate_path,
    ):
        files[path.name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    manifest = {
        "schema": "section12-HS-rank2-K-generator-Shamash-manifest-v1",
        "generator": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "files": files,
    }
    write_json_atomic(output_directory / "manifest.json", manifest)
    return certificate


def verify_artifacts(output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    certificate = json.loads(
        (output_directory / "generator_shamash_certificate.json").read_text(
            encoding="utf-8"
        )
    )
    manifest = json.loads(
        (output_directory / "manifest.json").read_text(encoding="utf-8")
    )
    if certificate.get("schema") != "section12-HS-rank2-K-generator-Shamash-v1":
        raise HSKGeneratorShamashError("the certificate schema changed")
    for name, record in manifest["files"].items():
        path = output_directory / name
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise HSKGeneratorShamashError("an emitted artifact hash changed")
    with np.load(output_directory / certificate["files"]["counit_systems"]) as data:
        systems = np.asarray(data["systems"], dtype=np.int64) % P
        rhs = np.asarray(data["right_hand_sides"], dtype=np.int64) % P
        solutions = np.asarray(data["canonical_solutions"], dtype=np.int64) % P
        expected_supports = int(certificate["counts"]["supports"])
        if systems.shape != (expected_supports, 68, 16):
            raise HSKGeneratorShamashError("the counit-system shape changed")
        if np.any(
            np.einsum("sij,sj->si", systems, solutions) % P - rhs
        ):
            raise HSKGeneratorShamashError(
                "a persisted counit solution failed substitution"
            )
    return certificate


def main():
    certificate = write_artifacts()
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
