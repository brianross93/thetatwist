"""Contract one certified K-generator comparison into rectangular tree bases.

This producer binds one standard K generator, rebuilds its finite overlap
algebra from the certified adapted jets, constructs its order-tree
Hartshorne--Serre presentation, and compares it with the independently
normalized length-98 target presentation.  The direction is

    target C0010  --->  source-C0010 restricted to the K overlap.

It is a restriction map, not an invertible chart-to-chart transition.

The rectangular maps are stored over the generator's product of dual-number
factors.  An entry is represented at each factor by
``value+tangent*eps``.  The producer checks the chain identities in every
factor; the dynamically sized invertible Hermite decomposition makes those
checks exhaustive on the finite overlap algebra.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_all_chart_presentations as all_charts,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_extension_from_surface_presentation as extension,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_rank2_k_generator_overlap_j1 as k_j1,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_rank2_k_generator_shamash as k_shamash,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_shamash_overlap_lifts as shamash,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_iw_hermite_border_basis as hermite,
)


P = 31
VARIABLE_COUNT = 8
EQUATION_COUNT = 6
SURFACE = 2
CHART = "C0010"
GENERATOR_EXPONENTS = (1, 0, 0, 0)
ZERO = (0,) * VARIABLE_COUNT
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/"
    "cm_bridge_hs_rank2_k_generator_tree_contraction"
)
TARGET_ARTIFACT = ROOT / (
    "results/section12_instantiation/"
    "cm_bridge_hs_extension_from_surface_presentation/"
    "Z_2_C0010_hs_extension.npz"
)
COUNUNIT_CACHE_SCHEMA = "section12-tree-contraction-counit-pq-cache-v2"


class HSKGeneratorTreeContractionError(RuntimeError):
    """An overlap tree, rectangular map, or localized identity changed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def array_sha256(value) -> str:
    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(value)).tobytes()
    ).hexdigest()


def write_json_atomic(path: Path, value) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def polynomial_arrays(matrix):
    return extension._flatten_polynomial_matrix(matrix)


def polynomial_from_arrays(rows, columns, coefficients, exponents):
    result = {}
    for row, column, coefficient, exponent in zip(
        rows, columns, coefficients, exponents
    ):
        key = (int(row), int(column))
        polynomial = result.setdefault(key, {})
        exponent = tuple(map(int, exponent))
        value = (int(polynomial.get(exponent, 0)) + int(coefficient)) % P
        if value:
            polynomial[exponent] = value
        else:
            polynomial.pop(exponent, None)
        if not polynomial:
            result.pop(key, None)
    return result


def load_stored_matrix(artifact, prefix):
    return polynomial_from_arrays(
        artifact[f"{prefix}_row"],
        artifact[f"{prefix}_column"],
        artifact[f"{prefix}_coefficient"],
        artifact[f"{prefix}_exponent"],
    )


def border_data(basis, multiplication):
    basis = [tuple(map(int, exponent)) for exponent in basis]
    basis_set = set(basis)
    normal_forms = {}
    for variable in range(VARIABLE_COUNT):
        for source, exponent in enumerate(basis):
            product = list(exponent)
            product[variable] += 1
            product = tuple(product)
            if product not in basis_set:
                normal_forms.setdefault(
                    product,
                    np.asarray(multiplication[variable, :, source], dtype=np.int64)
                    % P,
                )
    border = sorted(normal_forms)
    return border, np.stack([normal_forms[item] for item in border]) % P


def chart_data(basis, multiplication, equations):
    border, normal_forms = border_data(basis, multiplication)
    data = {
        "standard": [tuple(map(int, exponent)) for exponent in basis],
        "multiplication": np.asarray(multiplication, dtype=np.int64) % P,
        "border": border,
        "normal_forms": normal_forms,
        "surface_equations": equations,
    }
    quotient, tree, _internal, _border = all_charts.transition_quotient(data)
    data["_transition_quotient"] = quotient
    return data, tree


def build_tree(data, tree):
    dimension = len(data["standard"])
    with all_charts.extension_dimension(dimension):
        contraction = extension.build_tree_contraction(data, tree)
        d2, _cycles = extension.build_d2_tree(data, contraction)
    return contraction, d2


def evaluate_monomial_pair(exponent, values, tangents):
    value = 1
    tangent = 0
    for power, coordinate, derivative in zip(exponent, values, tangents):
        power = int(power)
        coordinate = int(coordinate) % P
        derivative = int(derivative) % P
        if power < 0:
            if not coordinate:
                raise HSKGeneratorTreeContractionError(
                    "a Laurent denominator vanished on the overlap"
                )
            coordinate_power = pow(coordinate, power, P)
            derivative_power = (
                power * pow(coordinate, power - 1, P) * derivative
            ) % P
        else:
            coordinate_power = pow(coordinate, power, P)
            derivative_power = (
                0
                if not power
                else power * pow(coordinate, power - 1, P) * derivative % P
            )
        tangent = (
            tangent * coordinate_power + value * derivative_power
        ) % P
        value = value * coordinate_power % P
    return value, tangent


def evaluate_polynomial_pair(polynomial, values, tangents):
    value = 0
    tangent = 0
    for exponent, coefficient in polynomial.items():
        monomial_value, monomial_tangent = evaluate_monomial_pair(
            exponent, values, tangents
        )
        value = (value + int(coefficient) * monomial_value) % P
        tangent = (tangent + int(coefficient) * monomial_tangent) % P
    return value, tangent


def evaluate_polynomial_matrix(matrix, shape, values, tangents):
    rows = []
    columns = []
    values_out = []
    tangents_out = []
    for (row, column), polynomial in matrix.items():
        value, tangent = evaluate_polynomial_pair(polynomial, values, tangents)
        if value or tangent:
            rows.append(int(row))
            columns.append(int(column))
            values_out.append(value)
            tangents_out.append(tangent)
    value_matrix = sparse.csr_matrix(
        (np.asarray(values_out, dtype=np.int64), (rows, columns)),
        shape=shape,
        dtype=np.int64,
    )
    tangent_matrix = sparse.csr_matrix(
        (np.asarray(tangents_out, dtype=np.int64), (rows, columns)),
        shape=shape,
        dtype=np.int64,
    )
    value_matrix.data %= P
    tangent_matrix.data %= P
    value_matrix.eliminate_zeros()
    tangent_matrix.eliminate_zeros()
    return value_matrix, tangent_matrix


def pair_add(left, right, scale=1):
    value = (left[0] + int(scale) * right[0]).tocsr()
    tangent = (left[1] + int(scale) * right[1]).tocsr()
    value.data %= P
    tangent.data %= P
    value.eliminate_zeros()
    tangent.eliminate_zeros()
    return value, tangent


def pair_multiply(left, right):
    left_value = left[0].astype(np.int64)
    left_tangent = left[1].astype(np.int64)
    right_value = right[0].astype(np.int64)
    right_tangent = right[1].astype(np.int64)
    value = (left_value @ right_value).tocsr()
    tangent = (
        left_tangent @ right_value + left_value @ right_tangent
    ).tocsr()
    value.data %= P
    tangent.data %= P
    value.eliminate_zeros()
    tangent.eliminate_zeros()
    return value, tangent


def pair_zero(pair):
    return pair[0].nnz == 0 and pair[1].nnz == 0


def pair_hash(pair):
    digest = hashlib.sha256()
    for matrix in pair:
        matrix = matrix.tocsr()
        for value in (matrix.indptr, matrix.indices, matrix.data % P):
            digest.update(np.ascontiguousarray(value).tobytes())
    return digest.hexdigest()


def constant_pair(matrix):
    matrix = sparse.csr_matrix(matrix, dtype=np.int64)
    matrix.data %= P
    matrix.eliminate_zeros()
    return matrix, sparse.csr_matrix(matrix.shape, dtype=np.int64)


def evaluate_block_polynomial(polynomial, comparison):
    return k_shamash.evaluate_block_multi(
        polynomial,
        comparison["edge_multiplication"],
        comparison["inverses"],
    )


def globalize_graph_maps_by_base(
    support_results,
    degree,
    target_basis,
    target_affine,
    base_affine,
    edge_hermite,
    edge_hermite_inverse,
):
    wedges = tuple(itertools.combinations(range(VARIABLE_COUNT), degree))
    wedge_index = {wedge: index for index, wedge in enumerate(wedges)}
    target_columns = np.stack(
        [hermite.hermite_column(exponent, target_affine) for exponent in target_basis],
        axis=1,
    ) % P
    edge_dimension = edge_hermite_inverse.shape[0]
    target_dimension = len(target_basis)
    block_keys = sorted(
        {
            key
            for result in support_results
            for key in result["graph_maps"][degree]
        }
    )
    recovery_checks = 0
    matrices = []
    for base_values, base_tangents in base_affine:
        value_rows = []
        value_columns = []
        value_coefficients = []
        tangent_rows = []
        tangent_columns = []
        tangent_coefficients = []
        for output_wedge, input_wedge in block_keys:
            coefficient_value = np.zeros(
                (edge_dimension, edge_dimension), dtype=np.int64
            )
            coefficient_tangent = np.zeros_like(coefficient_value)
            for coefficient_support, result in enumerate(support_results):
                polynomial = result["graph_maps"][degree].get(
                    (output_wedge, input_wedge), {}
                )
                local_value = np.zeros((2, 2), dtype=np.int64)
                local_tangent = np.zeros((2, 2), dtype=np.int64)
                for exponent, coefficient_matrix in polynomial.items():
                    monomial_value, monomial_tangent = evaluate_monomial_pair(
                        exponent, base_values, base_tangents
                    )
                    local_value = (
                        local_value
                        + monomial_value
                        * np.asarray(coefficient_matrix, dtype=np.int64)
                    ) % P
                    local_tangent = (
                        local_tangent
                        + monomial_tangent
                        * np.asarray(coefficient_matrix, dtype=np.int64)
                    ) % P
                start = 2 * coefficient_support
                coefficient_value[start : start + 2, start : start + 2] = (
                    local_value
                )
                coefficient_tangent[start : start + 2, start : start + 2] = (
                    local_tangent
                )
            value_jets = coefficient_value @ target_columns % P
            tangent_jets = coefficient_tangent @ target_columns % P
            value_block = edge_hermite_inverse @ value_jets % P
            tangent_block = edge_hermite_inverse @ tangent_jets % P
            if np.any((edge_hermite @ value_block - value_jets) % P) or np.any(
                (edge_hermite @ tangent_block - tangent_jets) % P
            ):
                raise HSKGeneratorTreeContractionError(
                    "a p/q-separated graph block failed Hermite replay"
                )
            for block, rows, columns, coefficients in (
                (
                    value_block,
                    value_rows,
                    value_columns,
                    value_coefficients,
                ),
                (
                    tangent_block,
                    tangent_rows,
                    tangent_columns,
                    tangent_coefficients,
                ),
            ):
                nonzero_rows, nonzero_columns = np.nonzero(block)
                rows.extend(
                    wedge_index[output_wedge] * edge_dimension + int(row)
                    for row in nonzero_rows
                )
                columns.extend(
                    wedge_index[input_wedge] * target_dimension + int(column)
                    for column in nonzero_columns
                )
                coefficients.extend(
                    int(block[row, column])
                    for row, column in zip(nonzero_rows, nonzero_columns)
                )
            recovery_checks += len(support_results)
        shape = (
            len(wedges) * edge_dimension,
            len(wedges) * target_dimension,
        )
        value_matrix = sparse.csr_matrix(
            (value_coefficients, (value_rows, value_columns)),
            shape=shape,
            dtype=np.int64,
        )
        tangent_matrix = sparse.csr_matrix(
            (tangent_coefficients, (tangent_rows, tangent_columns)),
            shape=shape,
            dtype=np.int64,
        )
        value_matrix.data %= P
        tangent_matrix.data %= P
        value_matrix.eliminate_zeros()
        tangent_matrix.eliminate_zeros()
        matrices.append((value_matrix, tangent_matrix))
    return matrices, recovery_checks, target_columns


def globalize_gamma_by_base(
    support_results,
    base_affine,
    edge_hermite,
    edge_hermite_inverse,
):
    wedges = tuple(itertools.combinations(range(VARIABLE_COUNT), 2))
    wedge_index = {wedge: index for index, wedge in enumerate(wedges)}
    edge_dimension = edge_hermite_inverse.shape[0]
    checks = 0
    matrices = []
    for base_values, base_tangents in base_affine:
        value_matrix = np.zeros(
            (len(wedges) * edge_dimension, EQUATION_COUNT), dtype=np.int64
        )
        tangent_matrix = np.zeros_like(value_matrix)
        for equation in range(EQUATION_COUNT):
            block_keys = sorted(
                {
                    output_wedge
                    for result in support_results
                    for (output_wedge, input_wedge) in result["gamma"][equation]
                    if input_wedge == ()
                }
            )
            for output_wedge in block_keys:
                value_jets = np.zeros(edge_dimension, dtype=np.int64)
                tangent_jets = np.zeros(edge_dimension, dtype=np.int64)
                for coefficient_support, result in enumerate(support_results):
                    polynomial = result["gamma"][equation].get(
                        (output_wedge, ()), {}
                    )
                    local_value = np.zeros((2, 2), dtype=np.int64)
                    local_tangent = np.zeros((2, 2), dtype=np.int64)
                    for exponent, coefficient_matrix in polynomial.items():
                        monomial_value, monomial_tangent = evaluate_monomial_pair(
                            exponent, base_values, base_tangents
                        )
                        local_value = (
                            local_value
                            + monomial_value
                            * np.asarray(coefficient_matrix, dtype=np.int64)
                        ) % P
                        local_tangent = (
                            local_tangent
                            + monomial_tangent
                            * np.asarray(coefficient_matrix, dtype=np.int64)
                        ) % P
                    start = 2 * coefficient_support
                    value_jets[start : start + 2] = local_value[:, 0]
                    tangent_jets[start : start + 2] = local_tangent[:, 0]
                value_coefficients = edge_hermite_inverse @ value_jets % P
                tangent_coefficients = edge_hermite_inverse @ tangent_jets % P
                if np.any(
                    (edge_hermite @ value_coefficients - value_jets) % P
                ) or np.any(
                    (edge_hermite @ tangent_coefficients - tangent_jets) % P
                ):
                    raise HSKGeneratorTreeContractionError(
                        "a p/q-separated Gamma column failed Hermite replay"
                    )
                start = wedge_index[output_wedge] * edge_dimension
                value_matrix[
                    start : start + edge_dimension, equation
                ] = value_coefficients
                tangent_matrix[
                    start : start + edge_dimension, equation
                ] = tangent_coefficients
                checks += len(support_results)
        matrices.append(
            (
                sparse.csr_matrix(value_matrix % P, dtype=np.int64),
                sparse.csr_matrix(tangent_matrix % P, dtype=np.int64),
            )
        )
    return matrices, checks


def globalize_local_block_maps_by_base(
    local_maps,
    degree_out,
    degree_in,
    input_columns,
    base_affine,
    output_hermite,
    output_hermite_inverse,
):
    """Globalize q-local block maps while p evaluates their exponents."""

    output_wedges = tuple(
        itertools.combinations(range(VARIABLE_COUNT), degree_out)
    )
    input_wedges = tuple(
        itertools.combinations(range(VARIABLE_COUNT), degree_in)
    )
    output_index = {
        wedge: index for index, wedge in enumerate(output_wedges)
    }
    input_index = {
        wedge: index for index, wedge in enumerate(input_wedges)
    }
    output_dimension = output_hermite_inverse.shape[0]
    input_dimension = input_columns.shape[1]
    block_keys = sorted({key for matrix in local_maps for key in matrix})
    matrices = []
    checks = 0
    for base_values, base_tangents in base_affine:
        value_rows = []
        value_columns = []
        value_data = []
        tangent_rows = []
        tangent_columns = []
        tangent_data = []
        for output_wedge, input_wedge in block_keys:
            coefficient_value = np.zeros(
                (output_dimension, output_dimension), dtype=np.int64
            )
            coefficient_tangent = np.zeros_like(coefficient_value)
            for coefficient_support, matrix in enumerate(local_maps):
                polynomial = matrix.get((output_wedge, input_wedge), {})
                local_value = np.zeros((2, 2), dtype=np.int64)
                local_tangent = np.zeros((2, 2), dtype=np.int64)
                for exponent, coefficient_matrix in polynomial.items():
                    monomial_value, monomial_tangent = evaluate_monomial_pair(
                        exponent, base_values, base_tangents
                    )
                    coefficient_matrix = np.asarray(
                        coefficient_matrix, dtype=np.int64
                    )
                    local_value = (
                        local_value + monomial_value * coefficient_matrix
                    ) % P
                    local_tangent = (
                        local_tangent + monomial_tangent * coefficient_matrix
                    ) % P
                start = 2 * coefficient_support
                coefficient_value[start : start + 2, start : start + 2] = (
                    local_value
                )
                coefficient_tangent[
                    start : start + 2, start : start + 2
                ] = local_tangent
            value_jets = coefficient_value @ input_columns % P
            tangent_jets = coefficient_tangent @ input_columns % P
            value_block = output_hermite_inverse @ value_jets % P
            tangent_block = output_hermite_inverse @ tangent_jets % P
            if np.any(
                (output_hermite @ value_block - value_jets) % P
            ) or np.any(
                (output_hermite @ tangent_block - tangent_jets) % P
            ):
                raise HSKGeneratorTreeContractionError(
                    "a p/q-separated block map failed Hermite replay"
                )
            for block, rows, columns, data in (
                (
                    value_block,
                    value_rows,
                    value_columns,
                    value_data,
                ),
                (
                    tangent_block,
                    tangent_rows,
                    tangent_columns,
                    tangent_data,
                ),
            ):
                nonzero_rows, nonzero_columns = np.nonzero(block)
                rows.extend(
                    output_index[output_wedge] * output_dimension + int(row)
                    for row in nonzero_rows
                )
                columns.extend(
                    input_index[input_wedge] * input_dimension + int(column)
                    for column in nonzero_columns
                )
                data.extend(
                    int(block[row, column])
                    for row, column in zip(nonzero_rows, nonzero_columns)
                )
            checks += len(local_maps)
        shape = (
            len(output_wedges) * output_dimension,
            len(input_wedges) * input_dimension,
        )
        value = sparse.csr_matrix(
            (value_data, (value_rows, value_columns)),
            shape=shape,
            dtype=np.int64,
        )
        tangent = sparse.csr_matrix(
            (tangent_data, (tangent_rows, tangent_columns)),
            shape=shape,
            dtype=np.int64,
        )
        value.data %= P
        tangent.data %= P
        value.eliminate_zeros()
        tangent.eliminate_zeros()
        matrices.append((value, tangent))
    return matrices, checks


def globalize_equation_maps_by_base(
    support_results,
    key,
    degree_out,
    degree_in,
    input_columns,
    base_affine,
    output_hermite,
    output_hermite_inverse,
):
    equation_matrices = []
    total_checks = 0
    for equation in range(EQUATION_COUNT):
        matrices, checks = globalize_local_block_maps_by_base(
            [result[key][equation] for result in support_results],
            degree_out,
            degree_in,
            input_columns,
            base_affine,
            output_hermite,
            output_hermite_inverse,
        )
        equation_matrices.append(matrices)
        total_checks += checks
    return equation_matrices, total_checks


def evaluate_row_polynomial_pair(polynomial, values, tangents):
    """Evaluate a Laurent row polynomial without collapsing its 1x2 block."""

    value = np.zeros((1, 2), dtype=np.int64)
    tangent = np.zeros((1, 2), dtype=np.int64)
    for exponent, coefficient_row in polynomial.items():
        monomial_value, monomial_tangent = evaluate_monomial_pair(
            exponent, values, tangents
        )
        coefficient_row = np.asarray(coefficient_row, dtype=np.int64) % P
        value = (value + monomial_value * coefficient_row) % P
        tangent = (tangent + monomial_tangent * coefficient_row) % P
    return value, tangent


def formal_source_lambda_rows(shamash_data, graph_data, adapted_affine):
    """Build formal source lambda rows before any p=q evaluation."""

    orientation = all_charts.orientation_map(SURFACE)
    orientation_scalar = int(
        shamash_data["normalization"][
            "adapted_source_orientation_scalar"
        ]
    )
    rows = []
    for point_row, adapted_jet in zip(
        graph_data["point_rows"], adapted_affine
    ):
        source_support = tuple(json.loads(point_row["source_support"]))
        delta0, delta1 = orientation[source_support]
        residue = [
            orientation_scalar * int(delta1) % P,
            orientation_scalar * int(delta0) % P,
        ]
        values, tangents = adapted_jet
        rows.append(
            k_shamash.lambda_transport._lambda_row(
                shamash_data["source_equations"],
                values,
                tangents,
                residue,
            )
        )
    return rows


def globalize_source_lambda_by_base(
    local_rows,
    adapted_affine,
    edge_hermite,
):
    """Apply p/q separation to the formal source lambda cochain."""

    pairs = tuple(itertools.combinations(range(VARIABLE_COUNT), 2))
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    edge_dimension = edge_hermite.shape[0]
    matrices = []
    checks = 0
    for base_values, base_tangents in adapted_affine:
        value_jets = np.zeros((len(pairs), edge_dimension), dtype=np.int64)
        tangent_jets = np.zeros_like(value_jets)
        for coefficient_support, row in enumerate(local_rows):
            start = 2 * coefficient_support
            for wedge, polynomial in row.items():
                value, tangent = evaluate_row_polynomial_pair(
                    polynomial, base_values, base_tangents
                )
                block = pair_index[wedge]
                value_jets[block, start : start + 2] = value
                tangent_jets[block, start : start + 2] = tangent
                checks += 1
        value_standard = value_jets @ edge_hermite % P
        tangent_standard = tangent_jets @ edge_hermite % P
        matrices.append(
            (
                sparse.csr_matrix(
                    value_standard.reshape(1, -1), dtype=np.int64
                ),
                sparse.csr_matrix(
                    tangent_standard.reshape(1, -1), dtype=np.int64
                ),
            )
        )
    return matrices, checks


def monomial_pairs_for_unique_exponents(unique_exponents, affine):
    values = np.zeros(
        (len(affine), len(unique_exponents)), dtype=np.int64
    )
    tangents = np.zeros_like(values)
    for base, (base_values, base_tangents) in enumerate(affine):
        for exponent_number, exponent in enumerate(unique_exponents):
            value, tangent = evaluate_monomial_pair(
                exponent, base_values, base_tangents
            )
            values[base, exponent_number] = value
            tangents[base, exponent_number] = tangent
    return values, tangents


def vectorized_serialized_row_by_base(
    columns,
    coefficients,
    exponents,
    column_count,
    affine,
):
    """Evaluate a serialized sparse polynomial row at all dual factors."""

    columns = np.asarray(columns, dtype=np.int64)
    coefficients = np.asarray(coefficients, dtype=np.int64) % P
    exponents = np.asarray(exponents, dtype=np.int16)
    minima = exponents.min(axis=0).astype(np.int64)
    spans = (
        exponents.max(axis=0).astype(np.int64) - minima + 1
    )
    multipliers = np.ones(exponents.shape[1], dtype=np.int64)
    for index in range(1, exponents.shape[1]):
        multipliers[index] = multipliers[index - 1] * spans[index - 1]
    if int(multipliers[-1]) * int(spans[-1]) >= 2**63:
        raise HSKGeneratorTreeContractionError(
            "the serialized exponent key exceeds int64"
        )
    packed = (
        (exponents.astype(np.int64) - minima) * multipliers
    ).sum(axis=1)
    _keys, first, inverse = np.unique(
        packed, return_index=True, return_inverse=True
    )
    unique_exponents = exponents[first]
    monomial_values, monomial_tangents = (
        monomial_pairs_for_unique_exponents(unique_exponents, affine)
    )
    matrices = []
    for base in range(len(affine)):
        value_data = (
            coefficients * monomial_values[base, inverse]
        ) % P
        tangent_data = (
            coefficients * monomial_tangents[base, inverse]
        ) % P
        # bincount is exact here: every integer sum is far below 2^53.
        value_row = np.bincount(
            columns,
            weights=value_data,
            minlength=column_count,
        ).astype(np.int64) % P
        tangent_row = np.bincount(
            columns,
            weights=tangent_data,
            minlength=column_count,
        ).astype(np.int64) % P
        matrices.append(
            (
                sparse.csr_matrix(value_row.reshape(1, -1)),
                sparse.csr_matrix(tangent_row.reshape(1, -1)),
            )
        )
    return matrices, int(len(unique_exponents))


def transform_target_lambda(
    lambda_jet_by_base,
    target_hermite,
    normalization_pairs,
):
    """Scale the output p dual pair while preserving every q/module column."""

    matrices = []
    target_dimension = target_hermite.shape[0]
    if len(lambda_jet_by_base) != len(normalization_pairs):
        raise HSKGeneratorTreeContractionError(
            "the target-lambda p rows and normalization ledger differ"
        )
    for (value_jet, tangent_jet), normalization_pair in zip(
        lambda_jet_by_base, normalization_pairs
    ):
        scale_value, scale_tangent = map(
            lambda item: int(item) % P, normalization_pair
        )
        value_blocks = (
            value_jet.toarray().reshape(28, target_dimension)
            @ target_hermite
        ) % P
        tangent_blocks = (
            tangent_jet.toarray().reshape(28, target_dimension)
            @ target_hermite
        ) % P
        normalized_value_blocks = scale_value * value_blocks % P
        normalized_tangent_blocks = (
            scale_tangent * value_blocks
            + scale_value * tangent_blocks
        ) % P
        value = np.zeros((1, 28 * target_dimension + EQUATION_COUNT), dtype=np.int64)
        tangent = np.zeros_like(value)
        value[0, : 28 * target_dimension] = (
            normalized_value_blocks.reshape(-1)
        )
        tangent[0, : 28 * target_dimension] = (
            normalized_tangent_blocks.reshape(-1)
        )
        matrices.append(
            (
                sparse.csr_matrix(value),
                sparse.csr_matrix(tangent),
            )
        )
    return matrices


def append_equation_zeros(row_pairs, equation_count=EQUATION_COUNT):
    result = []
    for value, tangent in row_pairs:
        zero = sparse.csr_matrix(
            (1, int(equation_count)), dtype=np.int64
        )
        result.append(
            (
                sparse.hstack((value, zero), format="csr"),
                sparse.hstack((tangent, zero), format="csr"),
            )
        )
    return result


def pair_rows_to_dense(row_pairs):
    return (
        np.stack([pair[0].toarray()[0] for pair in row_pairs]) % P,
        np.stack([pair[1].toarray()[0] for pair in row_pairs]) % P,
    )


def dense_to_pair_rows(values, tangents):
    return [
        (
            sparse.csr_matrix(value.reshape(1, -1), dtype=np.int64),
            sparse.csr_matrix(tangent.reshape(1, -1), dtype=np.int64),
        )
        for value, tangent in zip(values, tangents)
    ]


def pair_left_scalar(pair, scalar):
    value = (int(scalar[0]) * pair[0]).tocsr()
    tangent = (
        int(scalar[1]) * pair[0] + int(scalar[0]) * pair[1]
    ).tocsr()
    value.data %= P
    tangent.data %= P
    value.eliminate_zeros()
    tangent.eliminate_zeros()
    return value, tangent


def assemble_l0(weight, a_row, j1):
    value = sparse.bmat(
        [
            [
                sparse.csr_matrix([[int(weight[0])]], dtype=np.int64),
                (-a_row[0]).tocsr(),
            ],
            [
                sparse.csr_matrix((j1[0].shape[0], 1), dtype=np.int64),
                j1[0],
            ],
        ],
        format="csr",
    )
    tangent = sparse.bmat(
        [
            [
                sparse.csr_matrix([[int(weight[1])]], dtype=np.int64),
                (-a_row[1]).tocsr(),
            ],
            [
                sparse.csr_matrix((j1[1].shape[0], 1), dtype=np.int64),
                j1[1],
            ],
        ],
        format="csr",
    )
    value.data %= P
    tangent.data %= P
    value.eliminate_zeros()
    tangent.eliminate_zeros()
    return value, tangent


def exact_modular_rank_minor(matrix):
    """Return exact independent rows and pivot columns over F_31."""

    source = np.asarray(matrix, dtype=np.int64) % P
    work = source.copy()
    row_ids = np.arange(work.shape[0], dtype=np.int64)
    pivot_columns = []
    pivot_rows = []
    row = 0
    for column in range(work.shape[1]):
        candidates = np.flatnonzero(work[row:, column])
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        if pivot != row:
            work[[row, pivot]] = work[[pivot, row]]
            row_ids[[row, pivot]] = row_ids[[pivot, row]]
        inverse = pow(int(work[row, column]), -1, P)
        work[row, column:] = work[row, column:] * inverse % P
        below = np.flatnonzero(work[row + 1 :, column]) + row + 1
        if len(below):
            multipliers = work[below, column].copy()
            work[below, column:] = (
                work[below, column:]
                - multipliers[:, None] * work[row, column:]
            ) % P
        pivot_columns.append(column)
        pivot_rows.append(int(row_ids[row]))
        row += 1
        if row == work.shape[0]:
            break
    pivot_rows = np.asarray(pivot_rows, dtype=np.int64)
    pivot_columns = np.asarray(pivot_columns, dtype=np.int64)
    minor = source[np.ix_(pivot_rows, pivot_columns)]
    inverse = hermite.inverse_mod(minor)
    return pivot_rows, pivot_columns, inverse


def rowspace_witness(basis, inverse, pivot_columns, residual):
    """Return y with basis*y=0 and residual*y != 0."""

    pivot_set = set(map(int, pivot_columns))
    for column in np.flatnonzero(residual):
        column = int(column)
        if column in pivot_set:
            continue
        coefficients = (
            -inverse @ basis[:, column].reshape(-1, 1)
        ) % P
        witness = np.zeros((basis.shape[1], 1), dtype=np.int64)
        witness[pivot_columns, 0] = coefficients[:, 0]
        witness[column, 0] = 1
        null_residual = basis @ witness % P
        pairing = int(
            (np.asarray(residual).reshape(1, -1) @ witness % P).item()
        )
        if not np.any(null_residual) and pairing:
            return witness, pairing
    raise HSKGeneratorTreeContractionError(
        "failed to extract a rowspace nonmembership witness"
    )


def persist_rowspace_obstruction(kind, support, witness, pairing, ledger):
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    witness_path = OUTPUT_DIRECTORY / "counit_rowspace_obstruction.npz"
    np.savez_compressed(
        witness_path,
        field=np.asarray([P], dtype=np.uint8),
        support=np.asarray([support], dtype=np.uint16),
        kind=np.asarray(kind),
        witness=np.asarray(witness, dtype=np.uint8),
        nonzero_pairing=np.asarray([pairing], dtype=np.uint8),
    )
    record = dict(ledger)
    record.update(
        {
            "schema": "section12-tree-contraction-counit-obstruction-v1",
            "support": int(support),
            "component": kind,
            "nonzero_pairing_mod_31": int(pairing),
            "witness_file": witness_path.name,
            "witness_sha256": sha256(witness_path),
        }
    )
    write_json_atomic(
        OUTPUT_DIRECTORY / "counit_rowspace_obstruction.json", record
    )


def solve_row_pair(
    matrix_pair,
    rhs_pair,
    *,
    support,
    preferred_minor=None,
):
    """Solve ``a*matrix=rhs`` over one retained dual overlap factor.

    This is a direct exact solve.  The chosen structural minor is inverted
    modulo 31 and both the value and tangent equations are substituted back
    into every column.
    """

    d0 = matrix_pair[0].tocsr().astype(np.int64)
    d1 = matrix_pair[1].tocsr().astype(np.int64)
    r0 = rhs_pair[0].toarray().astype(np.int64) % P
    r1 = rhs_pair[1].toarray().astype(np.int64) % P
    d0_dense = d0.toarray().astype(np.int64) % P
    d1_dense = d1.toarray().astype(np.int64) % P
    pivot_rows = pivot_columns = inverse = None
    if preferred_minor is not None:
        candidate_rows, candidate_columns = preferred_minor
        candidate = d0_dense[np.ix_(candidate_rows, candidate_columns)]
        try:
            candidate_inverse = hermite.inverse_mod(candidate)
        except Exception:
            candidate_inverse = None
        if candidate_inverse is not None:
            pivot_rows = np.asarray(candidate_rows, dtype=np.int64)
            pivot_columns = np.asarray(candidate_columns, dtype=np.int64)
            inverse = candidate_inverse
    if inverse is None:
        pivot_rows, pivot_columns, inverse = exact_modular_rank_minor(
            d0_dense
        )
    rank = len(pivot_rows)
    basis = d0_dense[pivot_rows, :]
    coefficients0 = r0[:, pivot_columns] @ inverse % P
    a0 = np.zeros((1, d0.shape[0]), dtype=np.int64)
    a0[:, pivot_rows] = coefficients0
    value_residual = (a0 @ d0_dense - r0) % P
    if np.any(value_residual):
        witness, pairing = rowspace_witness(
            basis, inverse, pivot_columns, value_residual
        )
        ledger = {
            "value_rank": int(rank),
            "free_directions": int(d0.shape[0] - rank),
            "value_membership": False,
            "tangent_membership": None,
            "pivot_rows_sha256": array_sha256(pivot_rows),
            "pivot_columns_sha256": array_sha256(pivot_columns),
        }
        persist_rowspace_obstruction(
            "value", support, witness, pairing, ledger
        )
        raise HSKGeneratorTreeContractionError(
            "the whole-overlap counit value row is outside rowspace(D2); "
            "an exact witness was persisted"
        )
    tangent_rhs = (
        r1 - a0 @ d1_dense
    ) % P
    coefficients1 = tangent_rhs[:, pivot_columns] @ inverse % P
    a1 = np.zeros((1, d0.shape[0]), dtype=np.int64)
    a1[:, pivot_rows] = coefficients1
    tangent_residual = (
        a1 @ d0_dense - tangent_rhs
    ) % P
    if np.any(tangent_residual):
        witness, pairing = rowspace_witness(
            basis, inverse, pivot_columns, tangent_residual
        )
        ledger = {
            "value_rank": int(rank),
            "free_directions": int(d0.shape[0] - rank),
            "value_membership": True,
            "tangent_membership": False,
            "pivot_rows_sha256": array_sha256(pivot_rows),
            "pivot_columns_sha256": array_sha256(pivot_columns),
        }
        persist_rowspace_obstruction(
            "tangent", support, witness, pairing, ledger
        )
        raise HSKGeneratorTreeContractionError(
            "the whole-overlap counit tangent row is outside rowspace(D2); "
            "an exact witness was persisted"
        )
    result = (
        sparse.csr_matrix(a0, dtype=np.int64),
        sparse.csr_matrix(a1, dtype=np.int64),
    )
    replay = pair_add(pair_multiply(result, matrix_pair), rhs_pair, scale=-1)
    if not pair_zero(replay):
        raise HSKGeneratorTreeContractionError(
            "the exact whole-overlap counit solve failed full substitution"
        )
    return result, {
        "value_rank": int(rank),
        "free_directions": int(d0.shape[0] - rank),
        "canonical_free_variables_zero": True,
        "pivot_rows_sha256": array_sha256(pivot_rows),
        "pivot_columns_sha256": array_sha256(pivot_columns),
        "pivot_minor_inverse_sha256": array_sha256(inverse),
        "value_membership": True,
        "tangent_membership": True,
        "full_replay_sha256": pair_hash(replay),
        "_minor": (pivot_rows, pivot_columns),
    }


def raw_d1_pair(multiplication, values, tangents):
    dimension = int(np.asarray(multiplication).shape[1])
    identity = sparse.eye(dimension, format="csr", dtype=np.int64)
    value_blocks = []
    tangent_blocks = []
    for variable in range(VARIABLE_COUNT):
        value_block = (
            int(values[variable]) * identity
            - sparse.csr_matrix(
                np.asarray(multiplication[variable], dtype=np.int64)
            )
        )
        value_block.data %= P
        tangent_block = int(tangents[variable]) * identity
        tangent_block.data %= P
        value_blocks.append(value_block)
        tangent_blocks.append(tangent_block)
    value = sparse.hstack(value_blocks, format="csr")
    tangent = sparse.hstack(tangent_blocks, format="csr")
    value.data %= P
    tangent.data %= P
    value.eliminate_zeros()
    tangent.eliminate_zeros()
    return value, tangent


def raw_d3_pair(multiplication, values, tangents):
    dimension = int(np.asarray(multiplication).shape[1])
    pairs = tuple(itertools.combinations(range(VARIABLE_COUNT), 2))
    triples = tuple(itertools.combinations(range(VARIABLE_COUNT), 3))
    pair_numbers = {pair: index for index, pair in enumerate(pairs)}
    identity = sparse.eye(dimension, format="csr", dtype=np.int64)
    value_blocks = [[None for _ in triples] for _ in pairs]
    tangent_blocks = [[None for _ in triples] for _ in pairs]
    for triple_number, (a, b, c) in enumerate(triples):
        for pair, variable, sign in (
            ((b, c), a, 1),
            ((a, c), b, -1),
            ((a, b), c, 1),
        ):
            value = (
                int(values[variable]) * identity
                - sparse.csr_matrix(
                    np.asarray(multiplication[variable], dtype=np.int64)
                )
            )
            tangent = int(tangents[variable]) * identity
            if sign == -1:
                value = -value
                tangent = -tangent
            value.data %= P
            tangent.data %= P
            value_blocks[pair_numbers[pair]][triple_number] = value
            tangent_blocks[pair_numbers[pair]][triple_number] = tangent
    value = sparse.bmat(value_blocks, format="csr")
    tangent = sparse.bmat(tangent_blocks, format="csr")
    value.data %= P
    tangent.data %= P
    value.eliminate_zeros()
    tangent.eliminate_zeros()
    return value, tangent


def generator_row_matrix(generators):
    matrix = {}
    for column, polynomial in enumerate(generators):
        if polynomial:
            matrix[(0, int(column))] = dict(polynomial)
    return matrix


def assemble_j3(graph_j3, omega_columns, b_j1_columns):
    value_blocks = [[None for _ in range(7)] for _ in range(7)]
    tangent_blocks = [[None for _ in range(7)] for _ in range(7)]
    value_blocks[0][0] = graph_j3[0]
    tangent_blocks[0][0] = graph_j3[1]
    for equation in range(EQUATION_COUNT):
        value_blocks[0][equation + 1] = omega_columns[equation][0]
        tangent_blocks[0][equation + 1] = omega_columns[equation][1]
        value_blocks[equation + 1][equation + 1] = (
            b_j1_columns[equation][0]
        )
        tangent_blocks[equation + 1][equation + 1] = (
            b_j1_columns[equation][1]
        )
    value = sparse.bmat(value_blocks, format="csr")
    tangent = sparse.bmat(tangent_blocks, format="csr")
    value.data %= P
    tangent.data %= P
    value.eliminate_zeros()
    tangent.eliminate_zeros()
    return value, tangent


def csr_records(matrix):
    matrix = matrix.tocoo()
    order = np.lexsort((matrix.col, matrix.row))
    return np.column_stack(
        (
            np.asarray(matrix.row)[order],
            np.asarray(matrix.col)[order],
            np.asarray(matrix.data)[order] % P,
        )
    ).astype(np.int32)


def build():
    if not TARGET_ARTIFACT.is_file():
        raise HSKGeneratorTreeContractionError(
            "the independently normalized C0010 tree artifact is missing"
        )
    shamash_data = k_shamash.build()
    graph_data = k_j1.build_data()
    if tuple(shamash_data["record"]["exponents"]) != GENERATOR_EXPONENTS:
        raise HSKGeneratorTreeContractionError("the K generator changed")

    adapted_affine = []
    target_affine = []
    for point_row in graph_data["point_rows"]:
        source_values = tuple(json.loads(point_row["source_affine_values"]))
        source_tangents = tuple(json.loads(point_row["source_affine_tangent"]))
        adapted_affine.append(
            k_shamash.adapt_jet(
                source_values,
                source_tangents,
                shamash_data["factor_data"],
            )
        )
        target_affine.append(
            (
                tuple(json.loads(point_row["target_affine_values"])),
                tuple(json.loads(point_row["target_affine_tangent"])),
            )
        )

    edge_basis, edge_hermite, tested, final_degree = hermite.select_standard_monomials(
        adapted_affine
    )
    edge_multiplication, _edge_border = hermite.border_multiplication(
        edge_basis, edge_hermite, adapted_affine
    )
    edge_inverse = hermite.inverse_mod(edge_hermite)
    edge_data, edge_tree = chart_data(
        edge_basis, edge_multiplication, shamash_data["source_equations"]
    )
    edge_contraction, edge_d2 = build_tree(edge_data, edge_tree)
    edge_i1 = extension.i1_matrix(edge_contraction["i1"])

    target_data, target_tree = chart_data(
        graph_data["target_basis"],
        graph_data["target_multiplication"],
        k_shamash.shamash.surface_equations(
            k_shamash.iw_solver.load_problem_context(), SURFACE, CHART
        )[0],
    )
    target_contraction, target_d2_rebuilt = build_tree(target_data, target_tree)
    with np.load(TARGET_ARTIFACT) as artifact:
        target_d2 = load_stored_matrix(artifact, "D2_tree")
        target_i1 = load_stored_matrix(artifact, "tree_K1_inclusion_i1")
        stored_p1 = np.asarray(artifact["tree_projection_p1"], dtype=np.int64) % P
        stored_ranks = tuple(map(int, artifact["tree_complex_ranks"]))
        target_hermite = np.asarray(
            artifact["hermite_standard_to_jets"], dtype=np.int64
        ) % P
        target_lambda_rows = np.asarray(
            artifact["lambda_K2_jet_row"], dtype=np.int64
        )
        target_lambda_columns = np.asarray(
            artifact["lambda_K2_jet_column"], dtype=np.int64
        )
        target_lambda_coefficients = np.asarray(
            artifact["lambda_K2_jet_coefficient"], dtype=np.int64
        )
        target_lambda_exponents = np.asarray(
            artifact["lambda_K2_jet_exponent"], dtype=np.int16
        )
    if target_d2 != target_d2_rebuilt:
        raise HSKGeneratorTreeContractionError(
            "the independently rebuilt target D2 differs from the stored target"
        )
    if not np.array_equal(stored_p1, target_contraction["p1"] % P):
        raise HSKGeneratorTreeContractionError(
            "the independently rebuilt target tree ordering changed"
        )

    edge_dimension = len(edge_basis)
    target_dimension = len(graph_data["target_basis"])
    edge_f1 = 7 * edge_dimension + 1
    edge_f2 = 28 * edge_dimension + EQUATION_COUNT
    target_f1 = 7 * target_dimension + 1
    target_f2 = 28 * target_dimension + EQUATION_COUNT
    if edge_dimension != 2 * len(adapted_affine):
        raise HSKGeneratorTreeContractionError(
            "the overlap is not the product of its retained dual factors"
        )
    if target_dimension != 98:
        raise HSKGeneratorTreeContractionError(
            "the independently fixed C0010 target length changed"
        )
    if stored_ranks[:3] != (1, target_f1, target_f2):
        raise HSKGeneratorTreeContractionError("the stored target ranks changed")

    graph_j1_by_base, graph1_checks, target_columns = globalize_graph_maps_by_base(
        shamash_data["support_results"],
        1,
        graph_data["target_basis"],
        target_affine,
        adapted_affine,
        edge_hermite,
        edge_inverse,
    )
    graph_j2_by_base, graph2_checks, _ = globalize_graph_maps_by_base(
        shamash_data["support_results"],
        2,
        graph_data["target_basis"],
        target_affine,
        adapted_affine,
        edge_hermite,
        edge_inverse,
    )
    graph_j3_by_base, graph3_checks, _ = globalize_graph_maps_by_base(
        shamash_data["support_results"],
        3,
        graph_data["target_basis"],
        target_affine,
        adapted_affine,
        edge_hermite,
        edge_inverse,
    )
    gamma_by_base, gamma_checks = globalize_gamma_by_base(
        shamash_data["support_results"],
        adapted_affine,
        edge_hermite,
        edge_inverse,
    )
    omega_by_equation, omega_checks = globalize_equation_maps_by_base(
        shamash_data["support_results"],
        "omega",
        3,
        1,
        target_columns,
        adapted_affine,
        edge_hermite,
        edge_inverse,
    )
    edge_h1_by_equation, edge_h1_checks = (
        globalize_equation_maps_by_base(
            shamash_data["support_results"],
            "edge_h1",
            2,
            1,
            edge_hermite,
            adapted_affine,
            edge_hermite,
            edge_inverse,
        )
    )
    left_mixed_by_equation = []
    left_mixed_checks = 0
    for equation in range(EQUATION_COUNT):
        local_maps = []
        for result in shamash_data["support_results"]:
            left = shamash.add_block_maps(
                shamash.compose_block_maps(
                    result["graph_maps"][2], result["target_h1"][equation]
                ),
                shamash.compose_block_maps(
                    result["gamma"][equation],
                    result["target_differentials"][1],
                ),
            )
            local_maps.append(left)
        matrices, checks = globalize_local_block_maps_by_base(
            local_maps,
            2,
            1,
            target_columns,
            adapted_affine,
            edge_hermite,
            edge_inverse,
        )
        left_mixed_by_equation.append(matrices)
        left_mixed_checks += checks

    p1_edge = sparse.csr_matrix(edge_contraction["p1"] % P, dtype=np.int64)
    graph_j0 = sparse.csr_matrix(
        edge_inverse @ target_columns % P, dtype=np.int64
    )
    edge_unit = edge_basis.index(ZERO)
    target_unit = graph_data["target_basis"].index(ZERO)
    target_unit_column = sparse.csr_matrix(
        (
            np.asarray([1], dtype=np.int64),
            (
                np.asarray([target_unit], dtype=np.int64),
                np.asarray([0], dtype=np.int64),
            ),
        ),
        shape=(target_dimension, 1),
    )
    expected_edge_unit = sparse.csr_matrix(
        (
            np.asarray([1], dtype=np.int64),
            (
                np.asarray([edge_unit], dtype=np.int64),
                np.asarray([0], dtype=np.int64),
            ),
        ),
        shape=(edge_dimension, 1),
    )
    unit_value_residual = (
        graph_j0 @ target_unit_column - expected_edge_unit
    ).tocsr()
    unit_value_residual.data %= P
    unit_value_residual.eliminate_zeros()
    unit_tangent_residual = sparse.csr_matrix(
        (edge_dimension, 1), dtype=np.int64
    )
    if unit_value_residual.nnz:
        raise HSKGeneratorTreeContractionError(
            "the degree-zero graph map does not preserve the unit"
        )

    cache_key = hashlib.sha256()
    cache_key.update(COUNUNIT_CACHE_SCHEMA.encode("ascii"))
    cache_key.update(sha256(TARGET_ARTIFACT).encode("ascii"))
    for value in (
        edge_hermite,
        target_hermite,
        np.asarray(adapted_affine, dtype=np.int16),
        np.asarray(target_affine, dtype=np.int16),
    ):
        cache_key.update(np.ascontiguousarray(value).tobytes())
    cache_key.update(
        np.ascontiguousarray(
            np.asarray(
                shamash_data["target_lambda_normalization_pairs"],
                dtype=np.uint8,
            )
        ).tobytes()
    )
    cache_key = cache_key.hexdigest()
    counit_cache_path = OUTPUT_DIRECTORY / "counit_pq_cache.npz"
    cache_loaded = False
    if counit_cache_path.is_file():
        with np.load(counit_cache_path) as cache:
            cached_key = str(np.asarray(cache["cache_key"]).item())
            if cached_key == cache_key:
                source_lambda_by_base = dense_to_pair_rows(
                    cache["source_lambda_value"],
                    cache["source_lambda_tangent"],
                )
                target_lambda_by_base = dense_to_pair_rows(
                    cache["target_lambda_value"],
                    cache["target_lambda_tangent"],
                )
                source_lambda_checks = int(cache["source_lambda_checks"][0])
                target_lambda_unique_exponents = int(
                    cache["target_lambda_unique_exponents"][0]
                )
                cache_loaded = True
    if not cache_loaded:
        source_formal_lambda = formal_source_lambda_rows(
            shamash_data, graph_data, adapted_affine
        )
        source_lambda_by_base, source_lambda_checks = (
            globalize_source_lambda_by_base(
                source_formal_lambda, adapted_affine, edge_hermite
            )
        )
        source_lambda_by_base = append_equation_zeros(source_lambda_by_base)
        if np.any(target_lambda_rows):
            raise HSKGeneratorTreeContractionError(
                "the stored target lambda is no longer a row matrix"
            )
        target_lambda_jet_by_base, target_lambda_unique_exponents = (
            vectorized_serialized_row_by_base(
                target_lambda_columns,
                target_lambda_coefficients,
                target_lambda_exponents,
                28 * target_dimension,
                target_affine,
            )
        )
        target_lambda_by_base = transform_target_lambda(
            target_lambda_jet_by_base,
            target_hermite,
            shamash_data["target_lambda_normalization_pairs"],
        )
        source_value, source_tangent = pair_rows_to_dense(
            source_lambda_by_base
        )
        target_value, target_tangent = pair_rows_to_dense(
            target_lambda_by_base
        )
        OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            counit_cache_path,
            cache_key=np.asarray(cache_key),
            source_lambda_value=np.asarray(source_value, dtype=np.uint8),
            source_lambda_tangent=np.asarray(source_tangent, dtype=np.uint8),
            target_lambda_value=np.asarray(target_value, dtype=np.uint8),
            target_lambda_tangent=np.asarray(target_tangent, dtype=np.uint8),
            source_lambda_checks=np.asarray(
                [source_lambda_checks], dtype=np.uint32
            ),
            target_lambda_unique_exponents=np.asarray(
                [target_lambda_unique_exponents], dtype=np.uint32
            ),
        )

    support_outputs = []
    identity_hashes = []
    d1_identity_hashes = []
    j1_inclusion_identity_hashes = []
    top_identity_hashes = []
    full_a_identity_hashes = []
    a_solve_ledgers = []
    a_rows = []
    a_minor_template = None
    graph_j3_identity_hashes = []
    generator_identity_hashes = []
    edge_generator_definition_hashes = []
    target_generator_definition_hashes = []
    mixed_j3_identity_hashes = []
    contracted_mixed_j3_identity_hashes = []
    lower_j3_identity_hashes = []
    j3_factor_hashes = []
    for support, (
        adapted_jet,
        target_jet,
        comparison,
        graph_j1,
            graph_j2,
            gamma,
            lambda_edge,
            lambda_target,
            graph_j3,
        ) in enumerate(
        zip(
            adapted_affine,
            target_affine,
            shamash_data["support_results"],
            graph_j1_by_base,
            graph_j2_by_base,
            gamma_by_base,
            source_lambda_by_base,
            target_lambda_by_base,
            graph_j3_by_base,
        )
    ):
        adapted_values, adapted_tangents = adapted_jet
        target_values, target_tangents = target_jet
        d1_edge = raw_d1_pair(
            edge_multiplication, adapted_values, adapted_tangents
        )
        d1_target = raw_d1_pair(
            graph_data["target_multiplication"],
            target_values,
            target_tangents,
        )
        d1_residual = pair_add(
            pair_multiply(d1_edge, graph_j1),
            pair_multiply(constant_pair(graph_j0), d1_target),
            scale=-1,
        )
        if not pair_zero(d1_residual):
            value_coo = d1_residual[0].tocoo()
            tangent_coo = d1_residual[1].tocoo()
            raise HSKGeneratorTreeContractionError(
                f"the p/q-separated raw D1 identity failed at base {support}: "
                f"({d1_residual[0].nnz},{d1_residual[1].nnz}); "
                f"value={list(zip(value_coo.row.tolist(), value_coo.col.tolist(), (value_coo.data % P).tolist()))}; "
                f"tangent={list(zip(tangent_coo.row.tolist(), tangent_coo.col.tolist(), (tangent_coo.data % P).tolist()))}"
            )
        d1_identity_hashes.append(pair_hash(d1_residual))
        i1_target = evaluate_polynomial_matrix(
            target_i1,
            (VARIABLE_COUNT * target_dimension, target_f1),
            target_values,
            target_tangents,
        )
        j1_naive = pair_multiply(
            pair_multiply(constant_pair(p1_edge), graph_j1), i1_target
        )
        i1_edge = evaluate_polynomial_matrix(
            edge_i1,
            (VARIABLE_COUNT * edge_dimension, edge_f1),
            adapted_values,
            adapted_tangents,
        )
        j1_replay = pair_add(
            pair_multiply(i1_edge, j1_naive),
            pair_multiply(graph_j1, i1_target),
            scale=-1,
        )
        if not pair_zero(j1_replay):
            raise HSKGeneratorTreeContractionError(
                f"the tree J1 inclusion replay failed at base {support}: "
                f"({j1_replay[0].nnz},{j1_replay[1].nnz})"
            )
        j1_inclusion_identity_hashes.append(pair_hash(j1_replay))

        d3_edge = raw_d3_pair(
            edge_multiplication, adapted_values, adapted_tangents
        )
        d3_target = raw_d3_pair(
            graph_data["target_multiplication"],
            target_values,
            target_tangents,
        )
        graph_j3_residual = pair_add(
            pair_multiply(d3_edge, graph_j3),
            pair_multiply(graph_j2, d3_target),
            scale=-1,
        )
        if not pair_zero(graph_j3_residual):
            raise HSKGeneratorTreeContractionError(
                f"the rectangular raw J3 graph identity failed at {support}: "
                f"({graph_j3_residual[0].nnz},"
                f"{graph_j3_residual[1].nnz})"
            )
        graph_j3_identity_hashes.append(pair_hash(graph_j3_residual))

        j2_value = sparse.lil_matrix((edge_f2, target_f2), dtype=np.int64)
        j2_tangent = sparse.lil_matrix((edge_f2, target_f2), dtype=np.int64)
        j2_value[: 28 * edge_dimension, : 28 * target_dimension] = graph_j2[0]
        j2_tangent[: 28 * edge_dimension, : 28 * target_dimension] = graph_j2[1]
        j2_value[: 28 * edge_dimension, 28 * target_dimension :] = gamma[0]
        j2_tangent[: 28 * edge_dimension, 28 * target_dimension :] = gamma[1]
        for equation, entry in enumerate(shamash_data["b_entries"]):
            value, tangent = evaluate_polynomial_pair(
                entry["polynomial"], adapted_values, adapted_tangents
            )
            j2_value[28 * edge_dimension + equation, 28 * target_dimension + equation] = value
            j2_tangent[28 * edge_dimension + equation, 28 * target_dimension + equation] = tangent
        j2 = (j2_value.tocsr(), j2_tangent.tocsr())
        j2[0].data %= P
        j2[1].data %= P

        d2_edge = evaluate_polynomial_matrix(
            edge_d2,
            (edge_f1, edge_f2),
            adapted_values,
            adapted_tangents,
        )
        d2_target = evaluate_polynomial_matrix(
            target_d2,
            (target_f1, target_f2),
            target_values,
            target_tangents,
        )
        j1 = j1_naive
        residual = pair_add(
            pair_multiply(d2_edge, j2),
            pair_multiply(j1, d2_target),
            scale=-1,
        )
        if not pair_zero(residual):
            graph_value = residual[0][:, : 28 * target_dimension]
            equation_value = residual[0][:, 28 * target_dimension :]
            graph_tangent = residual[1][:, : 28 * target_dimension]
            equation_tangent = residual[1][:, 28 * target_dimension :]
            raise HSKGeneratorTreeContractionError(
                f"D2_edge J2-J1 D2_target is nonzero at support {support}: "
                f"value_nnz={residual[0].nnz}, tangent_nnz={residual[1].nnz}; "
                f"graph=({graph_value.nnz},{graph_tangent.nnz}), "
                f"equations=({equation_value.nnz},{equation_tangent.nnz}); "
                f"J1_replay=({j1_replay[0].nnz},{j1_replay[1].nnz})"
            )

        target_generator = evaluate_polynomial_matrix(
            generator_row_matrix(target_contraction["generators"]),
            (1, target_f1),
            target_values,
            target_tangents,
        )
        edge_generator = evaluate_polynomial_matrix(
            generator_row_matrix(edge_contraction["generators"]),
            (1, edge_f1),
            adapted_values,
            adapted_tangents,
        )
        edge_unit_column = sparse.csr_matrix(
            (
                np.asarray([1], dtype=np.int64),
                (
                    np.asarray([edge_unit], dtype=np.int64),
                    np.asarray([0], dtype=np.int64),
                ),
            ),
            shape=(edge_dimension, 1),
        )
        target_unit_column_definition = sparse.csr_matrix(
            (
                np.asarray([1], dtype=np.int64),
                (
                    np.asarray([target_unit], dtype=np.int64),
                    np.asarray([0], dtype=np.int64),
                ),
            ),
            shape=(target_dimension, 1),
        )
        edge_generator_definition = pair_add(
            pair_multiply(d1_edge, i1_edge),
            pair_multiply(
                constant_pair(edge_unit_column), edge_generator
            ),
            scale=-1,
        )
        target_generator_definition = pair_add(
            pair_multiply(d1_target, i1_target),
            pair_multiply(
                constant_pair(target_unit_column_definition),
                target_generator,
            ),
            scale=-1,
        )
        if not pair_zero(edge_generator_definition) or not pair_zero(
            target_generator_definition
        ):
            raise HSKGeneratorTreeContractionError(
                "a tree generator stopped agreeing with d1*i1"
            )
        edge_generator_definition_hashes.append(
            pair_hash(edge_generator_definition)
        )
        target_generator_definition_hashes.append(
            pair_hash(target_generator_definition)
        )
        generator_residual = pair_add(
            pair_multiply(edge_generator, j1),
            target_generator,
            scale=-1,
        )
        if not pair_zero(generator_residual):
            raise HSKGeneratorTreeContractionError(
                f"the tree generator/J1 identity failed at {support}: "
                f"({generator_residual[0].nnz},"
                f"{generator_residual[1].nnz})"
            )
        generator_identity_hashes.append(pair_hash(generator_residual))

        omega_columns = []
        b_j1_columns = []
        b_scalar_pairs = []
        mixed_support_hash = hashlib.sha256()
        contracted_support_hash = hashlib.sha256()
        lower_support_hash = hashlib.sha256()
        for equation, entry in enumerate(shamash_data["b_entries"]):
            omega = omega_by_equation[equation][support]
            omega_column = pair_multiply(omega, i1_target)
            omega_columns.append(omega_column)
            b_value, b_tangent = evaluate_polynomial_pair(
                entry["polynomial"], adapted_values, adapted_tangents
            )
            b_pair = (b_value, b_tangent)
            b_scalar_pairs.append(b_pair)
            b_j1 = pair_left_scalar(j1, b_pair)
            b_j1_columns.append(b_j1)
            right_top = pair_add(
                pair_multiply(d3_edge, omega),
                pair_left_scalar(
                    pair_multiply(
                        edge_h1_by_equation[equation][support],
                        graph_j1,
                    ),
                    b_pair,
                ),
            )
            mixed_residual = pair_add(
                left_mixed_by_equation[equation][support],
                right_top,
                scale=-1,
            )
            if not pair_zero(mixed_residual):
                raise HSKGeneratorTreeContractionError(
                    f"the factored mixed J3 identity failed at support "
                    f"{support}, equation {equation}: "
                    f"({mixed_residual[0].nnz},"
                    f"{mixed_residual[1].nnz})"
                )
            h1_edge_i1_edge = pair_multiply(
                edge_h1_by_equation[equation][support], i1_edge
            )
            actual_upper_left = pair_add(
                pair_multiply(d3_edge, omega_column),
                pair_multiply(h1_edge_i1_edge, b_j1),
            )
            actual_upper_right = pair_multiply(
                left_mixed_by_equation[equation][support], i1_target
            )
            contracted_residual = pair_add(
                actual_upper_left, actual_upper_right, scale=-1
            )
            if not pair_zero(contracted_residual):
                raise HSKGeneratorTreeContractionError(
                    "the actual contracted upper D3/J3 block failed"
                )
            mixed_support_hash.update(pair_hash(mixed_residual).encode("ascii"))
            contracted_support_hash.update(
                pair_hash(contracted_residual).encode("ascii")
            )
            lower_residual = pair_add(
                pair_multiply(edge_generator, b_j1),
                pair_left_scalar(target_generator, b_pair),
                scale=-1,
            )
            if not pair_zero(lower_residual):
                raise HSKGeneratorTreeContractionError(
                    "the lower mixed J3 block failed"
                )
            lower_support_hash.update(
                pair_hash(lower_residual).encode("ascii")
            )
        j3 = assemble_j3(graph_j3, omega_columns, b_j1_columns)
        if j3[0].shape != (
            56 * edge_dimension + EQUATION_COUNT * edge_f1,
            56 * target_dimension + EQUATION_COUNT * target_f1,
        ):
            raise HSKGeneratorTreeContractionError(
                "the rectangular J3 shape changed"
            )
        mixed_j3_identity_hashes.append(mixed_support_hash.hexdigest())
        contracted_mixed_j3_identity_hashes.append(
            contracted_support_hash.hexdigest()
        )
        lower_j3_identity_hashes.append(lower_support_hash.hexdigest())
        j3_factor_hashes.append(pair_hash(j3))
        del j3

        # Recompute the whole p/q finite-overlap counit row.  Only after this
        # replay may the finite-overlap a row be certified.
        w_value, w_tangent = evaluate_monomial_pair(
            shamash_data["weight_exponent"], adapted_values, adapted_tangents
        )
        w_value = w_value * int(shamash_data["weight_scalar"]) % P
        w_tangent = w_tangent * int(shamash_data["weight_scalar"]) % P
        top_residual = pair_add(
            pair_multiply(lambda_edge, j2),
            pair_left_scalar(lambda_target, (w_value, w_tangent)),
            scale=-1,
        )
        top_identity_hashes.append(pair_hash(top_residual))
        if pair_zero(top_residual):
            a_value = sparse.csr_matrix((1, target_f1), dtype=np.int64)
            a_tangent = sparse.csr_matrix((1, target_f1), dtype=np.int64)
            a_ledger = {
                "value_rank": None,
                "free_directions": None,
                "value_membership": True,
                "tangent_membership": True,
                "zero_rhs": True,
            }
        else:
            (a_value, a_tangent), a_ledger = solve_row_pair(
                d2_target,
                top_residual,
                support=support,
                preferred_minor=a_minor_template,
            )
            a_minor_template = a_ledger.pop("_minor")
            a_ledger["zero_rhs"] = False
        a_rows.append((a_value, a_tangent))
        a_solve_ledgers.append(a_ledger)
        a_d2_residual = pair_add(
            top_residual,
            pair_multiply((a_value, a_tangent), d2_target),
            scale=-1,
        )
        if not pair_zero(a_d2_residual):
            raise HSKGeneratorTreeContractionError(
                "the certified finite-overlap a row failed its D2 replay"
            )
        full_a_identity_hashes.append(pair_hash(a_d2_residual))
        l0_value, l0_tangent = assemble_l0(
            (w_value, w_tangent),
            (a_value, a_tangent),
            j1,
        )
        support_outputs.append(
            {
                "J1": j1,
                "J2": j2,
                "J3_G3": graph_j3,
                "J3_Omega_i1": omega_columns,
                "J3_B": b_scalar_pairs,
                "L0": (l0_value, l0_tangent),
                "L1": j2,
                "a": (a_value, a_tangent),
                "D2_residual": residual,
                "top_residual": top_residual,
                "a_D2_residual": a_d2_residual,
                "naive_J1_raw_replay_residual": j1_replay,
            }
        )
        identity_hashes.append(pair_hash(residual))

    return {
        "shamash_data": shamash_data,
        "graph_data": graph_data,
        "adapted_affine": adapted_affine,
        "target_affine": target_affine,
        "edge_basis": edge_basis,
        "edge_hermite": edge_hermite,
        "edge_inverse": edge_inverse,
        "edge_multiplication": edge_multiplication,
        "target_lambda_normalization_pairs": np.asarray(
            shamash_data["target_lambda_normalization_pairs"],
            dtype=np.uint8,
        ),
        "residue_transport_pairs": np.asarray(
            shamash_data["residue_transport_pairs"], dtype=np.uint8
        ),
        "raw_counit_weight_pairs": np.asarray(
            shamash_data["raw_counit_weight_pairs"], dtype=np.uint8
        ),
        "orientation_coboundary_pairs": np.asarray(
            shamash_data["orientation_coboundary_pairs"],
            dtype=np.uint8,
        ),
        "edge_contraction": edge_contraction,
        "edge_d2": edge_d2,
        "edge_i1": edge_i1,
        "target_contraction": target_contraction,
        "target_d2": target_d2,
        "target_i1": target_i1,
        "graph_j1_by_base": graph_j1_by_base,
        "graph_j2_by_base": graph_j2_by_base,
        "gamma_by_base": gamma_by_base,
        "support_outputs": support_outputs,
        "identity_hashes": identity_hashes,
        "d1_identity_hashes": d1_identity_hashes,
        "j1_inclusion_identity_hashes": j1_inclusion_identity_hashes,
        "top_identity_hashes": top_identity_hashes,
        "full_a_identity_hashes": full_a_identity_hashes,
        "a_solve_ledgers": a_solve_ledgers,
        "graph_j3_identity_hashes": graph_j3_identity_hashes,
        "generator_identity_hashes": generator_identity_hashes,
        "edge_generator_definition_hashes": edge_generator_definition_hashes,
        "target_generator_definition_hashes": (
            target_generator_definition_hashes
        ),
        "mixed_j3_identity_hashes": mixed_j3_identity_hashes,
        "contracted_mixed_j3_identity_hashes": (
            contracted_mixed_j3_identity_hashes
        ),
        "lower_j3_identity_hashes": lower_j3_identity_hashes,
        "j3_factor_hashes": j3_factor_hashes,
        "unit_identity_hash": pair_hash(
            (unit_value_residual, unit_tangent_residual)
        ),
        "counts": {
            "overlap_supports": len(adapted_affine),
            "overlap_length": edge_dimension,
            "target_length": target_dimension,
            "overlap_tree_ranks": [1, edge_f1, edge_f2, 56 * edge_dimension + 6 * edge_f1],
            "target_tree_ranks": list(stored_ranks),
            "J1_shape": [edge_f1, target_f1],
            "J2_shape": [edge_f2, target_f2],
            "J3_shape": [
                56 * edge_dimension + EQUATION_COUNT * edge_f1,
                56 * target_dimension + EQUATION_COUNT * target_f1,
            ],
            "L0_shape": [edge_f1 + 1, target_f1 + 1],
            "graph_J1_recovery_checks": graph1_checks,
            "graph_J2_recovery_checks": graph2_checks,
            "graph_J3_recovery_checks": graph3_checks,
            "Gamma_recovery_checks": gamma_checks,
            "Omega_recovery_checks": omega_checks,
            "edge_H1_recovery_checks": edge_h1_checks,
            "left_mixed_recovery_checks": left_mixed_checks,
            "raw_graph_J3_identities": len(graph_j3_identity_hashes),
            "tree_generator_J1_identities": len(generator_identity_hashes),
            "edge_d1_i1_generator_definition_identities": len(
                edge_generator_definition_hashes
            ),
            "target_d1_i1_generator_definition_identities": len(
                target_generator_definition_hashes
            ),
            "upper_mixed_J3_identities": (
                EQUATION_COUNT * len(mixed_j3_identity_hashes)
            ),
            "contracted_upper_mixed_J3_identities": (
                EQUATION_COUNT * len(contracted_mixed_j3_identity_hashes)
            ),
            "actual_lower_mixed_J3_identities": (
                EQUATION_COUNT * len(lower_j3_identity_hashes)
            ),
            "factorized_D3_chain_block_identities": (
                (1 + 2 * EQUATION_COUNT) * len(graph_j3_identity_hashes)
            ),
            "formal_source_lambda_pq_checks": source_lambda_checks,
            "target_lambda_unique_exponents": (
                target_lambda_unique_exponents
            ),
            "target_lambda_normalization_pairs_sha256": (
                k_shamash.pair_ledger_sha256(
                    shamash_data["target_lambda_normalization_pairs"]
                )
            ),
            "tree_chain_identities": len(identity_hashes),
            "edge_standard_monomials_tested": tested,
            "edge_highest_selected_degree": final_degree,
        },
    }


def write_artifacts(output_directory=OUTPUT_DIRECTORY):
    global OUTPUT_DIRECTORY
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    previous_output_directory = OUTPUT_DIRECTORY
    OUTPUT_DIRECTORY = output_directory
    try:
        result = build()
    finally:
        OUTPUT_DIRECTORY = previous_output_directory
    maps_path = output_directory / "tree_contracted_support_maps.npz"
    hermite_path = output_directory / "overlap_adapted_hermite.npz"
    certificate_path = output_directory / "tree_contraction_certificate.json"
    counit_cache_path = output_directory / "counit_pq_cache.npz"

    payload = {
        "field": np.asarray([P], dtype=np.uint8),
        "generator_exponents": np.asarray(GENERATOR_EXPONENTS, dtype=np.uint8),
        "J1_shape": np.asarray(result["counts"]["J1_shape"], dtype=np.uint32),
        "J2_shape": np.asarray(result["counts"]["J2_shape"], dtype=np.uint32),
        "J3_shape": np.asarray(result["counts"]["J3_shape"], dtype=np.uint32),
        "L0_shape": np.asarray(result["counts"]["L0_shape"], dtype=np.uint32),
        "target_lambda_normalization_pairs": result[
            "target_lambda_normalization_pairs"
        ],
        "residue_transport_pairs": result["residue_transport_pairs"],
        "raw_counit_weight_pairs": result[
            "raw_counit_weight_pairs"
        ],
        "orientation_coboundary_pairs": result[
            "orientation_coboundary_pairs"
        ],
    }
    for support, output in enumerate(result["support_outputs"]):
        for name in ("J1", "J2", "J3_G3", "L0", "L1", "a"):
            for component, matrix in zip(("value", "tangent"), output[name]):
                payload[f"s{support:02d}_{name}_{component}"] = csr_records(matrix)
        for equation, omega_i1 in enumerate(output["J3_Omega_i1"]):
            for component, matrix in zip(("value", "tangent"), omega_i1):
                payload[
                    f"s{support:02d}_J3_Omega_i1_e{equation}_{component}"
                ] = csr_records(matrix)
        payload[f"s{support:02d}_J3_B_value"] = np.asarray(
            [pair[0] for pair in output["J3_B"]], dtype=np.uint8
        )
        payload[f"s{support:02d}_J3_B_tangent"] = np.asarray(
            [pair[1] for pair in output["J3_B"]], dtype=np.uint8
        )
    np.savez_compressed(maps_path, **payload)
    np.savez_compressed(
        hermite_path,
        field=np.asarray([P], dtype=np.uint8),
        adapted_standard_monomials=np.asarray(result["edge_basis"], dtype=np.uint8),
        adapted_hermite_matrix=np.asarray(result["edge_hermite"], dtype=np.uint8),
        adapted_hermite_inverse=np.asarray(result["edge_inverse"], dtype=np.uint8),
        adapted_multiplication_matrices=np.asarray(
            result["edge_multiplication"], dtype=np.uint8
        ),
        tree_edges=np.asarray(result["edge_contraction"]["tree_edges"], dtype=np.uint16),
        non_tree_edges=np.asarray(
            result["edge_contraction"]["non_tree_edges"], dtype=np.uint16
        ),
        tree_projection_p1=np.asarray(result["edge_contraction"]["p1"], dtype=np.uint8),
    )
    support_count = int(result["counts"]["overlap_supports"])
    edge_dimension = int(result["counts"]["overlap_length"])
    j1_shape = tuple(map(int, result["counts"]["J1_shape"]))
    j2_shape = tuple(map(int, result["counts"]["J2_shape"]))
    j3_shape = tuple(map(int, result["counts"]["J3_shape"]))
    l0_shape = tuple(map(int, result["counts"]["L0_shape"]))
    nonzero_a_supports = sum(
        bool(output["a"][0].nnz or output["a"][1].nnz)
        for output in result["support_outputs"]
    )
    target_rank_values = sorted(
        {int(ledger["value_rank"]) for ledger in result["a_solve_ledgers"]}
    )
    free_direction_values = sorted(
        {
            int(ledger["free_directions"])
            for ledger in result["a_solve_ledgers"]
        }
    )
    certificate = {
        "schema": "section12-HS-rank2-K-generator-tree-contraction-v3",
        "field": "F_31",
        "status": (
            "one-generator-rectangular-tree-J1-J2-J3-L0-exact-on-finite-overlap;"
            "whole-overlap-counit-needs-nonzero-a;"
            "strict-Laurent-counit-lift-not-inferred"
        ),
        "surface": "Z_2",
        "chart": CHART,
        "K_generator_exponents": list(GENERATOR_EXPONENTS),
        "map_direction": (
            "independently normalized target C0010 to source C0010 restricted "
            f"to the length-{edge_dimension} K overlap"
        ),
        "counts": result["counts"],
        "generator_normalization": {
            "support_pairs": support_count,
            "target_lambda_normalization_pairs_sha256": (
                k_shamash.pair_ledger_sha256(
                    result["target_lambda_normalization_pairs"]
                )
            ),
            "residue_transport_pairs_sha256": (
                k_shamash.pair_ledger_sha256(
                    result["residue_transport_pairs"]
                )
            ),
            "raw_counit_weight_pairs_sha256": (
                k_shamash.pair_ledger_sha256(
                    result["raw_counit_weight_pairs"]
                )
            ),
            "orientation_coboundary_pairs_sha256": (
                k_shamash.pair_ledger_sha256(
                    result["orientation_coboundary_pairs"]
                )
            ),
            "pair_ledger_hash_encoding": (
                "signed-little-endian-int64-C-order"
            ),
            "orientation_coboundary_law": (
                "rho=A*(u_source/u_target)*J8"
            ),
            "raw_counit_weight_law": "q=rho/det(B)",
            "target_normalization_law": "n=q/w_H",
            "strict_Laurent_lift": "OPEN",
        },
        "exact_checks": {
            "adapted_overlap_Hermite_isomorphism_rank": edge_dimension,
            "adapted_overlap_Hermite_inverse_exact": bool(
                np.array_equal(
                    result["edge_hermite"] @ result["edge_inverse"] % P,
                    np.eye(edge_dimension, dtype=np.int64),
                )
            ),
            "independently_rebuilt_target_D2_matches_stored_target": True,
            "degree_zero_unit_preserved_with_zero_tangent": True,
            "raw_D1_edge_J1_equals_J0_D1_target_supports": len(
                result["d1_identity_hashes"]
            ),
            "tree_inclusion_i1_edge_J1_equals_G1_i1_target_supports": len(
                result["j1_inclusion_identity_hashes"]
            ),
            "rectangular_D2_edge_J2_equals_J1_D2_target_supports": len(
                result["identity_hashes"]
            ),
            "all_tree_chain_residuals_zero": True,
            "whole_overlap_top_delta_nonzero_supports": sum(
                not pair_zero(output["top_residual"])
                for output in result["support_outputs"]
            ),
            "target_D2_value_rank_all_supports": sorted(
                {
                    int(ledger["value_rank"])
                    for ledger in result["a_solve_ledgers"]
                }
            ),
            "target_D2_free_directions_all_supports": sorted(
                {
                    int(ledger["free_directions"])
                    for ledger in result["a_solve_ledgers"]
                }
            ),
            "whole_overlap_delta_value_rowspace_memberships": sum(
                bool(ledger["value_membership"])
                for ledger in result["a_solve_ledgers"]
            ),
            "whole_overlap_delta_tangent_rowspace_memberships": sum(
                bool(ledger["tangent_membership"])
                for ledger in result["a_solve_ledgers"]
            ),
            "finite_overlap_lambda_edge_J2_minus_w_lambda_target_equals_a_D2_supports": len(
                result["full_a_identity_hashes"]
            ),
            "full_L0_A_target_equals_A_edge_L1_blockwise_supports": len(
                result["full_a_identity_hashes"]
            ),
            "rectangular_raw_graph_D3_edge_G3_equals_G2_D3_target_supports": len(
                result["graph_j3_identity_hashes"]
            ),
            "tree_generator_edge_J1_equals_target_supports": len(
                result["generator_identity_hashes"]
            ),
            "edge_tree_generator_equals_raw_d1_i1_supports": len(
                result["edge_generator_definition_hashes"]
            ),
            "target_tree_generator_equals_raw_d1_i1_supports": len(
                result["target_generator_definition_hashes"]
            ),
            "upper_mixed_D3_J3_block_identities": (
                EQUATION_COUNT * len(result["mixed_j3_identity_hashes"])
            ),
            "upper_mixed_identities_after_target_i1": (
                EQUATION_COUNT
                * len(result["contracted_mixed_j3_identity_hashes"])
            ),
            "actual_lower_mixed_D3_J3_block_identities": (
                EQUATION_COUNT * len(result["lower_j3_identity_hashes"])
            ),
            "factorized_rectangular_D3_J3_equals_J2_D3_block_identities": (
                (1 + 2 * EQUATION_COUNT)
                * len(result["graph_j3_identity_hashes"])
            ),
            "full_assembled_D3_products_materialized": False,
            "reason_full_D3_products_not_materialized": (
                "the exact 1+6+6 block factorization is replayed directly; "
                f"materializing {support_count} multi-gigabyte full products "
                "adds no "
                "independent identity"
            ),
            "rectangular_map_is_not_asserted_invertible": True,
        },
        "residual_hashes": {
            "degree_zero_unit": result["unit_identity_hash"],
            "D1_edge_J1_minus_J0_D1_target": result[
                "d1_identity_hashes"
            ],
            "i1_edge_J1_minus_G1_i1_target": result[
                "j1_inclusion_identity_hashes"
            ],
            "D2_edge_J2_minus_J1_D2_target": result["identity_hashes"],
            "lambda_edge_J2_minus_w_lambda_target": result[
                "top_identity_hashes"
            ],
            "lambda_edge_J2_minus_w_lambda_target_minus_a_D2_target": result[
                "full_a_identity_hashes"
            ],
            "raw_D3_edge_G3_minus_G2_D3_target": result[
                "graph_j3_identity_hashes"
            ],
            "tree_generator_edge_J1_minus_target": result[
                "generator_identity_hashes"
            ],
            "edge_raw_d1_i1_minus_tree_generator": result[
                "edge_generator_definition_hashes"
            ],
            "target_raw_d1_i1_minus_tree_generator": result[
                "target_generator_definition_hashes"
            ],
            "upper_mixed_D3_J3_blocks": result["mixed_j3_identity_hashes"],
            "upper_mixed_after_target_i1": result[
                "contracted_mixed_j3_identity_hashes"
            ],
            "actual_lower_mixed_D3_J3_blocks": result[
                "lower_j3_identity_hashes"
            ],
            "factorized_J3": result["j3_factor_hashes"],
        },
        "counit": {
            "old_diagonal_p_equals_q_solution": (
                "a=0 on each isolated diagonal support calculation"
            ),
            "whole_overlap_solution": (
                f"a is nonzero on {nonzero_a_supports} of {support_count} "
                "factors after exact p/q separation"
            ),
            "target_D2_value_rank": target_rank_values,
            "canonical_free_directions_set_to_zero": (
                free_direction_values
            ),
            "value_and_tangent_membership": (
                f"PASS on all {support_count} factors"
            ),
            "full_target_column_replay": (
                f"PASS on all {support_count} factors across "
                f"{j2_shape[1]} target columns"
            ),
            "a_value_nnz_range": [
                min(output["a"][0].nnz for output in result["support_outputs"]),
                max(output["a"][0].nnz for output in result["support_outputs"]),
            ],
            "a_tangent_nnz_range": [
                min(output["a"][1].nnz for output in result["support_outputs"]),
                max(output["a"][1].nnz for output in result["support_outputs"]),
            ],
            "strict_Laurent_tree_lift_status": "OPEN",
            "reason": (
                "the finite-overlap p/q solve proves an exact nonzero boundary "
                "row, but it does not by itself construct a strict "
                "Laurent/surface-ring lift"
            ),
        },
        "claim_boundary": {
            "proved": (
                f"for generator {GENERATOR_EXPONENTS}, actual rectangular J1 "
                f"({j1_shape[0]}x{j1_shape[1]}), J2/L1 "
                f"({j2_shape[0]}x{j2_shape[1]}), factorized J3 "
                f"({j3_shape[0]}x{j3_shape[1]}), and L0 "
                f"({l0_shape[0]}x{l0_shape[1]}) "
                "over every factor of the finite "
                f"length-{edge_dimension} overlap algebra; exact D1, "
                "tree-inclusion, D2, unit, "
                "counit/top, and factorized D3 identities; the whole-overlap counit uses a "
                "canonical a row with exact target-D2 rank and free directions "
                "recorded above"
            ),
            "not_claimed": (
                "a strict Laurent/surface-ring counit lift, other K generators, "
                "K coherence, cross-RHom/Yoneda, Question 11.4, or Hodge"
            ),
            "next_exact_construction": (
                "repeat the certified one-generator construction on the remaining "
                "K generators and then impose K coherence; the strict "
                "Laurent/surface-ring lift of a remains a separate scope"
            ),
        },
        "files": {
            "factorized_support_maps_J1_J2_J3_L0": maps_path.name,
            "adapted_overlap_Hermite": hermite_path.name,
            "counit_pq_rows_cache": counit_cache_path.name,
        },
        "inputs": {
            "supportwise_Shamash": {
                "path": str(Path(k_shamash.__file__).resolve().relative_to(ROOT)),
                "sha256": sha256(Path(k_shamash.__file__).resolve()),
            },
            "independent_target": {
                "path": str(TARGET_ARTIFACT.relative_to(ROOT)),
                "sha256": sha256(TARGET_ARTIFACT),
            },
        },
    }
    if (
        certificate["counts"][
            "target_lambda_normalization_pairs_sha256"
        ]
        != certificate["generator_normalization"][
            "target_lambda_normalization_pairs_sha256"
        ]
    ):
        raise HSKGeneratorTreeContractionError(
            "the target-normalization pair ledger has inconsistent hashes"
        )
    write_json_atomic(certificate_path, certificate)
    files = {}
    manifest_paths = [maps_path, hermite_path, certificate_path]
    if counit_cache_path.is_file():
        manifest_paths.append(counit_cache_path)
    for path in manifest_paths:
        files[path.name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    write_json_atomic(
        output_directory / "manifest.json",
        {
            "schema": "section12-HS-rank2-K-generator-tree-contraction-manifest-v3",
            "generator": {
                "path": str(Path(__file__).resolve().relative_to(ROOT)),
                "sha256": sha256(Path(__file__).resolve()),
            },
            "files": files,
        },
    )
    return certificate


def main():
    certificate = write_artifacts()
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
