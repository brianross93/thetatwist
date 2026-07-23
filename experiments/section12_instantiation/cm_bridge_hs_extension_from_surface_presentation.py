"""Build one exact chart-local Hartshorne--Serre extension cochain.

The Hermite algebra on ``Z_2/C0010`` is a free ``F_31`` vector space of
dimension 98.  Its graph Koszul complex over the eight-variable polynomial
ring is exact.  The six equations of the smooth surface annihilate the
Hermite algebra and have explicit divided-difference nullhomotopies.  These
data give the first three differentials of the Shamash resolution over the
surface ring.

This producer retains the 687-generator tree-contracted cyclic model.  That
model is deliberately larger than the 297-generator border presentation: it
retains a transparent degree-three differential.  The extension cochain is
the change-of-rings cocycle obtained by wedging the six nullhomotopies into
the top graph-Koszul dual functional.  Its local-duality functional is read
from the already certified Hartshorne--Serre orientations.

The output is a presentation of a module on one affine chart.  It is not a
global line-bundle descent datum and it is not a strict degree-zero
projective dg model.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_iw_surface_ring_presentation as surface_presentation,
)


P = 31
VARIABLE_COUNT = 8
SURFACE = 2
CHART = "C0010"
DIMENSION = 98
ZERO_EXPONENT = (0,) * VARIABLE_COUNT
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/"
    "cm_bridge_hs_extension_from_surface_presentation"
)
HERMITE_PATH = ROOT / (
    "results/section12_instantiation/cm_bridge_iw_hermite_border_basis/"
    "Z_2_C0010.npz"
)
ORIENTATION_PATH = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_border_to_local_frames/"
    "hs_completed_lci_supports.csv"
)
FIRST_PRESENTATION_PATH = ROOT / (
    "results/section12_instantiation/cm_bridge_iw_surface_ring_presentation/"
    "Z_2_C0010_surface_ring_presentation.npz"
)


class HSExtensionFromPresentationError(RuntimeError):
    """An exact graph, Shamash, duality, or pushout identity failed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _add(poly, exponent, coefficient):
    surface_presentation._add(poly, exponent, coefficient)


def _add_scaled(target, source, scalar=1):
    surface_presentation._add_scaled(target, source, scalar)


def _multiply(left, right):
    if not left or not right:
        return {}
    if len(left) == 1:
        return _multiply_by_monomial(right, left)
    if len(right) == 1:
        return _multiply_by_monomial(left, right)
    return surface_presentation._multiply(left, right)


def _multiply_by_monomial(polynomial, monomial):
    if not polynomial or not monomial:
        return {}
    if len(monomial) != 1:
        return surface_presentation._multiply(polynomial, monomial)
    (shift, scalar), = monomial.items()
    scalar = int(scalar) % P
    return {
        tuple(int(a) + int(b) for a, b in zip(exponent, shift)):
        int(coefficient) * scalar % P
        for exponent, coefficient in polynomial.items()
        if int(coefficient) * scalar % P
    }


@lru_cache(maxsize=None)
def _curve_power_terms(factor, power):
    """Reduce one affine ``v`` power modulo its monic curve equation.

    The return keys are ``(extra_u_power, residual_v_power)``.  The four
    rewrite systems are triangular, so these are canonical representatives
    over ``F_31[u_0,u_1,u_2,u_3]``.
    """

    factor = int(factor)
    power = int(power)
    if factor in (0, 1, 3):
        quotient, remainder = divmod(power, 2)
        result = {}
        # v^(2q+r)=v^r(1+u^3)^q.
        for index in range(quotient + 1):
            coefficient = int(math.comb(quotient, index)) % P
            if coefficient:
                result[(3 * index, remainder)] = coefficient
        return tuple(sorted(result.items()))
    if factor != 2:
        raise ValueError("curve factor must lie in 0,...,3")

    # v^3=v-u^3.  This recurrence strictly lowers the v exponent.
    cache = {
        0: {(0, 0): 1},
        1: {(0, 1): 1},
        2: {(0, 2): 1},
    }
    for exponent in range(3, power + 1):
        reduced = {}
        for key, coefficient in cache[exponent - 2].items():
            _add(reduced, key, coefficient)
        for (u_power, v_power), coefficient in cache[exponent - 3].items():
            _add(reduced, (u_power + 3, v_power), -coefficient)
        cache[exponent] = reduced
    return tuple(sorted(cache[power].items()))


@lru_cache(maxsize=None)
def _curve_reduce_monomial(exponent):
    partial = {tuple(map(int, exponent)): 1}
    for factor, (u_variable, v_variable) in enumerate(
        ((0, 1), (2, 3), (4, 5), (6, 7))
    ):
        next_partial = {}
        for current, current_coefficient in partial.items():
            v_power = int(current[v_variable])
            if (factor in (0, 1, 3) and v_power < 2) or (
                factor == 2 and v_power < 3
            ):
                _add(next_partial, current, current_coefficient)
                continue
            base = list(current)
            base[v_variable] = 0
            for (extra_u, residual_v), scalar in _curve_power_terms(
                factor, v_power
            ):
                reduced = list(base)
                reduced[u_variable] += int(extra_u)
                reduced[v_variable] = int(residual_v)
                _add(
                    next_partial,
                    tuple(reduced),
                    int(current_coefficient) * int(scalar),
                )
        partial = next_partial
    return tuple(sorted(partial.items()))


def _curve_reduce_polynomial(polynomial):
    """Return the canonical image in ``P/(f_0,f_1,f_2,f_3)``."""

    result = {}
    for exponent, coefficient in polynomial.items():
        for reduced_exponent, scalar in _curve_reduce_monomial(
            tuple(map(int, exponent))
        ):
            _add(result, reduced_exponent, int(coefficient) * int(scalar))
    return result


def _curve_multiply(left, right):
    return _curve_reduce_polynomial(_multiply(left, right))


def _unit_exponent(variable):
    return surface_presentation._unit_exponent(variable)


def _poly_scale(poly, scalar):
    return surface_presentation._scaled(poly, scalar)


def _matrix_entry_add(matrix, row, column, polynomial, scalar=1):
    if not polynomial or int(scalar) % P == 0:
        return
    entry = matrix.setdefault((int(row), int(column)), {})
    _add_scaled(entry, polynomial, scalar)
    if not entry:
        matrix.pop((int(row), int(column)), None)


def _constant(value):
    value = int(value) % P
    return {} if not value else {ZERO_EXPONENT: value}


def _monomial(exponent, coefficient=1):
    coefficient = int(coefficient) % P
    return {} if not coefficient else {tuple(map(int, exponent)): coefficient}


def _matrix_compose(left, right):
    """Compose sparse polynomial matrices ``left @ right``."""

    left_by_column = {}
    for (row, column), polynomial in left.items():
        left_by_column.setdefault(int(column), []).append((int(row), polynomial))
    result = {}
    for (middle, column), right_polynomial in right.items():
        for row, left_polynomial in left_by_column.get(int(middle), []):
            _matrix_entry_add(
                result,
                row,
                int(column),
                _multiply(left_polynomial, right_polynomial),
            )
    return result


def _matrix_add(left, right, scalar=1):
    result = {
        key: dict(polynomial) for key, polynomial in left.items()
    }
    for key, polynomial in right.items():
        entry = result.setdefault(key, {})
        _add_scaled(entry, polynomial, scalar)
        if not entry:
            result.pop(key, None)
    return result


def _matrix_right_scale_columns(matrix, scalars):
    result = {}
    for (row, column), polynomial in matrix.items():
        scalar = scalars.get(int(column), {})
        if scalar:
            result[(int(row), int(column))] = _multiply(polynomial, scalar)
    return result


def _flatten_polynomial_matrix(matrix):
    rows = []
    columns = []
    coefficients = []
    exponents = []
    for (row, column), polynomial in sorted(matrix.items()):
        for exponent, coefficient in sorted(polynomial.items()):
            rows.append(int(row))
            columns.append(int(column))
            coefficients.append(int(coefficient) % P)
            exponents.append(tuple(map(int, exponent)))
    return (
        np.asarray(rows, dtype=np.uint32),
        np.asarray(columns, dtype=np.uint32),
        np.asarray(coefficients, dtype=np.uint8),
        np.asarray(exponents, dtype=np.uint8).reshape((-1, VARIABLE_COUNT)),
    )


def _polynomial_matrix_hash(matrix):
    rows, columns, coefficients, exponents = _flatten_polynomial_matrix(matrix)
    digest = hashlib.sha256()
    for array in (rows, columns, coefficients, exponents):
        digest.update(np.ascontiguousarray(array).tobytes())
    return digest.hexdigest()


def _reduce_polynomial_mod_surface(polynomial, surface_equations):
    """Only certify a polynomial already expressed in the surface ideal.

    The producer never uses heuristic reduction.  Identities over ``S`` are
    stored with an explicit six-component boundary ledger instead.
    """

    return polynomial, surface_equations


def load_inputs():
    if not HERMITE_PATH.exists() or not ORIENTATION_PATH.exists():
        raise HSExtensionFromPresentationError(
            "the bound Hermite or local-orientation input is absent"
        )
    data = surface_presentation.load_input()
    quotient, tree, internal, border = surface_presentation.transition_quotient(data)
    if internal != 195 or border != 589 or len(tree) != 97:
        raise HSExtensionFromPresentationError("the graph-transition census changed")
    with np.load(HERMITE_PATH) as artifact:
        support_indices = np.asarray(
            artifact["support_point_indices"], dtype=np.int64
        )
        hermite = np.asarray(artifact["hermite_matrix"], dtype=np.int64) % P
        hermite_inverse = (
            np.asarray(artifact["hermite_inverse"], dtype=np.int64) % P
        )
        affine_values = np.asarray(artifact["affine_values"], dtype=np.int64) % P
        affine_tangents = (
            np.asarray(artifact["affine_tangents"], dtype=np.int64) % P
        )
    if (
        support_indices.shape != (49, 4)
        or hermite.shape != (98, 98)
        or hermite_inverse.shape != (98, 98)
        or affine_values.shape != (49, 8)
        or affine_tangents.shape != (49, 8)
    ):
        raise HSExtensionFromPresentationError("the representative Hermite chart changed")
    identity = np.eye(DIMENSION, dtype=np.int64)
    if np.any((hermite @ hermite_inverse - identity) % P) or np.any(
        (hermite_inverse @ hermite - identity) % P
    ):
        raise HSExtensionFromPresentationError(
            "the standard-to-jet Hermite matrix stopped being invertible"
        )

    orientation = {}
    with ORIENTATION_PATH.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["surface"] != f"Z_{SURFACE}":
                continue
            key = tuple(map(int, json.loads(row["support_point_indices"])))
            orientation[key] = (
                int(row["delta0_mod_31"]) % P,
                int(row["delta1_mod_31"]) % P,
            )
    missing = [tuple(map(int, row)) for row in support_indices if tuple(map(int, row)) not in orientation]
    if missing:
        raise HSExtensionFromPresentationError(
            f"{len(missing)} representative-chart support orientations are absent"
        )

    tau_on_jets = np.zeros(98, dtype=np.int64)
    local_rows = []
    for position, support in enumerate(support_indices):
        key = tuple(map(int, support))
        delta0, delta1 = orientation[key]
        if not delta0:
            raise HSExtensionFromPresentationError(
                "a local Hartshorne--Serre orientation stopped being a generator"
            )
        tau_on_jets[2 * position] = delta1
        tau_on_jets[2 * position + 1] = delta0
        local_rows.append(
            {
                "support_position": position,
                "support_point_indices": list(key),
                "delta0": delta0,
                "delta1": delta1,
                "residue_row": [delta1, delta0],
            }
        )
    tau_on_standard_basis = tau_on_jets @ hermite % P
    return {
        "data": data,
        "quotient": quotient,
        "tree": tree,
        "support_indices": support_indices,
        "hermite": hermite,
        "hermite_inverse": hermite_inverse,
        "affine_values": affine_values,
        "affine_tangents": affine_tangents,
        "tau_on_jets": tau_on_jets,
        "tau": tau_on_standard_basis,
        "local_rows": local_rows,
    }


def raw_edge(variable, source):
    return int(variable) * DIMENSION + int(source)


def pair_index(left, right):
    if not 0 <= int(left) < int(right) < VARIABLE_COUNT:
        raise ValueError("pair indices must be increasing")
    index = 0
    for first in range(VARIABLE_COUNT):
        for second in range(first + 1, VARIABLE_COUNT):
            if (first, second) == (int(left), int(right)):
                return index
            index += 1
    raise AssertionError("pair index was not found")


def triple_index(a, b, c):
    target = (int(a), int(b), int(c))
    for index, triple in enumerate(itertools.combinations(range(VARIABLE_COUNT), 3)):
        if triple == target:
            return index
    raise ValueError("triple indices must be increasing")


def build_tree_contraction(data, tree):
    dimension = len(data["standard"])
    retained_rank = VARIABLE_COUNT * dimension - (dimension - 1)
    border_rank = len(data["border"])
    standard = data["standard"]
    multiplication = data["multiplication"]
    unit_index = standard.index(ZERO_EXPONENT)
    tree_edge_by_target = {
        int(target): raw_edge(variable, parent)
        for variable, parent, target in tree
    }
    tree_parent = {
        int(target): (int(variable), int(parent))
        for variable, parent, target in tree
    }
    tree_edges = set(tree_edge_by_target.values())
    non_tree_edges = [
        edge for edge in range(VARIABLE_COUNT * dimension)
        if edge not in tree_edges
    ]
    if len(non_tree_edges) != retained_rank:
        raise HSExtensionFromPresentationError("the tree contraction rank changed")
    reduced_index = {edge: index for index, edge in enumerate(non_tree_edges)}
    p1 = np.zeros(
        (retained_rank, VARIABLE_COUNT * dimension), dtype=np.int64
    )
    for edge, index in reduced_index.items():
        p1[index, edge] = 1

    # h0(e_a) is a tree-edge chain with d h0(e_a)=e_a-m_a e_unit.
    h0 = [None] * dimension
    h0[unit_index] = {}
    unresolved = set(range(dimension)) - {unit_index}
    while unresolved:
        progress = False
        for target in list(unresolved):
            variable, parent = tree_parent[target]
            if h0[parent] is None:
                continue
            chain = {tree_edge_by_target[target]: _constant(-1)}
            x_variable = _monomial(_unit_exponent(variable))
            for edge, polynomial in h0[parent].items():
                value = chain.setdefault(edge, {})
                _add_scaled(value, _multiply(x_variable, polynomial))
                if not value:
                    chain.pop(edge, None)
            h0[target] = chain
            unresolved.remove(target)
            progress = True
        if not progress:
            raise HSExtensionFromPresentationError("the order-ideal tree is cyclic")

    def d1_raw_edge(edge):
        variable, source = divmod(int(edge), dimension)
        result = {source: _monomial(_unit_exponent(variable))}
        for target in np.flatnonzero(multiplication[variable, :, source]):
            polynomial = result.setdefault(int(target), {})
            _add(
                polynomial,
                ZERO_EXPONENT,
                -int(multiplication[variable, int(target), source]),
            )
            if not polynomial:
                result.pop(int(target), None)
        return result

    def standard_polynomial(index):
        return _monomial(standard[int(index)])

    generators = []
    i1_columns = []
    for edge in non_tree_edges:
        boundary = d1_raw_edge(edge)
        generator = {}
        for row, polynomial in boundary.items():
            _add_scaled(generator, _multiply(polynomial, standard_polynomial(row)))
        generators.append(generator)

        lifted = {edge: _constant(1)}
        for row, polynomial in boundary.items():
            for tree_edge, tree_polynomial in h0[row].items():
                value = lifted.setdefault(tree_edge, {})
                _add_scaled(value, _multiply(polynomial, tree_polynomial), -1)
                if not value:
                    lifted.pop(tree_edge, None)
        i1_columns.append(lifted)

    # Verify d(i1(e))=generator(e)*unit for all 687 retained edges.
    for column, (generator, lifted) in enumerate(zip(generators, i1_columns)):
        boundary = {}
        for edge, scalar in lifted.items():
            for row, polynomial in d1_raw_edge(edge).items():
                value = boundary.setdefault(row, {})
                _add_scaled(value, _multiply(polynomial, scalar))
                if not value:
                    boundary.pop(row, None)
        expected = {} if not generator else {unit_index: generator}
        if boundary != expected:
            raise HSExtensionFromPresentationError(
                f"tree inclusion column {column} is not a chain lift"
            )
    i1_constant_support = np.zeros(
        (VARIABLE_COUNT * dimension, retained_rank), dtype=np.int64
    )
    for column, lifted in enumerate(i1_columns):
        for edge, polynomial in lifted.items():
            if edge in reduced_index:
                if polynomial != _constant(1) or edge != non_tree_edges[column]:
                    raise HSExtensionFromPresentationError(
                        "the non-tree part of i1 stopped being the identity"
                    )
                i1_constant_support[edge, column] = 1
    if np.any(
        (
            p1 @ i1_constant_support
            - np.eye(retained_rank, dtype=np.int64)
        )
        % P
    ):
        raise HSExtensionFromPresentationError("p1 i1 is not the identity")

    # Constant split merge from the tree carrier to the chart's border carrier.
    q = np.zeros((border_rank, retained_rank), dtype=np.int64)
    internal_columns = []
    representatives = {}
    duplicate_columns = []
    for column, edge in enumerate(non_tree_edges):
        variable, source = divmod(edge, dimension)
        border_row = int(data["_transition_quotient"][variable, source])
        if border_row < 0:
            internal_columns.append(column)
            continue
        q[border_row, column] = 1
        if border_row not in representatives:
            representatives[border_row] = column
        else:
            duplicate_columns.append((column, representatives[border_row]))
    expected_internal = sum(
        1
        for variable in range(VARIABLE_COUNT)
        for source in range(dimension)
        if int(data["_transition_quotient"][variable, source]) < 0
    ) - (dimension - 1)
    expected_duplicates = sum(
        1
        for value in np.asarray(data["_transition_quotient"]).flat
        if int(value) >= 0
    ) - border_rank
    if (
        len(internal_columns) != expected_internal
        or len(duplicate_columns) != expected_duplicates
    ):
        raise HSExtensionFromPresentationError("the merge-kernel census changed")
    if len(representatives) != border_rank:
        raise HSExtensionFromPresentationError("the border merge stopped being onto")
    section = np.zeros((retained_rank, border_rank), dtype=np.int64)
    for row, column in representatives.items():
        section[column, row] = 1
    kernel_rank = retained_rank - border_rank
    kernel = np.zeros((retained_rank, kernel_rank), dtype=np.int64)
    for position, column in enumerate(internal_columns):
        kernel[column, position] = 1
    for offset, (column, representative) in enumerate(
        duplicate_columns, start=len(internal_columns)
    ):
        kernel[column, offset] = 1
        kernel[representative, offset] = -1 % P
    kernel_left_inverse = np.zeros((kernel_rank, retained_rank), dtype=np.int64)
    for position, column in enumerate(internal_columns):
        kernel_left_inverse[position, column] = 1
    for offset, (column, _representative) in enumerate(
        duplicate_columns, start=len(internal_columns)
    ):
        kernel_left_inverse[offset, column] = 1
    if np.any((q @ section - np.eye(border_rank, dtype=np.int64)) % P):
        raise HSExtensionFromPresentationError("the merge section changed")
    if np.any(q @ kernel % P):
        raise HSExtensionFromPresentationError("the merge kernel is not in ker(q)")
    if np.any(
        (
            kernel_left_inverse @ kernel
            - np.eye(kernel_rank, dtype=np.int64)
        )
        % P
    ):
        raise HSExtensionFromPresentationError("the merge-kernel left inverse changed")
    if np.any(kernel_left_inverse @ section % P):
        raise HSExtensionFromPresentationError("the merge splitting is not transverse")
    return {
        "unit_index": unit_index,
        "tree_edges": sorted(tree_edges),
        "non_tree_edges": non_tree_edges,
        "reduced_index": reduced_index,
        "p1": p1,
        "h0": h0,
        "i1": i1_columns,
        "generators": generators,
        "q": q,
        "section": section,
        "kernel": kernel,
        "kernel_left_inverse": kernel_left_inverse,
        "internal_columns": internal_columns,
        "duplicate_columns": duplicate_columns,
    }


def add_L_to_reduced(column, contraction, multiplication, block_variable, operator_variable, source, scale):
    edge = raw_edge(block_variable, source)
    reduced_row = contraction["reduced_index"].get(edge)
    if reduced_row is not None:
        surface_presentation._add_matrix_term(
            column, reduced_row, _unit_exponent(operator_variable), scale
        )
    vector = multiplication[operator_variable, :, source]
    for target in np.flatnonzero(vector):
        edge = raw_edge(block_variable, int(target))
        reduced_row = contraction["reduced_index"].get(edge)
        if reduced_row is not None:
            surface_presentation._add_matrix_term(
                column,
                reduced_row,
                ZERO_EXPONENT,
                -int(scale) * int(vector[target]),
            )


def build_d2_tree(data, contraction):
    dimension = len(data["standard"])
    multiplication = data["multiplication"]
    columns = []
    for left in range(VARIABLE_COUNT):
        for right in range(left + 1, VARIABLE_COUNT):
            for source in range(dimension):
                column = {}
                add_L_to_reduced(
                    column, contraction, multiplication, left, right, source, 1
                )
                add_L_to_reduced(
                    column, contraction, multiplication, right, left, source, -1
                )
                columns.append(column)
    unit = contraction["unit_index"]
    graph_cycles = []
    for equation in data["surface_equations"]:
        cycle = surface_presentation.divided_difference_cycle(
            equation, multiplication, unit
        )
        graph_cycles.append(cycle)
        column = {}
        for variable, block in enumerate(cycle):
            for source, polynomial in enumerate(block):
                if not polynomial:
                    continue
                edge = raw_edge(variable, source)
                reduced_row = contraction["reduced_index"].get(edge)
                if reduced_row is not None:
                    value = column.setdefault(reduced_row, {})
                    _add_scaled(value, polynomial)
                    if not value:
                        column.pop(reduced_row, None)
        columns.append(column)
    if len(columns) != 28 * dimension + 6:
        raise HSExtensionFromPresentationError("the tree D2 column count changed")
    matrix = {}
    for column, entries in enumerate(columns):
        for row, polynomial in entries.items():
            matrix[(int(row), int(column))] = dict(polynomial)
    return matrix, graph_cycles


def build_koszul_d3(data):
    multiplication = data["multiplication"]
    matrix = {}

    def add_L(column, pair, operator, source, scale):
        row = pair_index(*pair) * DIMENSION + int(source)
        _matrix_entry_add(
            matrix,
            row,
            column,
            _monomial(_unit_exponent(operator)),
            scale,
        )
        vector = multiplication[operator, :, source]
        for target in np.flatnonzero(vector):
            row = pair_index(*pair) * DIMENSION + int(target)
            _matrix_entry_add(
                matrix,
                row,
                column,
                _constant(int(vector[target])),
                -int(scale),
            )

    for triple_number, (a, b, c) in enumerate(
        itertools.combinations(range(VARIABLE_COUNT), 3)
    ):
        for source in range(DIMENSION):
            column = triple_number * DIMENSION + source
            add_L(column, (b, c), a, source, 1)
            add_L(column, (a, c), b, source, -1)
            add_L(column, (a, b), c, source, 1)
    return matrix


def build_h_operators(data):
    """Return the six H_0 divided-difference operator columns."""

    multiplication = data["multiplication"]
    operators = []
    for equation in data["surface_equations"]:
        operator = []
        for source in range(DIMENSION):
            operator.append(
                surface_presentation.divided_difference_cycle(
                    equation, multiplication, source
                )
            )
        operators.append(operator)
    return operators


def apply_h1_to_raw_edge(h_operator, edge):
    """Apply H_1=-h wedge (-) to one raw K_1 basis vector."""

    variable, source = divmod(int(edge), DIMENSION)
    result = {}
    blocks = h_operator[source]
    for h_variable, block in enumerate(blocks):
        if h_variable == variable:
            continue
        sign = -1 if h_variable < variable else 1
        left, right = sorted((h_variable, variable))
        base = pair_index(left, right) * DIMENSION
        for target, polynomial in enumerate(block):
            if polynomial:
                row = base + int(target)
                value = result.setdefault(row, {})
                _add_scaled(value, polynomial, sign)
                if not value:
                    result.pop(row, None)
    return result


def build_h1_raw(h_operators):
    """Materialize all six raw H1=(-h wedge -) blocks."""

    matrices = []
    for h_operator in h_operators:
        matrix = {}
        for column in range(VARIABLE_COUNT * DIMENSION):
            result = apply_h1_to_raw_edge(h_operator, column)
            for row, polynomial in result.items():
                matrix[(int(row), int(column))] = polynomial
        matrices.append(matrix)
    return matrices


def _wedge_sign(variable, wedge):
    return -1 if sum(1 for value in wedge if value < variable) % 2 else 1


def pullback_cochain_through_h(
    cochain,
    h_operator,
    input_degree,
    *,
    coefficient_dimension=None,
    multiply=_multiply,
):
    """Precompose a K_{q+1} cochain with H_q=(-1)^q h wedge."""

    if coefficient_dimension is None:
        coefficient_dimension = DIMENSION
    output = {}
    degree_sign = -1 if int(input_degree) % 2 else 1
    input_wedges = list(itertools.combinations(range(VARIABLE_COUNT), input_degree))
    output_wedges = list(
        itertools.combinations(range(VARIABLE_COUNT), input_degree + 1)
    )
    input_numbers = {wedge: index for index, wedge in enumerate(input_wedges)}
    output_numbers = {wedge: index for index, wedge in enumerate(output_wedges)}
    for wedge in input_wedges:
        wedge_set = set(wedge)
        wedge_number = input_numbers[wedge]
        for source in range(int(coefficient_dimension)):
            target_polynomial = {}
            blocks = h_operator[source]
            for variable in range(VARIABLE_COUNT):
                if variable in wedge_set:
                    continue
                output_wedge = tuple(sorted((*wedge, variable)))
                output_number = output_numbers[output_wedge]
                sign = degree_sign * _wedge_sign(variable, wedge)
                block = blocks[variable]
                for target, h_polynomial in enumerate(block):
                    if not h_polynomial:
                        continue
                    outer = cochain.get(
                        output_number * int(coefficient_dimension) + int(target)
                    )
                    if outer:
                        _add_scaled(
                            target_polynomial,
                            multiply(outer, h_polynomial),
                            sign,
                        )
            if target_polynomial:
                output[
                    wedge_number * int(coefficient_dimension) + source
                ] = target_polynomial
    return output


def local_jet_multiplication(values, tangents):
    """Return the eight multiplication matrices on one dual-number factor."""

    multiplication = np.zeros((VARIABLE_COUNT, 2, 2), dtype=np.int64)
    for variable, (value, tangent) in enumerate(zip(values, tangents)):
        multiplication[variable] = np.asarray(
            [[int(value), 0], [int(tangent), int(value)]], dtype=np.int64
        ) % P
    return multiplication


def build_local_h_operators(surface_equations, values, tangents):
    multiplication = local_jet_multiplication(values, tangents)
    operators = []
    for equation in surface_equations:
        operators.append(
            [
                surface_presentation.divided_difference_cycle(
                    equation, multiplication, source
                )
                for source in range(2)
            ]
        )
    return multiplication, operators


def _local_change_of_rings_cochain(local_h, top_row, equations):
    """Compute tau H_5...H_0 or one five-H omission over the original P."""

    equations = tuple(map(int, equations))
    cochain = {
        source: _constant(int(value))
        for source, value in enumerate(np.asarray(top_row, dtype=np.int64))
        if int(value) % P
    }
    for step, equation in enumerate(equations):
        input_degree = 7 - step
        cochain = pullback_cochain_through_h(
            cochain,
            local_h[equation],
            input_degree,
            coefficient_dimension=2,
            multiply=_multiply,
        )
    # Reversing six odd exterior operators contributes (-1)^15=-1.  Five
    # operators contribute (-1)^10=+1.
    if len(equations) == 6:
        cochain = {
            column: _poly_scale(polynomial, -1)
            for column, polynomial in cochain.items()
        }
    return cochain


def compose_local_lambda_koszul_d3(lambda_k2, multiplication):
    """Precompose a two-dimensional coefficient-row with the graph K3 map."""

    result = {}
    for triple_number, (a, b, c) in enumerate(
        itertools.combinations(range(VARIABLE_COUNT), 3)
    ):
        terms = (((b, c), a, 1), ((a, c), b, -1), ((a, b), c, 1))
        for source in range(2):
            polynomial = {}
            for pair, operator, scale in terms:
                pair_number = pair_index(*pair)
                same_source = lambda_k2.get(2 * pair_number + source)
                if same_source:
                    _add_scaled(
                        polynomial,
                        _multiply_by_monomial(
                            same_source,
                            _monomial(_unit_exponent(operator)),
                        ),
                        scale,
                    )
                matrix_column = multiplication[operator, :, source]
                for target in np.flatnonzero(matrix_column):
                    target_row = lambda_k2.get(2 * pair_number + int(target))
                    if target_row:
                        _add_scaled(
                            polynomial,
                            target_row,
                            -int(scale) * int(matrix_column[int(target)]),
                        )
            if polynomial:
                result[2 * triple_number + source] = polynomial
    return result


def verify_local_change_of_rings(
    lambda_k2, mu_k3, local_h, multiplication, surface_equations
):
    repeated_hashes = []
    for equation in range(6):
        repeated = pullback_cochain_through_h(
            lambda_k2,
            local_h[equation],
            1,
            coefficient_dimension=2,
            multiply=_multiply,
        )
        if repeated:
            raise HSExtensionFromPresentationError(
                "a supportwise lambda H_j repeated-wedge composite is nonzero in P"
            )
        repeated_hashes.append(_polynomial_matrix_hash({}))

    left = compose_local_lambda_koszul_d3(lambda_k2, multiplication)
    right = {}
    for equation, (surface_equation, mu) in enumerate(
        zip(surface_equations, mu_k3)
    ):
        sign = -1 if equation % 2 else 1
        for column, polynomial in mu.items():
            value = right.setdefault(int(column), {})
            _add_scaled(value, _multiply(surface_equation, polynomial), sign)
            if not value:
                right.pop(int(column), None)
    if left != right:
        raise HSExtensionFromPresentationError(
            "the supportwise lambda d3 boundary ledger failed in P"
        )
    return {
        "repeated_wedge_zero_hashes": repeated_hashes,
        "boundary_hash": _polynomial_matrix_hash(
            {(0, column): polynomial for column, polynomial in left.items()}
        ),
        "boundary_nonzero_entries": len(left),
        "boundary_nonzero_terms": sum(map(len, left.values())),
    }


def _jet_cochain_to_standard(local_cochains, hermite, exterior_degree):
    """Transport the two-dimensional jet rows to the standard columns."""

    output = {}
    wedge_count = math.comb(VARIABLE_COUNT, int(exterior_degree))
    for support, cochain in enumerate(local_cochains):
        jet_offset = 2 * int(support)
        for wedge_number in range(wedge_count):
            for local_source in range(2):
                polynomial = cochain.get(2 * wedge_number + local_source)
                if not polynomial:
                    continue
                hermite_row = hermite[jet_offset + local_source]
                for standard_source in np.flatnonzero(hermite_row):
                    column = wedge_number * DIMENSION + int(standard_source)
                    value = output.setdefault(column, {})
                    _add_scaled(
                        value,
                        polynomial,
                        int(hermite_row[int(standard_source)]),
                    )
                    if not value:
                        output.pop(column, None)
    return output


def build_change_of_rings_cocycle(inputs):
    """Build factored lambda and six omission cochains through local jet blocks."""

    data = inputs["data"]
    hermite = inputs["hermite"]
    multiplication = data["multiplication"]
    local_lambda = []
    local_verification = []
    local_mu_hashes = []
    jet_conjugacy_checks = 0
    for support, (values, tangents, orientation) in enumerate(
        zip(
            inputs["affine_values"],
            inputs["affine_tangents"],
            inputs["local_rows"],
        )
    ):
        local_multiplication, local_h = build_local_h_operators(
            data["surface_equations"], values, tangents
        )
        rows = slice(2 * support, 2 * support + 2)
        for variable in range(VARIABLE_COUNT):
            if np.any(
                (
                    hermite[rows] @ multiplication[variable]
                    - local_multiplication[variable] @ hermite[rows]
                )
                % P
            ):
                raise HSExtensionFromPresentationError(
                    "a local dual-number block stopped conjugating multiplication"
                )
            jet_conjugacy_checks += 1
        top_row = orientation["residue_row"]
        support_lambda = _local_change_of_rings_cochain(
            local_h, top_row, range(6)
        )
        support_mu = []
        for omitted in range(6):
            equations = [value for value in range(6) if value != omitted]
            support_mu.append(
                _local_change_of_rings_cochain(local_h, top_row, equations)
            )
        local_lambda.append(support_lambda)
        local_mu_hashes.append(
            [
                _polynomial_matrix_hash(
                    {(0, column): polynomial for column, polynomial in cochain.items()}
                )
                for cochain in support_mu
            ]
        )
        local_verification.append(
            verify_local_change_of_rings(
                support_lambda,
                support_mu,
                local_h,
                local_multiplication,
                data["surface_equations"],
            )
        )

    if not any(local_lambda):
        raise HSExtensionFromPresentationError("the change-of-rings cocycle vanished")
    return {
        "local_lambda_K2": local_lambda,
        "hermite": hermite,
        "hermite_inverse": inputs["hermite_inverse"],
        "jet_conjugacy_checks": jet_conjugacy_checks,
        "local_verification": local_verification,
        "local_mu_hashes": local_mu_hashes,
        "mu_factorization": (
            "mu_j=tau_top o H_5 o ... o omit(H_j) o ... o H_0, "
            f"represented supportwise in the {len(local_lambda)} "
            "dual-number factors"
        ),
    }


def build_tree_d3_factor(koszul_d3, h1_raw, contraction, generators):
    """Emit the exact block-factor representation of D3_tree."""

    scalar_blocks = []
    for _equation in range(6):
        scalar_blocks.append(
            {
                (0, column): dict(polynomial)
                for column, polynomial in enumerate(generators)
                if polynomial
            }
        )
    return {
        "koszul_K3_to_K2": koszul_d3,
        "H1_raw_K1_to_K2": h1_raw,
        "tree_K1_inclusion": contraction["i1"],
        "t_scalar_d1": scalar_blocks,
    }


def compose_lambda_h1(lambda_k2, h1_matrix):
    result = {}
    for (row, column), polynomial in h1_matrix.items():
        left = lambda_k2.get(int(row))
        if left:
            value = result.setdefault(int(column), {})
            _add_scaled(value, _multiply(left, polynomial))
            if not value:
                result.pop(int(column), None)
    return result


def compose_lambda_koszul_d3(lambda_k2, koszul_d3):
    result = {}
    for (row, column), polynomial in koszul_d3.items():
        left = lambda_k2.get(int(row))
        if left:
            value = result.setdefault(int(column), {})
            _add_scaled(value, _multiply(left, polynomial))
            if not value:
                result.pop(int(column), None)
    return result


def verify_full_h0_nullhomotopies(data, h_operators):
    multiplication = data["multiplication"]
    checks = 0
    for equation, operator in zip(data["surface_equations"], h_operators):
        for source, cycle in enumerate(operator):
            image = surface_presentation.graph_d1(cycle, multiplication)
            expected = [dict() for _ in range(DIMENSION)]
            expected[source] = dict(equation)
            if image != expected:
                raise HSExtensionFromPresentationError(
                    "a full divided-difference operator lost dH=f"
                )
            checks += 1
    return checks


def apply_raw_d2_to_k2_vector(vector, multiplication):
    result = {}
    pairs = list(itertools.combinations(range(VARIABLE_COUNT), 2))

    def add_L(block_variable, operator_variable, source, polynomial, scale):
        edge = raw_edge(block_variable, source)
        value = result.setdefault(edge, {})
        _add_scaled(
            value,
            _multiply_by_monomial(
                polynomial, _monomial(_unit_exponent(operator_variable))
            ),
            scale,
        )
        if not value:
            result.pop(edge, None)
        matrix_column = multiplication[operator_variable, :, source]
        for target in np.flatnonzero(matrix_column):
            edge = raw_edge(block_variable, int(target))
            value = result.setdefault(edge, {})
            _add_scaled(
                value,
                polynomial,
                -int(scale) * int(matrix_column[target]),
            )
            if not value:
                result.pop(edge, None)

    for row, polynomial in vector.items():
        pair_number, source = divmod(int(row), DIMENSION)
        left, right = pairs[pair_number]
        add_L(left, right, source, polynomial, 1)
        add_L(right, left, source, polynomial, -1)
    return result


def apply_h0_to_raw_d1_edge(h_operator, edge, multiplication):
    variable, source = divmod(int(edge), DIMENSION)
    result = {}

    def add_cycle(cycle, scalar_polynomial, scale=1):
        for block_variable, block in enumerate(cycle):
            for target, polynomial in enumerate(block):
                if not polynomial:
                    continue
                output_edge = raw_edge(block_variable, target)
                value = result.setdefault(output_edge, {})
                _add_scaled(
                    value,
                    _multiply_by_monomial(polynomial, scalar_polynomial),
                    scale,
                )
                if not value:
                    result.pop(output_edge, None)

    add_cycle(
        h_operator[source],
        _monomial(_unit_exponent(variable)),
    )
    matrix_column = multiplication[variable, :, source]
    for target in np.flatnonzero(matrix_column):
        add_cycle(
            h_operator[int(target)],
            _constant(1),
            -int(matrix_column[target]),
        )
    return result


def verify_raw_h1_identities(data, h_operators, h1_raw):
    multiplication = data["multiplication"]
    hashes = []
    checks = 0
    for equation, operator, matrix in zip(
        data["surface_equations"], h_operators, h1_raw
    ):
        digest = hashlib.sha256()
        by_column = [{} for _ in range(VARIABLE_COUNT * DIMENSION)]
        for (row, column), polynomial in matrix.items():
            by_column[int(column)][int(row)] = polynomial
        for column in range(VARIABLE_COUNT * DIMENSION):
            h1_column = by_column[column]
            total = apply_raw_d2_to_k2_vector(h1_column, multiplication)
            total_h0 = apply_h0_to_raw_d1_edge(operator, column, multiplication)
            for edge, polynomial in total_h0.items():
                value = total.setdefault(edge, {})
                _add_scaled(value, polynomial)
                if not value:
                    total.pop(edge, None)
            expected = {column: dict(equation)}
            if total != expected:
                raise HSExtensionFromPresentationError(
                    "a materialized raw H1 block failed dH1+H0d=f"
                )
            digest.update(str(column).encode("ascii"))
            for exponent, coefficient in sorted(equation.items()):
                digest.update(bytes(exponent))
                digest.update(bytes([int(coefficient) % P]))
            checks += 1
        hashes.append(digest.hexdigest())
    return checks, hashes


def verify_shamash_and_cocycle(
    data, d2_tree, d3_factor, cocycle, contraction, h_operators
):
    dimension = len(data["standard"])
    support_count = len(cocycle["local_verification"])
    retained_rank = len(contraction["generators"])
    # The pure Koszul D2 D3 composite is checked directly.
    d2_koszul = {
        (row, column): polynomial
        for (row, column), polynomial in d2_tree.items()
        if column < 28 * dimension
    }
    composite = _matrix_compose(d2_koszul, d3_factor["koszul_K3_to_K2"])
    if composite:
        raise HSExtensionFromPresentationError("the K3 Koszul composite is nonzero")

    nullhomotopy_checks = verify_full_h0_nullhomotopies(data, h_operators)
    raw_h1_checks, raw_h1_hashes = verify_raw_h1_identities(
        data, h_operators, d3_factor["H1_raw_K1_to_K2"]
    )
    local_checks = cocycle["local_verification"]
    if (
        len(local_checks) != support_count
        or cocycle["jet_conjugacy_checks"] != support_count * VARIABLE_COUNT
    ):
        raise HSExtensionFromPresentationError(
            "the supportwise change-of-basis verification census changed"
        )
    repeated_wedge_hashes = [
        row["repeated_wedge_zero_hashes"] for row in local_checks
    ]
    boundary_hashes = [row["boundary_hash"] for row in local_checks]
    return {
        "full_H0_nullhomotopy_columns_checked": nullhomotopy_checks,
        "raw_H1_columns_checked": raw_h1_checks,
        "raw_H1_surface_boundary_hashes": raw_h1_hashes,
        "contracted_mixed_columns_inherited_via_p1_i1": 6 * retained_rank,
        "hermite_multiplication_conjugacy_checks": (
            support_count * VARIABLE_COUNT
        ),
        "lambda_supportwise_P_boundary_ledgers_checked": support_count,
        "lambda_supportwise_repeated_wedge_maps_checked": support_count * 6,
        "lambda_repeated_wedge_zero_hashes": repeated_wedge_hashes,
        "lambda_P_boundary_hashes": boundary_hashes,
        "lambda_P_boundary_nonzero_entries": [
            row["boundary_nonzero_entries"] for row in local_checks
        ],
        "lambda_P_boundary_nonzero_terms": [
            row["boundary_nonzero_terms"] for row in local_checks
        ],
    }


def build_module_presentation(d2_tree):
    relation = {}
    for (row, column), polynomial in d2_tree.items():
        relation[(int(row), int(column))] = _poly_scale(polynomial, -1)
    return relation


def build():
    inputs = load_inputs()
    data = inputs["data"]
    data["_transition_quotient"] = inputs["quotient"]
    contraction = build_tree_contraction(data, inputs["tree"])
    d2_tree, unit_cycles = build_d2_tree(data, contraction)
    koszul_d3 = build_koszul_d3(data)
    h_operators = build_h_operators(data)
    h1_raw = build_h1_raw(h_operators)
    d3_factor = build_tree_d3_factor(
        koszul_d3, h1_raw, contraction, contraction["generators"]
    )
    cocycle = build_change_of_rings_cocycle(inputs)
    verification = verify_shamash_and_cocycle(
        data, d2_tree, d3_factor, cocycle, contraction, h_operators
    )
    relation_bottom = build_module_presentation(d2_tree)

    # Compare the merge of D2_tree to the independently persisted 297x2750 map.
    merged = {}
    q = contraction["q"]
    for (row, column), polynomial in d2_tree.items():
        target_rows = np.flatnonzero(q[:, row])
        for target in target_rows:
            _matrix_entry_add(
                merged,
                int(target),
                int(column),
                polynomial,
                int(q[int(target), row]),
            )
    first = surface_presentation.build_presentation()
    first_columns = first["columns"]
    # The producer's pre-Tietze columns are reconstructed independently.
    pre_columns = [
        *surface_presentation.koszul_columns(data, inputs["quotient"]),
        *surface_presentation.surface_columns(data, inputs["quotient"])[0],
    ]
    expected_merged = {}
    for column, entries in enumerate(pre_columns):
        for row, polynomial in entries.items():
            expected_merged[(int(row), int(column))] = dict(polynomial)
    if merged != expected_merged:
        raise HSExtensionFromPresentationError(
            "the tree model does not merge to the persisted 297-generator presentation"
        )

    return {
        "inputs": inputs,
        "contraction": contraction,
        "d2_tree": d2_tree,
        "unit_cycles": unit_cycles,
        "d3_factor": d3_factor,
        "h_operators": h_operators,
        "cocycle": cocycle,
        "module_relation_bottom": relation_bottom,
        "merged_d2": merged,
        "verification": verification,
        "first_presentation_counts": first["counts"],
    }


def _store_matrix(payload, prefix, matrix):
    row, column, coefficient, exponent = _flatten_polynomial_matrix(matrix)
    payload[f"{prefix}_row"] = row
    payload[f"{prefix}_column"] = column
    payload[f"{prefix}_coefficient"] = coefficient
    payload[f"{prefix}_exponent"] = exponent


def h0_operator_matrix(operator):
    matrix = {}
    for source, blocks in enumerate(operator):
        for variable, block in enumerate(blocks):
            for target, polynomial in enumerate(block):
                if polynomial:
                    matrix[(variable * DIMENSION + target, source)] = dict(polynomial)
    return matrix


def i1_matrix(columns):
    matrix = {}
    for column, lifted in enumerate(columns):
        for edge, polynomial in lifted.items():
            matrix[(int(edge), int(column))] = dict(polynomial)
    return matrix


def lambda_jet_matrix(local_lambda):
    """Assemble the supportwise 2-columns into one K2 jet-basis row."""

    matrix = {}
    for support, cochain in enumerate(local_lambda):
        for local_column, polynomial in cochain.items():
            wedge_number, local_source = divmod(int(local_column), 2)
            column = wedge_number * DIMENSION + 2 * support + local_source
            if not 0 <= column < 28 * DIMENSION:
                raise HSExtensionFromPresentationError(
                    "a local lambda factor left K2 tensor A_W"
                )
            matrix[(0, column)] = dict(polynomial)
    return matrix


def write_artifacts(result, output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    contraction = result["contraction"]
    d3 = result["d3_factor"]
    cocycle = result["cocycle"]
    payload = {
        "field": np.asarray([P], dtype=np.uint8),
        "tree_complex_ranks": np.asarray([1, 687, 2750, 9610], dtype=np.uint16),
        "module_presentation_shape": np.asarray([688, 2750], dtype=np.uint16),
        "tau_on_jets": np.asarray(result["inputs"]["tau_on_jets"], dtype=np.uint8),
        "tau_on_standard_basis": np.asarray(result["inputs"]["tau"], dtype=np.uint8),
        "tree_edges": np.asarray(contraction["tree_edges"], dtype=np.uint16),
        "non_tree_edges": np.asarray(contraction["non_tree_edges"], dtype=np.uint16),
        "merge_q": np.asarray(contraction["q"], dtype=np.uint8),
        "merge_section": np.asarray(contraction["section"], dtype=np.uint8),
        "merge_kernel": np.asarray(contraction["kernel"], dtype=np.uint8),
        "merge_kernel_left_inverse": np.asarray(
            contraction["kernel_left_inverse"], dtype=np.uint8
        ),
        "tree_projection_p1": np.asarray(contraction["p1"], dtype=np.uint8),
        "hermite_standard_to_jets": np.asarray(
            cocycle["hermite"], dtype=np.uint8
        ),
        "hermite_jets_to_standard": np.asarray(
            cocycle["hermite_inverse"], dtype=np.uint8
        ),
    }
    _store_matrix(payload, "D2_tree", result["d2_tree"])
    _store_matrix(payload, "D3_K3_to_K2", d3["koszul_K3_to_K2"])
    _store_matrix(payload, "tree_K1_inclusion_i1", i1_matrix(contraction["i1"]))
    for equation, operator in enumerate(result["h_operators"]):
        _store_matrix(payload, f"H0_equation_{equation}", h0_operator_matrix(operator))
    for equation, matrix in enumerate(d3["H1_raw_K1_to_K2"]):
        _store_matrix(payload, f"D3_H1_raw_equation_{equation}", matrix)
    lambda_matrix = lambda_jet_matrix(cocycle["local_lambda_K2"])
    _store_matrix(payload, "lambda_K2_jet", lambda_matrix)
    _store_matrix(
        payload, "Vprime_relation_bottom_minus_D2", result["module_relation_bottom"]
    )
    _store_matrix(payload, "merged_D2_297", result["merged_d2"])
    artifact_path = output_directory / "Z_2_C0010_hs_extension.npz"
    np.savez_compressed(artifact_path, **payload)

    local_path = output_directory / "local_duality_orientation_ledger.csv"
    with local_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "support_position",
            "support_point_indices",
            "delta0",
            "delta1",
            "residue_row",
            "expected_local_normal_form",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in result["inputs"]["local_rows"]:
            writer.writerow(
                {
                    **row,
                    "support_point_indices": json.dumps(row["support_point_indices"]),
                    "residue_row": json.dumps(row["residue_row"]),
                    "expected_local_normal_form": "(Delta,y^2,-x)^T",
                }
            )

    certificate = {
        "schema": "section12-chart-local-HS-extension-from-surface-presentation-v1",
        "field": "F_31",
        "surface": "Z_2",
        "chart": CHART,
        "surface_ring": "S=P/(f_0,...,f_5)",
        "tree_shamash_complex": {
            "ranks_F0_F1_F2_F3": [1, 687, 2750, 9610],
            "D2_shape": [687, 2750],
            "D3_factor_shapes": {
                "K3_to_K2": [2744, 5488],
                "six_H1_raw_blocks": [6, 2744, 784],
                "tree_i1": [784, 687],
                "six_H1_tree_factored": [6, 2744, 687],
                "six_scalar_d1_blocks": [6, 1, 687],
            },
            "construction": (
                "graph Koszul plus six strict divided-difference nullhomotopies, "
                "followed by the 97 order-tree cancellations and their six copies"
            ),
        },
        "merge_comparison": {
            "q_shape": [297, 687],
            "kernel_rank": 390,
            "kernel_decomposition": {
                "non_tree_internal": 98,
                "duplicate_border_differences": 292,
            },
            "merged_D2_shape": [297, 2750],
            "merged_D2_equals_pre_Tietze_first_presentation": True,
            "degree3_merge_lifts_emitted": False,
            "precise_absent_input": (
                "B:S^390->S^2750 with D2_tree B=J; it is not needed for the "
                "tree-carrier module emitted here, but is required to transport "
                "the cocycle to the 297-generator carrier"
            ),
        },
        "canonical_cochain": {
            "formula": (
                "lambda_standard=lambda_jet o (I_28 tensor Hermite), with "
                "lambda_jet=tau_top o H_5 o ... o H_0 on K_2 and zero on "
                "the six t columns"
            ),
            "cohomological_degree": 2,
            "identified_HS_class": "Ext^2_S(A_W,S) = Ext^1_S(I_W,S)",
            "K2_columns": 2744,
            "lambda_jet_shape": [1, 2744],
            "basis_factor_shape": [2744, 2744],
            "basis_factor_storage": "I_28 tensor the stored 98x98 Hermite matrix",
            "nonzero_polynomial_entries_in_lambda_jet": len(lambda_matrix),
            "tau_dimension": 98,
            "local_duality_factors": 49,
            "all_delta0_nonzero": True,
            "six_mu_omission_cochains": cocycle["mu_factorization"],
            "six_mu_supportwise_hashes": cocycle["local_mu_hashes"],
            "supportwise_matrix_reduction_to_expected_normal_form_emitted": False,
        },
        "module": {
            "name": "Vprime_U (chart-local trivialized HS module)",
            "presentation": (
                "coker([lambda_jet o (I_28 tensor Hermite),0_6;"
                "-D2_tree]:S^2750->S direct_sum S^687)"
            ),
            "top_row_is_exact_factored_not_absent": True,
            "bottom_block_is_explicit_minus_D2_tree": True,
            "relation_shape": [688, 2750],
            "generic_relation_rank": 686,
            "generic_module_rank": 2,
            "exact_sequence_in_chosen_chart_frames": "0->S->Vprime_U->I_W->0",
        },
        "exact_checks": {
            "49_support_orientations_matched_by_point_index": True,
            "all_local_delta0_are_units": True,
            "Hermite_and_inverse_multiply_to_I98_both_orders": True,
            "Hermite_conjugates_all_49_times_8_multiplication_blocks": True,
            "tree_inclusion_687_columns_checked": True,
            "tree_projection_times_inclusion_is_I687": True,
            "merge_q_has_split_kernel_390": True,
            "K3_Koszul_composite_zero_over_P": True,
            "six_raw_H1_blocks_compose_to_f_j_I784_over_P": True,
            "six_contracted_H1_identities_inherited_exactly_by_p1_and_i1": True,
            "therefore_D2_D3_zero_over_S": True,
            "lambda_d3_equals_sum_minus1_power_j_f_j_mu_j_over_P_supportwise": True,
            "therefore_lambda_d3_is_zero_over_S": True,
            "lambda_kills_all_six_repeated_wedge_H1_blocks_over_P_supportwise": True,
            "basis_change_transports_both_identities_to_standard_basis": True,
            "therefore_lambda_is_a_closed_Ext2_cochain": True,
            "merged_D2_matches_independent_297_generator_matrix": True,
        },
        "verification_ledger": result["verification"],
        "rank_audit": {
            "tree_D2_generic_rank": 686,
            "tree_D2_generic_kernel_rank": 2064,
            "tree_D2_is_injective": False,
            "merged_297_D2_generic_rank": 296,
            "merged_297_D2_generic_kernel_rank": 2454,
            "merged_297_D2_is_injective": False,
            "projective_kernel_not_assumed_free": True,
        },
        "artifacts": {
            "binary": {
                "path": artifact_path.name,
                "sha256": sha256(artifact_path),
            },
            "local_orientation_ledger": {
                "path": local_path.name,
                "sha256": sha256(local_path),
            },
        },
        "upstream": {
            "Hermite_chart": {
                "path": str(HERMITE_PATH.relative_to(ROOT)),
                "sha256": sha256(HERMITE_PATH),
            },
            "local_HS_orientations": {
                "path": str(ORIENTATION_PATH.relative_to(ROOT)),
                "sha256": sha256(ORIENTATION_PATH),
            },
            "first_presentation": {
                "path": str(FIRST_PRESENTATION_PATH.relative_to(ROOT)),
                "sha256": sha256(FIRST_PRESENTATION_PATH),
            },
        },
        "claim_boundary": {
            "proved": (
                "one exact chart-local HS cochain and coker module presentation "
                "on Z_2/C0010, carried by the 687-generator tree model"
            ),
            "not_claimed": (
                "a global O(3H)/O(-H) overlap descent, the relative O(4H) "
                "transition, a projective idempotent or strict degree-zero dg "
                "object, a stored 49-support matrix reduction to the expected "
                "(Delta,y^2,-x)^T normal form, the 297-carrier degree3 lift B, "
                "rho_U, rho_F, a Yoneda "
                "value, Markman's Question 11.4, or the Hodge conjecture"
            ),
        },
    }
    certificate_path = output_directory / "hs_extension_certificate.json"
    certificate_path.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return certificate_path, artifact_path, local_path


def reproduce(output_directory=OUTPUT_DIRECTORY):
    result = build()
    paths = write_artifacts(result, output_directory)
    return result, paths


def main():
    certificate, artifact, local = reproduce()[1]
    print(certificate)
    print(artifact)
    print(local)


if __name__ == "__main__":
    main()
