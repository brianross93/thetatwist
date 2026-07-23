"""Derive the supportwise residue/counit transport for every standard K edge.

The earlier ``e0`` calculation happened to compress the top-residue transport
to one scalar Laurent monomial.  That is not the invariant construction.  The
invariant construction is local Frobenius transport on each retained
dual-number factor:

    tau_source J8 = tau_target m_r.

Since the constant term of ``tau_target`` is a unit, multiplication by the
unknown dual number ``r`` is recovered uniquely.  The resulting value/tangent
pairs are then reconstructed as one exact element of the whole overlap
algebra.  This retains the two p/q coordinates; it never evaluates only on the
diagonal.

The same transport factors uniquely as

    r = det(B) * w_4H * n,

where ``w_4H`` is the fixed normalized O(4H) frame and ``n`` is the
supportwise normalization of the target residue/counit.  For ``e0``, ``n`` is
the old constant 25.  For the other generators it is a genuine overlap unit.
"""

from __future__ import annotations

import hashlib
import json
import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_all_chart_presentations as all_charts,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_rank2_k_generator_overlap_j1 as k_j1,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_rank2_k_generator_shamash as k_shamash,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_rank2_k_generator_tree_contraction_driver as driver,
)


P = 31
VARIABLE_COUNT = 8
EQUATION_COUNT = 6
SURFACE = 2
CHART = "C0010"
ZERO = (0,) * VARIABLE_COUNT


class HSKGeneratorResidueTransportError(RuntimeError):
    """A local Frobenius law or whole-overlap reconstruction failed."""


def array_sha256(value):
    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(value)).tobytes()
    ).hexdigest()


def pair_ledger_sha256(value):
    """Hash dual-number pair ledgers in one platform-independent encoding."""

    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(value, dtype="<i8")).tobytes()
    ).hexdigest()


def dual_matrix(pair):
    """Return multiplication by ``value + tangent*epsilon``."""

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
        raise HSKGeneratorResidueTransportError(
            "a required overlap normalization is not a unit"
        )
    inverse_value = pow(value, -1, P)
    return (
        inverse_value,
        (-tangent * inverse_value * inverse_value) % P,
    )


def dual_divide(numerator, denominator):
    return dual_multiply(numerator, dual_inverse(denominator))


def scalar_laurent_pair(scalar, exponent, values, tangents):
    value, tangent = k_shamash.hermite.monomial_jet(
        tuple(map(int, exponent)), values, tangents
    )
    return int(scalar) * int(value) % P, int(scalar) * int(tangent) % P


def pair_from_scalar_operator(operator):
    """Recover a dual number and certify that the operator is scalar."""

    operator = np.asarray(operator, dtype=np.int64) % P
    pair = int(operator[0, 0]), int(operator[1, 0])
    if not np.array_equal(operator, dual_matrix(pair)):
        raise HSKGeneratorResidueTransportError(
            "a top-wedge map is not multiplication by one dual number"
        )
    return pair


def frobenius_multiplier(target_residue, transported_residue):
    """Solve ``target_residue*m_r = transported_residue`` uniquely."""

    delta1, delta0 = map(lambda item: int(item) % P, target_residue)
    left0, left1 = map(lambda item: int(item) % P, transported_residue)
    if not delta0:
        raise HSKGeneratorResidueTransportError(
            "the target Frobenius constant term is not invertible"
        )
    value = left1 * pow(delta0, -1, P) % P
    tangent = (left0 - delta1 * value) * pow(delta0, -1, P) % P
    pair = value, tangent
    replay = (
        np.asarray(target_residue, dtype=np.int64).reshape(1, 2)
        @ dual_matrix(pair)
    ) % P
    if not np.array_equal(
        replay, np.asarray(transported_residue, dtype=np.int64).reshape(1, 2) % P
    ):
        raise HSKGeneratorResidueTransportError(
            "the recovered Frobenius multiplier did not replay"
        )
    return pair


@contextmanager
def bind_generator(generator):
    """Bind the selected edge only while its overlap data are constructed."""

    if generator not in driver.GENERATOR_EXPONENTS:
        raise HSKGeneratorResidueTransportError(
            f"unknown standard K generator {generator!r}"
        )
    previous_overlap = k_j1.GENERATOR_EXPONENTS
    previous_shamash = k_shamash.GENERATOR_EXPONENTS
    exponents = driver.GENERATOR_EXPONENTS[generator]
    k_j1.GENERATOR_EXPONENTS = exponents
    k_shamash.GENERATOR_EXPONENTS = exponents
    try:
        yield
    finally:
        k_j1.GENERATOR_EXPONENTS = previous_overlap
        k_shamash.GENERATOR_EXPONENTS = previous_shamash


def build_inputs(generator):
    """Build the common frame and regular-sequence data without e0 constants."""

    with bind_generator(generator):
        graph_input = k_j1.build_data()
    context = k_shamash.iw_solver.load_problem_context()
    canonical_equations, equation_labels = k_shamash.shamash.surface_equations(
        context, SURFACE, CHART
    )
    descent, _extension, _stable, _elements, _refined, _edge, theta, labels = (
        k_shamash.k_chain.load_inputs()
    )
    records = k_shamash.k_chain.group_records(descent, theta)
    exponents = driver.GENERATOR_EXPONENTS[generator]
    record = next(
        item for item in records if tuple(item["exponents"]) == exponents
    )
    if labels != [2, 3]:
        raise HSKGeneratorResidueTransportError(
            "the Z_2 theta labels changed"
        )
    for left, right in zip(
        record["matrices"], graph_input["record"]["factor_matrices"]
    ):
        if np.any((np.asarray(left) - np.asarray(right)) % P):
            raise HSKGeneratorResidueTransportError(
                "the overlap and chain K ledgers disagree"
            )

    bits = tuple(map(int, context["chart_bits"][CHART]))
    factor_data = []
    original_substitutions = [None] * VARIABLE_COUNT
    target_substitutions = [None] * VARIABLE_COUNT
    for factor, matrix in enumerate(record["matrices"]):
        frame, numerator_u, numerator_v = k_j1.affine_forms(
            matrix, bits[factor]
        )
        data = k_shamash.adapted_factor(
            frame, numerator_u, numerator_v, factor
        )
        factor_data.append(data)
        for local in range(2):
            original_substitutions[2 * factor + local] = data[
                "original_substitutions"
            ][local]
            target_substitutions[2 * factor + local] = data[
                "target_substitutions"
            ][local]

    source_equations = [
        k_shamash.substitute_scalar_polynomial(
            polynomial, original_substitutions
        )
        for polynomial in canonical_equations
    ]
    orientation_scalar = 1
    for data in factor_data:
        orientation_scalar = (
            orientation_scalar
            * k_shamash.determinant_mod(data["transform"])
        ) % P
    b_entries = k_shamash.equation_change_data(
        factor_data,
        record,
        source_equations,
        canonical_equations,
        target_substitutions,
    )
    inverse_variables = tuple(
        data["inverse_variable"]
        for data in factor_data
        if data["inverse_variable"] is not None
    )
    weight_exponent = tuple(
        4 if variable in inverse_variables else 0
        for variable in range(VARIABLE_COUNT)
    )
    weight_scalar = pow(int(record["H_normalization"]), 4, P)
    return {
        "generator": generator,
        "exponents": exponents,
        "graph_input": graph_input,
        "record": record,
        "factor_data": factor_data,
        "canonical_equations": canonical_equations,
        "equation_labels": equation_labels,
        "source_equations": source_equations,
        "target_substitutions": target_substitutions,
        "orientation_scalar": orientation_scalar,
        "b_entries": b_entries,
        "inverse_variables": inverse_variables,
        "weight_exponent": weight_exponent,
        "weight_scalar": weight_scalar,
    }


def assert_lambda_rows_equal(left, right, message):
    wedges = set(left) | set(right)
    for wedge in wedges:
        left_value = np.asarray(left.get(wedge, np.zeros((1, 2))), dtype=np.int64) % P
        right_value = np.asarray(
            right.get(wedge, np.zeros((1, 2))), dtype=np.int64
        ) % P
        if not np.array_equal(left_value, right_value):
            raise HSKGeneratorResidueTransportError(message)


def direct_lambda_matrices(h0_rows, residue, multiplication, inverses):
    row = k_shamash.evaluated_lambda_row(
        h0_rows, residue, multiplication, inverses
    )
    return {
        wedge: k_shamash.evaluate_row_multi(
            polynomial, multiplication, inverses
        )
        for wedge, polynomial in row.items()
    }


def formal_lambda_matrices(
    equations, values, tangents, residue, multiplication, inverses
):
    row = k_shamash.lambda_transport._lambda_row(
        equations, values, tangents, residue
    )
    return {
        wedge: k_shamash.evaluate_row_multi(
            polynomial, multiplication, inverses
        )
        for wedge, polynomial in row.items()
    }


def local_transport_support(
    adapted_values,
    adapted_tangents,
    target_values,
    target_tangents,
    source_equations,
    target_equations,
    factor_data,
):
    """Build only the top graph and H0 rows consumed by this derivation."""

    edge_multiplication = k_shamash.shamash.local_multiplication(
        adapted_values, adapted_tangents
    )
    target_multiplication = k_shamash.shamash.local_multiplication(
        target_values, target_tangents
    )
    inverses = {
        data["inverse_variable"]: k_shamash.shamash.inverse_2x2(
            edge_multiplication[data["inverse_variable"]]
        )
        for data in factor_data
        if data["inverse_variable"] is not None
    }
    graph_maps = k_shamash.graph_compounds_general(
        target_multiplication, factor_data
    )
    edge_h0 = []
    target_h0_raw = []
    for source_equation, target_equation in zip(
        source_equations, target_equations
    ):
        if np.any(
            k_shamash.shamash.evaluate_matrix_polynomial(
                source_equation, edge_multiplication
            )
        ):
            raise HSKGeneratorResidueTransportError(
                "an adapted source equation is nonzero"
            )
        if np.any(
            k_shamash.shamash.evaluate_matrix_polynomial(
                target_equation, target_multiplication
            )
        ):
            raise HSKGeneratorResidueTransportError(
                "a canonical target equation is nonzero"
            )
        edge_h0.append(
            k_shamash.shamash.divided_difference_matrix(
                source_equation, edge_multiplication
            )
        )
        target_h0_raw.append(
            k_shamash.shamash.divided_difference_matrix(
                target_equation, target_multiplication
            )
        )
    return {
        "edge_multiplication": edge_multiplication,
        "target_multiplication": target_multiplication,
        "inverses": inverses,
        "graph_maps": graph_maps,
        "edge_h0": edge_h0,
        "target_h0_raw": target_h0_raw,
    }


def reconstruct_overlap_element(
    pairs, affine, basis, hermite, inverse, multiplication
):
    """Reconstruct one quotient element and retain its value/tangent blocks."""

    if len(pairs) != len(affine):
        raise HSKGeneratorResidueTransportError(
            "the support pair and affine-factor counts differ"
        )
    jets = np.asarray(pairs, dtype=np.int64).reshape(-1) % P
    coefficients = inverse @ jets % P
    operator = k_shamash.hermite.coordinate_multiplication_from_vector(
        coefficients, basis, multiplication
    )
    operator_on_jets = hermite @ operator @ inverse % P
    expected = np.zeros_like(operator_on_jets)
    for position, pair in enumerate(pairs):
        expected[2 * position : 2 * position + 2, 2 * position : 2 * position + 2] = (
            dual_matrix(pair)
        )
    if not np.array_equal(operator_on_jets, expected):
        raise HSKGeneratorResidueTransportError(
            "the p/q-separated overlap reconstruction lost a value or tangent"
        )
    # Independent replay from the reconstructed standard coefficients.
    if not np.array_equal(hermite @ coefficients % P, jets):
        raise HSKGeneratorResidueTransportError(
            "the reconstructed overlap element did not replay its Hermite jets"
        )
    return {
        "jets": jets,
        "coefficients": coefficients,
        "operator": operator,
        "operator_on_jets": operator_on_jets,
        "unit": all(int(pair[0]) % P for pair in pairs),
    }


def build(generator):
    """Derive and certify the residue/counit normalization for one generator."""

    inputs = build_inputs(generator)
    orientation = all_charts.orientation_map(SURFACE)
    support_rows = []
    adapted_affine = []
    residue_pairs = []
    orientation_coboundary_pairs = []
    normalization_pairs = []
    raw_weight_pairs = []
    b_determinant_pairs = []
    weight_pairs = []
    top_graph_pairs = []
    source_residues = []
    target_residues = []
    lambda_checks = 0

    for support_position, point_row in enumerate(
        inputs["graph_input"]["point_rows"]
    ):
        source_support = tuple(json.loads(point_row["source_support"]))
        target_support = tuple(json.loads(point_row["target_support"]))
        source_values = tuple(json.loads(point_row["source_affine_values"]))
        source_tangents = tuple(
            json.loads(point_row["source_affine_tangent"])
        )
        target_values = tuple(json.loads(point_row["target_affine_values"]))
        target_tangents = tuple(
            json.loads(point_row["target_affine_tangent"])
        )
        adapted_values, adapted_tangents = k_shamash.adapt_jet(
            source_values, source_tangents, inputs["factor_data"]
        )
        adapted_affine.append((adapted_values, adapted_tangents))
        comparison = local_transport_support(
            adapted_values,
            adapted_tangents,
            target_values,
            target_tangents,
            inputs["source_equations"],
            inputs["canonical_equations"],
            inputs["factor_data"],
        )

        source_delta0, source_delta1 = orientation[source_support]
        target_delta0, target_delta1 = orientation[target_support]
        source_residue = (
            inputs["orientation_scalar"] * int(source_delta1) % P,
            inputs["orientation_scalar"] * int(source_delta0) % P,
        )
        target_residue = (
            int(target_delta1) % P,
            int(target_delta0) % P,
        )
        source_residues.append(source_residue)
        target_residues.append(target_residue)

        top_wedge = tuple(range(VARIABLE_COUNT))
        top_graph = comparison["graph_maps"][8][(top_wedge, top_wedge)]
        top_graph_operator = k_shamash.evaluate_block_multi(
            top_graph,
            comparison["edge_multiplication"],
            comparison["inverses"],
        )
        top_graph_pair = pair_from_scalar_operator(top_graph_operator)
        transported_residue = (
            np.asarray(source_residue, dtype=np.int64).reshape(1, 2)
            @ top_graph_operator
        )[0] % P
        residue_pair = frobenius_multiplier(
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
            source_orientation_pair,
            target_orientation_pair,
        )
        formula_residue_pair = dual_multiply(
            dual_multiply(
                (inputs["orientation_scalar"], 0),
                orientation_coboundary_pair,
            ),
            top_graph_pair,
        )
        if formula_residue_pair != residue_pair:
            raise HSKGeneratorResidueTransportError(
                "the orientation-coboundary formula for residue transport failed"
            )

        b_determinant_pair = (1, 0)
        for entry in inputs["b_entries"]:
            entry_pair = scalar_laurent_pair(
                entry["scalar"],
                entry["exponent"],
                adapted_values,
                adapted_tangents,
            )
            b_determinant_pair = dual_multiply(
                b_determinant_pair, entry_pair
            )
        weight_pair = scalar_laurent_pair(
            inputs["weight_scalar"],
            inputs["weight_exponent"],
            adapted_values,
            adapted_tangents,
        )
        raw_weight_pair = dual_divide(
            residue_pair,
            b_determinant_pair,
        )
        normalization_pair = dual_divide(
            raw_weight_pair,
            weight_pair,
        )
        if dual_multiply(
            dual_multiply(b_determinant_pair, weight_pair),
            normalization_pair,
        ) != residue_pair:
            raise HSKGeneratorResidueTransportError(
                "the residue/B/O(4H)/normalization factorization failed"
            )

        source_direct = direct_lambda_matrices(
            comparison["edge_h0"],
            source_residue,
            comparison["edge_multiplication"],
            comparison["inverses"],
        )
        source_formal = formal_lambda_matrices(
            inputs["source_equations"],
            adapted_values,
            adapted_tangents,
            source_residue,
            comparison["edge_multiplication"],
            comparison["inverses"],
        )
        assert_lambda_rows_equal(
            source_direct,
            source_formal,
            "the formal source lambda row failed exact local replay",
        )
        target_raw = direct_lambda_matrices(
            comparison["target_h0_raw"],
            target_residue,
            comparison["target_multiplication"],
            {},
        )
        target_formal = formal_lambda_matrices(
            inputs["canonical_equations"],
            target_values,
            target_tangents,
            target_residue,
            comparison["target_multiplication"],
            {},
        )
        assert_lambda_rows_equal(
            target_raw,
            target_formal,
            "the formal target lambda row failed exact local replay",
        )
        normalized_residue = (
            np.asarray(target_residue, dtype=np.int64).reshape(1, 2)
            @ dual_matrix(normalization_pair)
        )[0] % P
        target_normalized = direct_lambda_matrices(
            comparison["target_h0_raw"],
            normalized_residue,
            comparison["target_multiplication"],
            {},
        )
        expected_normalized = {
            wedge: matrix @ dual_matrix(normalization_pair) % P
            for wedge, matrix in target_raw.items()
        }
        assert_lambda_rows_equal(
            target_normalized,
            expected_normalized,
            "the supportwise target normalization failed lambda replay",
        )
        lambda_checks += 3

        residue_pairs.append(residue_pair)
        orientation_coboundary_pairs.append(
            orientation_coboundary_pair
        )
        normalization_pairs.append(normalization_pair)
        raw_weight_pairs.append(raw_weight_pair)
        b_determinant_pairs.append(b_determinant_pair)
        weight_pairs.append(weight_pair)
        top_graph_pairs.append(top_graph_pair)
        support_rows.append(
            {
                "support_position": support_position,
                "source_support": list(source_support),
                "target_support": list(target_support),
                "top_graph_pair": list(top_graph_pair),
                "orientation_coboundary_pair": list(
                    orientation_coboundary_pair
                ),
                "residue_transport_pair": list(residue_pair),
                "B_determinant_pair": list(b_determinant_pair),
                "raw_counit_weight_pair": list(raw_weight_pair),
                "O4H_weight_pair": list(weight_pair),
                "target_normalization_pair": list(normalization_pair),
                "lambda_replay_checks": 3,
            }
        )

    basis, hermite, tested, final_degree = (
        k_shamash.hermite.select_standard_monomials(adapted_affine)
    )
    multiplication, _border = k_shamash.hermite.border_multiplication(
        basis, hermite, adapted_affine
    )
    inverse = k_shamash.hermite.inverse_mod(hermite)
    reconstructed = {
        "top_graph": reconstruct_overlap_element(
            top_graph_pairs,
            adapted_affine,
            basis,
            hermite,
            inverse,
            multiplication,
        ),
        "residue_transport": reconstruct_overlap_element(
            residue_pairs,
            adapted_affine,
            basis,
            hermite,
            inverse,
            multiplication,
        ),
        "orientation_coboundary": reconstruct_overlap_element(
            orientation_coboundary_pairs,
            adapted_affine,
            basis,
            hermite,
            inverse,
            multiplication,
        ),
        "B_determinant": reconstruct_overlap_element(
            b_determinant_pairs,
            adapted_affine,
            basis,
            hermite,
            inverse,
            multiplication,
        ),
        "raw_counit_weight": reconstruct_overlap_element(
            raw_weight_pairs,
            adapted_affine,
            basis,
            hermite,
            inverse,
            multiplication,
        ),
        "O4H_weight": reconstruct_overlap_element(
            weight_pairs,
            adapted_affine,
            basis,
            hermite,
            inverse,
            multiplication,
        ),
        "target_normalization": reconstruct_overlap_element(
            normalization_pairs,
            adapted_affine,
            basis,
            hermite,
            inverse,
            multiplication,
        ),
    }
    for name in (
        "residue_transport",
        "orientation_coboundary",
        "B_determinant",
        "raw_counit_weight",
        "O4H_weight",
        "target_normalization",
    ):
        if not reconstructed[name]["unit"]:
            raise HSKGeneratorResidueTransportError(
                f"{name} is not a unit on the whole overlap"
            )

    source_tau_local = np.asarray(source_residues, dtype=np.int64).reshape(1, -1) % P
    target_tau_local = np.asarray(target_residues, dtype=np.int64).reshape(1, -1) % P
    source_tau_standard = source_tau_local @ hermite % P
    target_tau_standard = target_tau_local @ hermite % P
    if not np.array_equal(
        source_tau_standard @ reconstructed["top_graph"]["operator"] % P,
        target_tau_standard
        @ reconstructed["residue_transport"]["operator"]
        % P,
    ):
        raise HSKGeneratorResidueTransportError(
            "the whole-overlap top-residue identity failed after p/q separation"
        )
    if not np.array_equal(
        reconstructed["residue_transport"]["operator"],
        reconstructed["B_determinant"]["operator"]
        @ reconstructed["raw_counit_weight"]["operator"]
        % P,
    ):
        raise HSKGeneratorResidueTransportError(
            "the whole-overlap raw counit factorization failed"
        )
    if not np.array_equal(
        reconstructed["raw_counit_weight"]["operator"],
        reconstructed["O4H_weight"]["operator"]
        @ reconstructed["target_normalization"]["operator"]
        % P,
    ):
        raise HSKGeneratorResidueTransportError(
            "the whole-overlap normalization factorization failed"
        )

    normalization_constant = (
        len(set(normalization_pairs)) == 1
        and int(normalization_pairs[0][1]) % P == 0
    )
    return {
        "generator": generator,
        "K_generator_exponents": list(inputs["exponents"]),
        "support_rows": support_rows,
        "pairs": {
            "top_graph": top_graph_pairs,
            "residue_transport": residue_pairs,
            "orientation_coboundary": orientation_coboundary_pairs,
            "B_determinant": b_determinant_pairs,
            "raw_counit_weight": raw_weight_pairs,
            "O4H_weight": weight_pairs,
            "target_normalization": normalization_pairs,
        },
        "reconstructed": reconstructed,
        "certificate": {
            "schema": "section12-HS-rank2-K-residue-transport-v1",
            "generator": generator,
            "K_generator_exponents": list(inputs["exponents"]),
            "overlap_factors": len(support_rows),
            "overlap_length": 2 * len(support_rows),
            "adapted_orientation_scalar": inputs["orientation_scalar"],
            "inverse_variables": list(inputs["inverse_variables"]),
            "O4H_weight_scalar": inputs["weight_scalar"],
            "O4H_weight_exponent": list(inputs["weight_exponent"]),
            "local_top_residue_replays": len(support_rows),
            "local_lambda_replays": lambda_checks,
            "p_q_value_tangent_reconstructions": len(reconstructed),
            "whole_overlap_top_residue_identity": True,
            "orientation_coboundary_formula_exact": True,
            "whole_overlap_normalization_factorization": True,
            "target_normalization_is_unit": bool(
                reconstructed["target_normalization"]["unit"]
            ),
            "raw_counit_weight_is_unit": bool(
                reconstructed["raw_counit_weight"]["unit"]
            ),
            "target_normalization_is_constant": normalization_constant,
            "target_normalization_constant_pair": (
                list(normalization_pairs[0])
                if normalization_constant
                else None
            ),
            "target_normalization_pairs_sha256": pair_ledger_sha256(
                normalization_pairs
            ),
            "raw_counit_weight_pairs_sha256": pair_ledger_sha256(
                raw_weight_pairs
            ),
            "residue_transport_pairs_sha256": pair_ledger_sha256(
                residue_pairs
            ),
            "orientation_coboundary_pairs_sha256": pair_ledger_sha256(
                orientation_coboundary_pairs
            ),
            "pair_ledger_hash_encoding": "signed-little-endian-int64-C-order",
            "orientation_coboundary_distinct_pairs": len(
                set(orientation_coboundary_pairs)
            ),
            "overlap_basis_size": len(basis),
            "overlap_monomials_tested": tested,
            "overlap_final_degree": final_degree,
            "scope": (
                "exact supportwise residue/counit transport and p/q-separated "
                "whole-overlap reconstruction only; no tree contraction, K "
                "coherence, cross-RHom/Yoneda, Question 11.4, or Hodge claim"
            ),
            "normalization_convention": (
                "raw q=r/det(B) is invariant; n=q/w_H uses the fixed normalized "
                "O(4H) frame. The e0 value n=25 is recovered, not imposed."
            ),
        },
    }


def compact_certificate(generator):
    return build(generator)["certificate"]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--generator",
        required=True,
        choices=tuple(driver.GENERATOR_EXPONENTS),
    )
    args = parser.parse_args()
    print(json.dumps(compact_certificate(args.generator), indent=2, sort_keys=True))
