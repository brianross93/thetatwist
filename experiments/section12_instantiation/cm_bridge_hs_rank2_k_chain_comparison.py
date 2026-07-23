"""Bind the exact K-orbit data for the Z_2 rank-two presentation.

The repository has one whole-ring presentation on ``Z_2/C0010`` and an
exact projective action of ``K=(Z/3)^4`` on the four plane-cubic factors.
This producer answers a narrower question than the eventual descent solve:
which parts of a K-linearized comparison are already determined by those
inputs, and which matrices still have to be computed?

For every one of the 81 group elements it emits the rational affine ring
map on the translated copy of ``C0010``.  It also verifies the two theta
equations defining ``Z_2`` are eigenvectors, records their characters, and
constructs the normalized ``O(4H)`` frame circuit.  The normalization is
the one fixed by the theta eigenline in ``descent_section``.  With

    ell_k = product_i (k^* F_i / F_i),
    n_k   = inverse(raw theta_0 character),

the normalized H factor is ``n_k*ell_k`` and the relative extension factor
is ``(n_k*ell_k)^4``.  All 81^2 group products are checked exactly.

The existing files do *not* contain a second, independently normalized
whole-ring tree presentation on a nonidentity translated chart, nor the
tree/Shamash comparison matrices ``J1,J2``.  Therefore emitting
``L0=I,L1=I`` after defining ``A_k=phi_k(A_e)`` would be a tautology, not a
K-linearization.  This producer does not do that.  It emits the exact
matrix contract whose solution would certify

    L0_k A_e = A_k L1_k

and the 324 typed Cayley-edge jobs on which group coherence must be checked.
No supportwise or finite-quotient identity is substituted for that
whole-ring comparison.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
P = 31
BASE_CHART = "C0010"
BASE_BITS = (0, 0, 1, 0)
SURFACE = "Z_2"
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_rank2_k_chain_comparison"
)

DESCENT = ROOT / (
    "results/section12_instantiation/descent_section/"
    "descent_section_certificate.json"
)
K_ELEMENTS = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit/"
    "K_degree2_action_elements.csv"
)
K_STABLE_CHARTS = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit/"
    "K_stable_product_chart_refinement.csv"
)
K_REFINED_INTERSECTIONS = ROOT / (
    "results/section12_instantiation/cm_bridge_k_equivariant_edge_audit/"
    "k_refined_chart_intersections.csv"
)
K_REFINED_EDGE_PULLBACKS = ROOT / (
    "results/section12_instantiation/cm_bridge_k_equivariant_edge_audit/"
    "k_refined_edge_pullbacks.csv"
)
THETA_COEFFICIENTS = ROOT / (
    "results/section12_instantiation/independent_incidence_fiber/"
    "candidate_section_coefficients.csv"
)
LOCAL_SEEDS = ROOT / (
    "results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/"
    "six_hs_local_seeds.csv"
)
EXTENSION_CERTIFICATE = ROOT / (
    "results/section12_instantiation/"
    "cm_bridge_hs_extension_from_surface_presentation/"
    "hs_extension_certificate.json"
)
EXTENSION_ARTIFACT = ROOT / (
    "results/section12_instantiation/"
    "cm_bridge_hs_extension_from_surface_presentation/"
    "Z_2_C0010_hs_extension.npz"
)


class HSRank2KChainComparisonError(RuntimeError):
    """An exact orbit input or a group/coherence identity changed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError("cannot infer a CSV schema from no rows")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def matrix_multiply(left, right):
    return [
        [
            sum(
                int(left[row][middle]) * int(right[middle][column])
                for middle in range(3)
            )
            % P
            for column in range(3)
        ]
        for row in range(3)
    ]


def matrix_power(matrix, exponent):
    result = [[int(row == column) for column in range(3)] for row in range(3)]
    factor = [[int(value) % P for value in row] for row in matrix]
    exponent = int(exponent)
    while exponent:
        if exponent & 1:
            result = matrix_multiply(result, factor)
        factor = matrix_multiply(factor, factor)
        exponent //= 2
    return result


def projective_scalar(left, right):
    left_flat = [int(value) % P for row in left for value in row]
    right_flat = [int(value) % P for row in right for value in row]
    pivot = next((index for index, value in enumerate(right_flat) if value), None)
    if pivot is None:
        raise HSRank2KChainComparisonError("a projective matrix is zero")
    scalar = left_flat[pivot] * pow(right_flat[pivot], -1, P) % P
    if not scalar or any(
        (a - scalar * b) % P for a, b in zip(left_flat, right_flat)
    ):
        raise HSRank2KChainComparisonError("two group matrices are not projective")
    return scalar


def canonical_projective(vector):
    vector = [int(value) % P for value in vector]
    pivot = next((value for value in vector if value), None)
    if pivot is None:
        raise HSRank2KChainComparisonError("a frame pullback vanished")
    inverse = pow(pivot, -1, P)
    return [value * inverse % P for value in vector]


def polynomial_add(target, source, scale=1):
    for exponent, coefficient in source.items():
        value = (target.get(exponent, 0) + int(scale) * int(coefficient)) % P
        if value:
            target[exponent] = value
        else:
            target.pop(exponent, None)


def polynomial_multiply(left, right):
    result = {}
    for a, ca in left.items():
        for b, cb in right.items():
            exponent = tuple(int(x) + int(y) for x, y in zip(a, b))
            polynomial_add(result, {exponent: int(ca) * int(cb)})
    return result


def polynomial_power(polynomial, exponent):
    result = {(0, 0, 0): 1}
    for _ in range(int(exponent)):
        result = polynomial_multiply(result, polynomial)
    return result


def cubic_preservation_scalar(matrix):
    """Return c in F(M[Z,X,Y])=c F(Z,X,Y)."""

    linear = []
    for row in matrix:
        linear.append(
            {
                exponent: int(coefficient) % P
                for exponent, coefficient in zip(
                    ((1, 0, 0), (0, 1, 0), (0, 0, 1)), row
                )
                if int(coefficient) % P
            }
        )
    transformed = {}
    polynomial_add(
        transformed,
        polynomial_multiply(linear[0], polynomial_power(linear[2], 2)),
    )
    polynomial_add(transformed, polynomial_power(linear[1], 3), -1)
    polynomial_add(transformed, polynomial_power(linear[0], 3), -1)
    base = {(1, 0, 2): 1, (0, 3, 0): -1 % P, (3, 0, 0): -1 % P}
    scalar = transformed.get((1, 0, 2), 0) % P
    if not scalar or transformed != {
        exponent: scalar * coefficient % P for exponent, coefficient in base.items()
    }:
        raise HSRank2KChainComparisonError(
            "a combined factor matrix stopped preserving the plane cubic"
        )
    return scalar


def apply_tensor_degree_one(vector, factor_matrices):
    """Apply pullback to the 3^4 coefficient tensor."""

    result = list(map(lambda value: int(value) % P, vector))
    if len(result) != 81:
        raise ValueError("a theta coefficient vector must have length 81")
    for factor, matrix in enumerate(factor_matrices):
        # Coefficient pullback is M^T when point coordinates transform by M.
        operator = [[int(matrix[source][target]) % P for source in range(3)] for target in range(3)]
        stride = 3 ** (3 - factor)
        block = 3 * stride
        transformed = [0] * 81
        for base in range(0, 81, block):
            for offset in range(stride):
                for target in range(3):
                    transformed[base + target * stride + offset] = sum(
                        operator[target][source]
                        * result[base + source * stride + offset]
                        for source in range(3)
                    ) % P
        result = transformed
    return result


def vector_projective_scalar(left, right):
    pivot = next((index for index, value in enumerate(right) if int(value) % P), None)
    if pivot is None:
        raise HSRank2KChainComparisonError("a theta vector is zero")
    scalar = int(left[pivot]) * pow(int(right[pivot]), -1, P) % P
    if not scalar or any(
        (int(a) - scalar * int(b)) % P for a, b in zip(left, right)
    ):
        raise HSRank2KChainComparisonError("a theta section lost its eigenline")
    return scalar


def affine_form(row, source_bit):
    """Return coefficients of 1,u,v in the source chart frame."""

    z, x, y = map(lambda value: int(value) % P, row)
    return [z, x, y] if int(source_bit) == 0 else [y, x, z]


def local_form_as_ledger(row, source_bit):
    coefficients = affine_form(row, source_bit)
    basis = ("1", "x", "y") if int(source_bit) == 0 else (
        "y^-1",
        "x*y^-1",
        "1",
    )
    return [
        {"basis_monomial": monomial, "coefficient_mod_31": coefficient}
        for monomial, coefficient in zip(basis, map(lambda value: int(value) % P, row))
        if coefficient
    ]


def load_inputs():
    required = (
        DESCENT,
        K_ELEMENTS,
        K_STABLE_CHARTS,
        K_REFINED_INTERSECTIONS,
        K_REFINED_EDGE_PULLBACKS,
        THETA_COEFFICIENTS,
        LOCAL_SEEDS,
        EXTENSION_CERTIFICATE,
        EXTENSION_ARTIFACT,
    )
    for path in required:
        if not path.is_file():
            raise HSRank2KChainComparisonError(f"missing exact input: {path}")
    descent = load_json(DESCENT)
    extension = load_json(EXTENSION_CERTIFICATE)
    stable = read_csv(K_STABLE_CHARTS)
    elements = read_csv(K_ELEMENTS)
    refined = read_csv(K_REFINED_INTERSECTIONS)
    edge_pullbacks = read_csv(K_REFINED_EDGE_PULLBACKS)
    if descent.get("status") != "certified":
        raise HSRank2KChainComparisonError("the descent certificate changed")
    if descent.get("checks", {}).get("translation_group_law_checks") != 81:
        raise HSRank2KChainComparisonError("the factor translation law changed")
    if len(elements) != 81 or len(stable) != 1296:
        raise HSRank2KChainComparisonError("the 81-element chart orbit is incomplete")
    if len(refined) != 1024 or len(edge_pullbacks) != 2048:
        raise HSRank2KChainComparisonError("the generator refinement ledger changed")
    module = extension.get("module", {})
    if module.get("relation_shape") != [688, 2750]:
        raise HSRank2KChainComparisonError("the rank-two relation shape changed")
    checks = extension.get("exact_checks", {})
    if not checks.get("therefore_D2_D3_zero_over_S") or not checks.get(
        "therefore_lambda_is_a_closed_Ext2_cochain"
    ):
        raise HSRank2KChainComparisonError("the whole-ring D3/lambda check changed")

    theta = {label: [0] * 81 for label in range(6)}
    for row in read_csv(THETA_COEFFICIENTS):
        theta[int(row["label"])][int(row["theta_index"], 3)] = int(
            row["coefficient_mod_31"]
        ) % P
    if any(not any(vector) for vector in theta.values()):
        raise HSRank2KChainComparisonError("a theta section is empty")

    seed = next(
        (row for row in read_csv(LOCAL_SEEDS) if row.get("surface") == SURFACE),
        None,
    )
    if seed is None:
        raise HSRank2KChainComparisonError("the Z_2 seed is absent")
    labels = list(map(int, json.loads(seed["surface_section_labels"])))
    if labels != [2, 3]:
        raise HSRank2KChainComparisonError(
            "the exact Z_2 theta-equation labels changed"
        )
    return descent, extension, stable, elements, refined, edge_pullbacks, theta, labels


def group_records(descent, theta):
    rows = [list(map(int, row)) for row in descent["kernel_rows"]]
    translations = {
        tuple(map(int, key.split(","))): [
            [int(value) % P for value in row] for row in matrix
        ]
        for key, matrix in descent["translation_matrices"].items()
    }
    eigenvalues = list(map(int, descent["eigenline"]["eigenvalues"]))
    exponents_list = list(itertools.product(range(3), repeat=4))
    records = []
    for index, exponents in enumerate(exponents_list):
        flat = [
            sum(exponents[generator] * rows[generator][column] for generator in range(4))
            % 3
            for column in range(8)
        ]
        matrices = []
        keys = []
        for factor in range(4):
            matrix = [[int(row == column) for column in range(3)] for row in range(3)]
            for generator, exponent in enumerate(exponents):
                key = (
                    rows[generator][2 * factor],
                    rows[generator][2 * factor + 1],
                )
                matrix = matrix_multiply(
                    matrix, matrix_power(translations[key], exponent)
                )
            matrices.append(matrix)
            keys.append((flat[2 * factor], flat[2 * factor + 1]))

        theta_characters = {}
        for label, vector in theta.items():
            theta_characters[label] = vector_projective_scalar(
                apply_tensor_degree_one(vector, matrices), vector
            )
        expected_theta0 = 1
        for exponent, eigenvalue in zip(exponents, eigenvalues):
            expected_theta0 = expected_theta0 * pow(eigenvalue, exponent, P) % P
        if theta_characters[0] != expected_theta0:
            raise HSRank2KChainComparisonError(
                "the combined theta_0 character disagrees with descent_section"
            )
        normalization = pow(theta_characters[0], -1, P)
        records.append(
            {
                "index": index,
                "exponents": tuple(exponents),
                "flat": tuple(flat),
                "keys": tuple(keys),
                "matrices": matrices,
                "theta_characters": theta_characters,
                "H_normalization": normalization,
                "curve_scalars": [cubic_preservation_scalar(matrix) for matrix in matrices],
            }
        )
    return records


def verify_cached_orbit(records, stable, elements, refined):
    stable_lookup = {
        (tuple(json.loads(row["K_exponents"])), row["base_chart"]): json.loads(
            row["pulled_back_factor_linear_forms_Z_X_Y"]
        )
        for row in stable
    }
    if len(stable_lookup) != 1296:
        raise HSRank2KChainComparisonError("the cached stable chart keys repeat")
    checks = 0
    for record in records:
        for bits in itertools.product((0, 1), repeat=4):
            chart = "C" + "".join(map(str, bits))
            expected = [
                canonical_projective(
                    record["matrices"][factor][2 if bit else 0]
                )
                for factor, bit in enumerate(bits)
            ]
            if stable_lookup[(record["exponents"], chart)] != expected:
                raise HSRank2KChainComparisonError(
                    "the recomputed K-stable frame orbit changed"
                )
            checks += 1

    element_lookup = {
        tuple(json.loads(row["K_generator_exponents"])): int(row["K_element_index"])
        for row in elements
    }
    if any(element_lookup[record["exponents"]] != record["index"] for record in records):
        raise HSRank2KChainComparisonError("the cached K element ordering changed")

    refined_lookup = {
        (
            int(row["kernel_generator"]),
            row["source_chart"],
            row["target_chart"],
        ): json.loads(row["localized_nonvanishing_forms"])
        for row in refined
    }
    generator_checks = 0
    for generator in range(4):
        exponents = tuple(int(index == generator) for index in range(4))
        record = records[element_lookup[exponents]]
        expected = [
            local_form_as_ledger(
                record["matrices"][factor][2 if BASE_BITS[factor] else 0],
                BASE_BITS[factor],
            )
            for factor in range(4)
        ]
        if refined_lookup[(generator, BASE_CHART, BASE_CHART)] != expected:
            raise HSRank2KChainComparisonError(
                "a generator-localization row disagrees with the group matrices"
            )
        generator_checks += 1
    return checks, generator_checks


def verify_group_coherence(records):
    lookup = {record["exponents"]: record for record in records}
    product_checks = 0
    factor_map_checks = 0
    normalized_H_checks = 0
    normalized_O4_checks = 0
    normalized_surface_character_checks = 0
    for left in records:
        for right in records:
            target_exp = tuple(
                (a + b) % 3
                for a, b in zip(left["exponents"], right["exponents"])
            )
            target = lookup[target_exp]
            factor_scalars = []
            for left_matrix, right_matrix, target_matrix in zip(
                left["matrices"], right["matrices"], target["matrices"]
            ):
                scalar = projective_scalar(
                    matrix_multiply(left_matrix, right_matrix), target_matrix
                )
                factor_scalars.append(scalar)
                factor_map_checks += 1
            product_scalar = 1
            for scalar in factor_scalars:
                product_scalar = product_scalar * scalar % P
            if (
                left["H_normalization"]
                * right["H_normalization"]
                * product_scalar
                - target["H_normalization"]
            ) % P:
                raise HSRank2KChainComparisonError(
                    "the normalized H frame action lost the group law"
                )
            normalized_H_checks += 1
            if (
                pow(left["H_normalization"], 4, P)
                * pow(right["H_normalization"], 4, P)
                * pow(product_scalar, 4, P)
                - pow(target["H_normalization"], 4, P)
            ) % P:
                raise HSRank2KChainComparisonError(
                    "the normalized O(4H) factor lost the group law"
                )
            normalized_O4_checks += 1
            for label in (2, 3):
                left_relative = (
                    left["theta_characters"][label]
                    * left["H_normalization"]
                    % P
                )
                right_relative = (
                    right["theta_characters"][label]
                    * right["H_normalization"]
                    % P
                )
                target_relative = (
                    target["theta_characters"][label]
                    * target["H_normalization"]
                    % P
                )
                if left_relative * right_relative % P != target_relative:
                    raise HSRank2KChainComparisonError(
                        "a normalized Z_2 theta character lost the group law"
                    )
                normalized_surface_character_checks += 1
            product_checks += 1
    return {
        "group_product_checks": product_checks,
        "projective_factor_matrix_checks": factor_map_checks,
        "normalized_H_cocycle_checks": normalized_H_checks,
        "normalized_O4H_cocycle_checks": normalized_O4_checks,
        "normalized_Z2_theta_character_checks": normalized_surface_character_checks,
    }


def orbit_ring_rows(records):
    rows = []
    for record in records:
        for factor, (matrix, source_bit, target_bit) in enumerate(
            zip(record["matrices"], BASE_BITS, BASE_BITS)
        ):
            frame_index = 2 if target_bit else 0
            other_index = 0 if target_bit else 2
            rows.append(
                {
                    "K_element_index": record["index"],
                    "K_exponents": json.dumps(list(record["exponents"])),
                    "orbit_chart": f"K{record['index']:02d}_{BASE_CHART}",
                    "factor": factor,
                    "source_chart_bit": source_bit,
                    "target_chart_bit": target_bit,
                    "factor_point_matrix_Z_X_Y": json.dumps(matrix, separators=(",", ":")),
                    "target_frame_pullback_homogeneous_Z_X_Y": json.dumps(
                        matrix[frame_index], separators=(",", ":")
                    ),
                    "target_frame_pullback_in_source_1_u_v": json.dumps(
                        affine_form(matrix[frame_index], source_bit), separators=(",", ":")
                    ),
                    "target_u_numerator_in_source_1_u_v": json.dumps(
                        affine_form(matrix[1], source_bit), separators=(",", ":")
                    ),
                    "target_v_numerator_in_source_1_u_v": json.dumps(
                        affine_form(matrix[other_index], source_bit), separators=(",", ":")
                    ),
                    "target_affine_substitution": (
                        "u'=u_numerator/frame_pullback;"
                        "v'=v_numerator/frame_pullback"
                    ),
                    "plane_cubic_pullback_scalar_mod_31": record["curve_scalars"][factor],
                }
            )
    return rows


def surface_rows(records, labels):
    rows = []
    for record in records:
        raw0 = record["theta_characters"][0]
        rows.append(
            {
                "K_element_index": record["index"],
                "K_exponents": json.dumps(list(record["exponents"])),
                "surface": SURFACE,
                "theta_equation_labels": json.dumps(labels),
                "raw_theta0_character_mod_31": raw0,
                "H_linearization_normalization_mod_31": record["H_normalization"],
                "theta2_raw_character_mod_31": record["theta_characters"][2],
                "theta3_raw_character_mod_31": record["theta_characters"][3],
                "theta2_normalized_character_mod_31": (
                    record["theta_characters"][2] * record["H_normalization"] % P
                ),
                "theta3_normalized_character_mod_31": (
                    record["theta_characters"][3] * record["H_normalization"] % P
                ),
                "four_curve_pullback_scalars_mod_31": json.dumps(
                    record["curve_scalars"]
                ),
                "surface_preserved": True,
            }
        )
    return rows


def cayley_edge_rows(records):
    lookup = {record["exponents"]: record for record in records}
    rows = []
    for source in records:
        for generator in range(4):
            target_exp = tuple(
                (value + int(index == generator)) % 3
                for index, value in enumerate(source["exponents"])
            )
            target = lookup[target_exp]
            scalar = (
                target["H_normalization"]
                * pow(source["H_normalization"], -1, P)
                % P
            )
            source_frames = [
                affine_form(
                    source["matrices"][factor][2 if bit else 0], bit
                )
                for factor, bit in enumerate(BASE_BITS)
            ]
            target_frames = [
                affine_form(
                    target["matrices"][factor][2 if bit else 0], bit
                )
                for factor, bit in enumerate(BASE_BITS)
            ]
            rows.append(
                {
                    "job_id": len(rows),
                    "source_K_index": source["index"],
                    "target_K_index": target["index"],
                    "generator": generator,
                    "source_K_exponents": json.dumps(list(source["exponents"])),
                    "target_K_exponents": json.dumps(list(target["exponents"])),
                    "localized_overlap": (
                        "invert all source and target factor-frame forms listed"
                    ),
                    "source_frame_forms_in_base_1_u_v": json.dumps(
                        source_frames, separators=(",", ":")
                    ),
                    "target_frame_forms_in_base_1_u_v": json.dumps(
                        target_frames, separators=(",", ":")
                    ),
                    "normalized_H_scalar_mod_31": scalar,
                    "normalized_H_factor_circuit": (
                        "scalar*product(target_frame_forms)/"
                        "product(source_frame_forms)"
                    ),
                    "normalized_O4H_scalar_mod_31": pow(scalar, 4, P),
                    "normalized_O4H_factor_circuit": (
                        "normalized_H_factor_circuit^4"
                    ),
                    "whole_ring_J1_emitted": False,
                    "whole_ring_J2_emitted": False,
                    "whole_ring_lambda_correction_emitted": False,
                    "L0_A_source_equals_A_target_L1_certified": False,
                }
            )
    if len(rows) != 324:
        raise HSRank2KChainComparisonError("the K Cayley edge census changed")
    return rows


def readiness_rows(records):
    lookup = {record["exponents"]: record for record in records}
    rows = []
    for generator in range(4):
        exponents = tuple(int(index == generator) for index in range(4))
        record = lookup[exponents]
        rows.append(
            {
                "kernel_generator": generator,
                "target_K_index": record["index"],
                "rational_ring_map": "exact_and_emitted",
                "six_surface_equation_transport": "exact_and_emitted",
                "O4H_factor": "exact_normalized_circuit_emitted",
                "source_A_688x2750": "exact_bound_npz",
                "independently_normalized_target_A_688x2750": "not_emitted",
                "J1_687x687": "not_emitted",
                "J2_2750x2750": "not_emitted",
                "lambda_correction_a_1x687": "not_emitted",
                "comparison_verdict": "not_materializable_from_cached_matrices_alone",
            }
        )
    return rows


def build(output_directory: Path = OUTPUT_DIRECTORY):
    (
        descent,
        extension,
        stable,
        elements,
        refined,
        edge_pullbacks,
        theta,
        labels,
    ) = load_inputs()
    records = group_records(descent, theta)
    stable_checks, generator_refinement_checks = verify_cached_orbit(
        records, stable, elements, refined
    )
    coherence = verify_group_coherence(records)
    ring_rows = orbit_ring_rows(records)
    equation_rows = surface_rows(records, labels)
    edge_rows = cayley_edge_rows(records)
    ready_rows = readiness_rows(records)

    output_directory.mkdir(parents=True, exist_ok=True)
    artifacts = {}
    for name, rows in (
        ("z2_k_orbit_ring_maps.csv", ring_rows),
        ("z2_surface_equation_transport.csv", equation_rows),
        ("z2_k_cayley_chain_comparison_jobs.csv", edge_rows),
        ("generator_comparison_readiness.csv", ready_rows),
    ):
        path = output_directory / name
        write_csv(path, rows)
        artifacts[name] = path

    contract = {
        "schema": "section12-Z2-rank2-K-chain-comparison-input-contract-v1",
        "field": "F_31",
        "base_presentation": {
            "surface": SURFACE,
            "chart": BASE_CHART,
            "ring": (
                "S=F_31[u0,v0,...,u3,v3]/(four affine cubic equations,"
                " theta_2,theta_3)"
            ),
            "A": "[lambda; -D2]:S^2750 -> S^688",
            "A_shape": [688, 2750],
            "D3_source_rank": 9610,
            "kernel_rank_on_rank_two_locus": 2064,
            "lambda_closed": True,
        },
        "already_exact": {
            "K_group": "(Z/3)^4 with 81 exact projective factor actions",
            "orbit_ring_maps": 81,
            "factor_affine_substitutions": len(ring_rows),
            "Z2_equation_character_rows": len(equation_rows),
            "typed_generator_edges": len(edge_rows),
            "O4H_factor": (
                "w_k=(n_k*product_i(k^*F_i/F_i))^4, where n_k is the "
                "inverse theta_0 character; this normalized factor obeys the "
                "K group law"
            ),
        },
        "comparison_equations": {
            "target_identity": "L0_k A_source = A_target L1_k",
            "bottom_square": "J1_k D2_source = D2_target J2_k",
            "counit_square": (
                "lambda_target J2_k - w_k lambda_source = a_k D2_source"
            ),
            "L0_block": "[[w_k,-a_k],[0,J1_k]]",
            "L1_block": "J2_k",
            "shapes": {
                "J1_k": [687, 687],
                "J2_k": [2750, 2750],
                "a_k": [1, 687],
                "L0_k": [688, 688],
                "L1_k": [2750, 2750],
            },
        },
        "smallest_unemitted_exact_data": [
            {
                "name": "translated_tree_presentations",
                "required": (
                    "an independently normalized A_target on each generator "
                    "overlap, obtained by transporting the raw graph/Shamash "
                    "carrier and then applying the stored tree contraction"
                ),
                "why_not_phi_A": (
                    "defining A_target=phi_k(A_source) makes identity matrices "
                    "a tautological base-change comparison and does not choose "
                    "a K-linearization of the rank-two module"
                ),
            },
            {
                "name": "tree_chain_comparison",
                "required": (
                    "sparse localized J1_k and J2_k satisfying the bottom "
                    "presentation square and compatible with D3"
                ),
                "derivation_inputs_already_present": (
                    "the emitted rational coordinate substitutions, the six "
                    "diagonal equation characters, the full graph exterior "
                    "incidence, and the stored tree inclusion/projection"
                ),
            },
            {
                "name": "one_counit_correction_per_typed_edge",
                "required": (
                    "one row a_k solving the displayed 2750-column equation; "
                    "the six surface-equation columns are constraints on that "
                    "same row, not six independent repairs"
                ),
            },
            {
                "name": "K_group_coherence",
                "required": (
                    "compare every typed pair composite (or, equivalently, all "
                    "generator commutation squares and order-three cycles); "
                    "emit equality when strict and an explicit chain homotopy "
                    "when the tree comparison is only homotopy-natural"
                ),
                "Cayley_edges": 324,
                "generator_square_cells": 81 * 6,
                "distinct_generator_order_three_cycles": 81 * 4 // 3,
                "all_typed_pair_composites": 81**2,
            },
        ],
        "required_sparse_output_schema": {
            "matrix_rows": (
                "source_K_index,target_K_index,generator,matrix_name,row,column,"
                "coefficient_mod_31,laurent_numerator_exponents,"
                "frame_denominator_powers"
            ),
            "allowed_matrix_names": ["A_target", "J1", "J2", "a", "L0", "L1"],
            "identity_receipt": (
                "one zero-residual hash for L0*A_source-A_target*L1 on every "
                "typed edge, plus D3 compatibility, commutation-square, and "
                "order-three-cycle coherence hashes"
            ),
        },
        "verdict": {
            "one_nonidentity_generator_L0_L1_currently_certified": False,
            "reason": (
                "the ring maps, equation characters, and O(4H) factor are exact, "
                "but no nonidentity whole-ring J1/J2 or counit-correction row is "
                "cached; the source presentation alone does not determine a "
                "chosen equivariant module isomorphism"
            ),
            "next_computation_is_well_posed": True,
        },
    }
    contract_path = output_directory / "k_chain_comparison_input_contract.json"
    contract_path.write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    artifacts[contract_path.name] = contract_path

    provenance = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (
            DESCENT,
            K_ELEMENTS,
            K_STABLE_CHARTS,
            K_REFINED_INTERSECTIONS,
            K_REFINED_EDGE_PULLBACKS,
            THETA_COEFFICIENTS,
            LOCAL_SEEDS,
            EXTENSION_CERTIFICATE,
            EXTENSION_ARTIFACT,
        )
    }
    certificate = {
        "schema": "section12-Z2-rank2-K-chain-comparison-audit-v1",
        "status": (
            "81-orbit-ring-maps-and-Z2-equation-characters-certified;"
            "normalized-O4H-group-cocycle-certified;324-chain-jobs-indexed;"
            "nonidentity-whole-ring-L0-L1-open"
        ),
        "field": "F_31",
        "counts": {
            "K_elements": len(records),
            "K_stable_chart_rows_replayed": stable_checks,
            "generator_refinement_rows_replayed": generator_refinement_checks,
            "factor_affine_map_rows": len(ring_rows),
            "Z2_equation_rows": len(equation_rows),
            "typed_Cayley_edges": len(edge_rows),
            **coherence,
        },
        "exact_line_factor": {
            "H_normalization_source": "inverse raw theta_0 eigencharacter",
            "H_factor": "n_k*product_i(k^*F_i/F_i)",
            "relative_extension_factor": "H_factor^4",
            "group_coherent": True,
        },
        "whole_ring_boundary": {
            "source_A_D3_lambda_bound": True,
            "translated_ring_maps_bound": True,
            "translated_surface_equation_characters_bound": True,
            "nonidentity_J1_J2_bound": False,
            "nonidentity_counit_correction_bound": False,
            "nonidentity_L0_L1_identity_bound": False,
            "supportwise_evaluation_used_as_substitute": False,
            "tautological_phi_A_identity_used_as_linearization": False,
        },
        "provenance": provenance,
        "artifacts": {},
    }
    for name, path in artifacts.items():
        certificate["artifacts"][name] = sha256(path)
    certificate_path = output_directory / "k_chain_comparison_certificate.json"
    certificate_path.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    artifacts[certificate_path.name] = certificate_path

    manifest = {name: sha256(path) for name, path in sorted(artifacts.items())}
    manifest_path = output_directory / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return certificate_path


if __name__ == "__main__":
    print(build())
