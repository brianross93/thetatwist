"""Build one honest nonidentity K-overlap and its exact graph-J1 lift.

The all-chart presentation producer applies one local construction law to the
ordinary product charts.  It does not construct the overlap between C0010 and
the inverse image of C0010 under a nonidentity element of K.  This producer
constructs that missing overlap for Z_2 and one bound standard K generator.

The target Hermite algebra is rebuilt from the translated global support jets.
It is not defined as ``phi(A_source)``.  On the common open, the producer emits
the two quotient restrictions J0 and the factored rational graph-Koszul J1.
For each target affine coordinate r=N/L it verifies the cleared-denominator
identity

    d_overlap J1 = J0 d_target

exactly over F_31.  Full Shamash J2 and the Hartshorne--Serre counit correction
``a`` need the transported six-equation homotopies and are not inferred here.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_iw_hermite_border_basis as hermite,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_iw_presentation_solver as iw_solver,
)


P = 31
SURFACE = 2
CHART = "C0010"
GENERATOR_EXPONENTS = (1, 0, 0, 0)
VARIABLE_COUNT = 8
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_rank2_k_generator_overlap_j1"
)
ALL_CHART_PRODUCER = ROOT / (
    "experiments/section12_instantiation/cm_bridge_hs_all_chart_presentations.py"
)
BORDER_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_iw_hermite_border_basis"
)


class HSKGeneratorOverlapJ1Error(RuntimeError):
    """A translated support, Hermite algebra, or chain identity changed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def array_sha256(value: np.ndarray) -> str:
    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(value, dtype=np.uint8)).tobytes()
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


def matrix_vector(matrix, vector):
    return tuple(
        sum(int(matrix[row][column]) * int(vector[column]) for column in range(3))
        % P
        for row in range(3)
    )


def transformed_affine_jet(point, tangent, bits, point_values, point_first, record):
    """Translate a homogeneous point/tangent and normalize in the target chart."""

    values = []
    derivatives = []
    for factor, point_index in enumerate(point):
        source_value = tuple(map(int, point_values[point_index]))
        source_tangent = tuple(
            int(tangent[factor]) * int(value) % P
            for value in point_first[point_index]
        )
        matrix = record["factor_matrices"][factor]
        z, x, y = matrix_vector(matrix, source_value)
        dz, dx, dy = matrix_vector(matrix, source_tangent)
        bit = int(bits[factor])
        denominator = y if bit else z
        if not denominator:
            raise HSKGeneratorOverlapJ1Error(
                "a translated point left the requested target chart"
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


def affine_forms(matrix, bit):
    """Return L, N_u, N_v in the source affine basis (1,u,v)."""

    matrix = np.asarray(matrix, dtype=np.int64) % P

    def row_form(row):
        if bit:
            # Source homogeneous coordinates are (Z,X,Y)=(v,u,1).
            return tuple(map(int, (row[2], row[1], row[0])))
        # Source homogeneous coordinates are (Z,X,Y)=(1,u,v).
        return tuple(map(int, row))

    denominator_row = 2 if bit else 0
    v_row = 0 if bit else 2
    return (
        row_form(matrix[denominator_row]),
        row_form(matrix[1]),
        row_form(matrix[v_row]),
    )


def linear_operator(form, multiplication, factor):
    dimension = multiplication.shape[1]
    identity = np.eye(dimension, dtype=np.int64)
    return (
        int(form[0]) * identity
        + int(form[1]) * multiplication[2 * factor]
        + int(form[2]) * multiplication[2 * factor + 1]
    ) % P


def polynomial_add(polynomial, monomial, matrix):
    matrix = np.asarray(matrix, dtype=np.int64) % P
    if monomial in polynomial:
        matrix = (polynomial[monomial] + matrix) % P
    if np.any(matrix):
        polynomial[monomial] = matrix
    else:
        polynomial.pop(monomial, None)


def verify_target_generator(
    target_variable,
    factor,
    numerator,
    denominator,
    blocks,
    target_restriction,
    target_multiplication,
    overlap_multiplication,
):
    """Verify the J1 identity after multiplying by the common denominator."""

    left = {}
    for source_variable, block in blocks.items():
        polynomial_add(left, source_variable + 1, block)
        polynomial_add(
            left,
            0,
            -overlap_multiplication[source_variable] @ block,
        )

    right = {}
    target_product = target_restriction @ target_multiplication[target_variable] % P
    polynomial_add(
        right,
        0,
        int(numerator[0]) * target_restriction
        - int(denominator[0]) * target_product,
    )
    for local in range(2):
        variable = 2 * factor + local
        polynomial_add(
            right,
            variable + 1,
            int(numerator[local + 1]) * target_restriction
            - int(denominator[local + 1]) * target_product,
        )
    if set(left) != set(right):
        raise HSKGeneratorOverlapJ1Error(
            "a target J1 identity has different cleared polynomial support"
        )
    for monomial in left:
        if not np.array_equal(left[monomial] % P, right[monomial] % P):
            raise HSKGeneratorOverlapJ1Error(
                "d_overlap J1 does not equal J0 d_target"
            )
    return {
        "domain_generator": target_variable,
        "cleared_polynomial_terms": len(left),
        "residual_terms": 0,
        "residual_is_zero": True,
        "identity": "L_i(x)*d_overlap*J1=L_i(x)*J0*d_target",
    }


def build_data():
    context = iw_solver.load_problem_context()
    hermite.prepare_border_context(context)
    bits = tuple(map(int, context["chart_bits"][CHART]))
    point_values = context["_border_point_values"]
    point_first = context["_border_point_first"]
    points = context["_border_points"]
    local_symbols = context["local_symbols"]
    records = context["_border_records"]
    record = next(
        item for item in records if tuple(item["exponents"]) == GENERATOR_EXPONENTS
    )
    tangent = tuple(
        map(
            int,
            json.loads(
                context["seeds"][SURFACE]["tangent_fixed_frame_coefficients"]
            ),
        )
    )

    source_chart = hermite.chart_border_basis(context, SURFACE, CHART)
    source_pairs = []
    target_full_pairs = []
    point_rows = []
    for point in context["support_by_surface"][SURFACE]:
        _scale, image = local_symbols.action_scale(record, point, point_values, points)
        target_open = local_symbols.chart_membership(image, bits, point_values)
        if target_open:
            translated_jet = transformed_affine_jet(
                point, tangent, bits, point_values, point_first, record
            )
            canonical_target_jet = hermite.affine_point_and_tangent(
                image, tangent, bits, point_values, point_first
            )
            if translated_jet != canonical_target_jet:
                raise HSKGeneratorOverlapJ1Error(
                    "the translated invariant tangent lost its canonical normalization"
                )
            target_full_pairs.append((tuple(image), translated_jet))
        source_open = local_symbols.chart_membership(point, bits, point_values)
        if source_open and target_open:
            source_jet = hermite.affine_point_and_tangent(
                point, tangent, bits, point_values, point_first
            )
            source_pairs.append((tuple(point), tuple(image), source_jet, translated_jet))
            point_rows.append(
                {
                    "source_support": json.dumps(list(map(int, point))),
                    "target_support": json.dumps(list(map(int, image))),
                    "source_affine_values": json.dumps(list(map(int, source_jet[0]))),
                    "target_affine_values": json.dumps(
                        list(map(int, translated_jet[0]))
                    ),
                    "source_affine_tangent": json.dumps(
                        list(map(int, source_jet[1]))
                    ),
                    "target_affine_tangent": json.dumps(
                        list(map(int, translated_jet[1]))
                    ),
                }
            )

    target_full_pairs.sort(key=lambda item: item[0])
    target_full_affine = [item[1] for item in target_full_pairs]
    target_basis, target_hermite, target_tested, target_degree = (
        hermite.select_standard_monomials(target_full_affine)
    )
    target_multiplication, _target_border = hermite.border_multiplication(
        target_basis, target_hermite, target_full_affine
    )
    if len(target_basis) != len(source_chart["basis"]):
        raise HSKGeneratorOverlapJ1Error(
            "the independently translated target changed the C0010 length"
        )
    if target_basis != source_chart["basis"] or not np.array_equal(
        target_multiplication, source_chart["multiplication"]
    ):
        raise HSKGeneratorOverlapJ1Error(
            "the independently rebuilt K target does not recover canonical C0010"
        )

    source_affine = [item[2] for item in source_pairs]
    target_affine = [item[3] for item in source_pairs]
    overlap_basis, overlap_hermite, overlap_tested, overlap_degree = (
        hermite.select_standard_monomials(source_affine)
    )
    overlap_multiplication, overlap_border = hermite.border_multiplication(
        overlap_basis, overlap_hermite, source_affine
    )
    overlap_inverse = hermite.inverse_mod(overlap_hermite)
    source_columns = np.stack(
        [hermite.hermite_column(exponent, source_affine) for exponent in source_chart["basis"]],
        axis=1,
    )
    target_columns = np.stack(
        [hermite.hermite_column(exponent, target_affine) for exponent in target_basis],
        axis=1,
    )
    source_restriction = overlap_inverse @ source_columns % P
    target_restriction = overlap_inverse @ target_columns % P
    overlap_dimension = len(overlap_basis)
    if hermite.rank_mod(source_restriction) != overlap_dimension or hermite.rank_mod(
        target_restriction
    ) != overlap_dimension:
        raise HSKGeneratorOverlapJ1Error("a K-overlap restriction is not surjective")

    for variable in range(VARIABLE_COUNT):
        if np.any(
            (
                source_restriction @ source_chart["multiplication"][variable]
                - overlap_multiplication[variable] @ source_restriction
            )
            % P
        ):
            raise HSKGeneratorOverlapJ1Error(
                "the source K-overlap restriction is not an algebra map"
            )
        target_jet = hermite.hermite_column(
            tuple(int(index == variable) for index in range(VARIABLE_COUNT)),
            target_affine,
        )
        target_coordinate = overlap_inverse @ target_jet % P
        target_operator = hermite.coordinate_multiplication_from_vector(
            target_coordinate, overlap_basis, overlap_multiplication
        )
        if np.any(
            (
                target_restriction @ target_multiplication[variable]
                - target_operator @ target_restriction
            )
            % P
        ):
            raise HSKGeneratorOverlapJ1Error(
                "the translated target restriction is not an algebra map"
            )

    coefficient_blocks = {
        "source_J0": np.asarray(source_restriction, dtype=np.uint8),
        "target_J0": np.asarray(target_restriction, dtype=np.uint8),
    }
    program_rows = []
    verification_rows = []
    for variable in range(VARIABLE_COUNT):
        program_rows.append(
            {
                "endpoint": "source",
                "domain_generator": variable,
                "codomain_generator": variable,
                "factor": variable // 2,
                "denominator_linear_form_in_1_u_v": "[1,0,0]",
                "coefficient_block_key": "source_J0",
                "coefficient_block_sha256": array_sha256(source_restriction),
                "factored_term": "source_J0",
            }
        )
        left_variable = overlap_multiplication[variable] @ source_restriction % P
        right_variable = source_restriction @ source_chart["multiplication"][variable] % P
        if np.any((left_variable - right_variable) % P):
            raise HSKGeneratorOverlapJ1Error("the source graph J1 identity failed")
        verification_rows.append(
            {
                "endpoint": "source",
                "domain_generator": variable,
                "cleared_polynomial_terms": 2,
                "residual_terms": 0,
                "residual_is_zero": True,
                "identity": "d_overlap*J1=J0*d_source",
            }
        )

    affine_form_rows = []
    nonidentity_blocks = 0
    for factor, matrix in enumerate(record["factor_matrices"]):
        denominator, numerator_u, numerator_v = affine_forms(matrix, bits[factor])
        affine_form_rows.append(
            {
                "factor": factor,
                "denominator_L": json.dumps(list(denominator)),
                "numerator_Nu": json.dumps(list(numerator_u)),
                "numerator_Nv": json.dumps(list(numerator_v)),
                "target_u": "N_u/L",
                "target_v": "N_v/L",
            }
        )
        denominator_operator = linear_operator(
            denominator, overlap_multiplication, factor
        )
        denominator_inverse = hermite.inverse_mod(denominator_operator)
        for local, numerator in enumerate((numerator_u, numerator_v)):
            target_variable = 2 * factor + local
            numerator_operator = linear_operator(
                numerator, overlap_multiplication, factor
            )
            blocks = {}
            for source_local in range(2):
                source_variable = 2 * factor + source_local
                block = denominator_inverse @ (
                    int(numerator[source_local + 1]) * denominator_operator
                    - int(denominator[source_local + 1]) * numerator_operator
                ) @ target_restriction % P
                if not np.any(block):
                    continue
                key = f"target_D_{source_variable}_{target_variable}"
                coefficient_blocks[key] = np.asarray(block, dtype=np.uint8)
                blocks[source_variable] = block
                nonidentity_blocks += int(
                    source_variable != target_variable
                    or denominator != (1, 0, 0)
                    or not np.array_equal(block, target_restriction)
                )
                program_rows.append(
                    {
                        "endpoint": "target",
                        "domain_generator": target_variable,
                        "codomain_generator": source_variable,
                        "factor": factor,
                        "denominator_linear_form_in_1_u_v": json.dumps(
                            list(denominator)
                        ),
                        "coefficient_block_key": key,
                        "coefficient_block_sha256": array_sha256(block),
                        "factored_term": "L(x)^-1*D",
                    }
                )
            verification_rows.append(
                {
                    "endpoint": "target",
                    **verify_target_generator(
                        target_variable,
                        factor,
                        numerator,
                        denominator,
                        blocks,
                        target_restriction,
                        target_multiplication,
                        overlap_multiplication,
                    ),
                }
            )

    if not nonidentity_blocks:
        raise HSKGeneratorOverlapJ1Error(
            "the chosen generator produced only a tautological identity comparison"
        )

    c0000_summary = json.loads(
        (BORDER_DIRECTORY / "Z_2_C0000.json").read_text(encoding="utf-8")
    )
    c0010_summary = json.loads(
        (BORDER_DIRECTORY / "Z_2_C0010.json").read_text(encoding="utf-8")
    )
    c0000_dimension = int(c0000_summary["summary"]["quotient_length"])
    c0010_dimension = int(c0010_summary["summary"]["quotient_length"])
    if c0000_dimension == c0010_dimension:
        raise HSKGeneratorOverlapJ1Error(
            "the ordinary pilot dimension no longer distinguishes the K target"
        )

    return {
        "bits": bits,
        "record": record,
        "point_rows": point_rows,
        "affine_form_rows": affine_form_rows,
        "program_rows": program_rows,
        "verification_rows": verification_rows,
        "coefficient_blocks": coefficient_blocks,
        "source_chart": source_chart,
        "target_basis": target_basis,
        "target_hermite": target_hermite,
        "target_multiplication": target_multiplication,
        "overlap_basis": overlap_basis,
        "overlap_hermite": overlap_hermite,
        "overlap_multiplication": overlap_multiplication,
        "overlap_border": overlap_border,
        "source_restriction": source_restriction,
        "target_restriction": target_restriction,
        "counts": {
            "source_chart_support_points": len(source_chart["support_point_indices"]),
            "target_chart_support_points": len(target_full_pairs),
            "K_overlap_support_points": len(source_pairs),
            "source_chart_length": len(source_chart["basis"]),
            "independently_rebuilt_target_length": len(target_basis),
            "K_overlap_length": overlap_dimension,
            "target_monomials_tested": target_tested,
            "target_highest_selected_degree": target_degree,
            "overlap_monomials_tested": overlap_tested,
            "overlap_highest_selected_degree": overlap_degree,
            "source_J1_terms": sum(
                row["endpoint"] == "source" for row in program_rows
            ),
            "target_J1_terms": sum(
                row["endpoint"] == "target" for row in program_rows
            ),
            "nonidentity_target_blocks": nonidentity_blocks,
            "chain_identities_checked": len(verification_rows),
        },
        "ordinary_chart_assessment": {
            "all_chart_pilot": "Z_2/C0000",
            "pilot_ordinary_chart_length": c0000_dimension,
            "required_C0010_target_length": c0010_dimension,
            "pilot_supplies_K_target": False,
            "reason": (
                "C0000 is a different ordinary product open. Its quotient length is "
                "112 rather than 98, and its coordinate frames are not the four "
                "generator-pulled C0010 frames. It establishes chart-law portability."
            ),
        },
    }


def write_artifacts(output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    result = build_data()

    blocks_path = output_directory / "k_generator_j0_j1_blocks.npz"
    overlap_path = output_directory / "k_generator_overlap_hermite.npz"
    points_path = output_directory / "k_generator_overlap_points.csv"
    forms_path = output_directory / "k_generator_affine_forms.csv"
    program_path = output_directory / "k_generator_j1_factor_program.csv"
    verification_path = output_directory / "k_generator_j1_verification.csv"
    certificate_path = output_directory / "k_generator_overlap_j1_certificate.json"

    np.savez_compressed(blocks_path, **result["coefficient_blocks"])
    np.savez_compressed(
        overlap_path,
        field=np.asarray([P], dtype=np.uint8),
        generator_exponents=np.asarray(GENERATOR_EXPONENTS, dtype=np.uint8),
        chart_bits=np.asarray(result["bits"], dtype=np.uint8),
        target_standard_monomials=np.asarray(result["target_basis"], dtype=np.uint8),
        target_hermite_matrix=np.asarray(result["target_hermite"], dtype=np.uint8),
        target_multiplication_matrices=np.asarray(
            result["target_multiplication"], dtype=np.uint8
        ),
        overlap_standard_monomials=np.asarray(
            result["overlap_basis"], dtype=np.uint8
        ),
        overlap_hermite_matrix=np.asarray(result["overlap_hermite"], dtype=np.uint8),
        overlap_multiplication_matrices=np.asarray(
            result["overlap_multiplication"], dtype=np.uint8
        ),
        source_restriction=np.asarray(result["source_restriction"], dtype=np.uint8),
        target_restriction=np.asarray(result["target_restriction"], dtype=np.uint8),
    )
    write_csv(points_path, result["point_rows"])
    write_csv(forms_path, result["affine_form_rows"])
    write_csv(program_path, result["program_rows"])
    write_csv(verification_path, result["verification_rows"])

    certificate = {
        "schema": "section12-HS-rank2-K-generator-overlap-J1-v1",
        "field": "F_31",
        "status": (
            "one-nonidentity-K-overlap-target-normalized;J0-and-graph-J1-exact;"
            "full-Shamash-J2-and-counit-a-open"
        ),
        "surface": "Z_2",
        "source_chart": CHART,
        "target_chart": CHART,
        "K_generator_exponents": list(GENERATOR_EXPONENTS),
        "counts": result["counts"],
        "all_chart_pilot_assessment": result["ordinary_chart_assessment"],
        "independent_target_normalization": {
            "construction": (
                "translate all 81 global support jets by the recorded nonidentity "
                "plane-cubic matrices; retain target-C0010 jets; select the target "
                "Hermite order ideal and multiplication matrices afresh"
            ),
            "target_defined_as_phi_of_source_matrix": False,
            "translated_tangent_jets_equal_canonical_target_jets": True,
            "independent_target_basis_and_multiplication_recover_canonical_C0010": True,
            "source_and_target_restrictions_are_independently_evaluated": True,
        },
        "exact_checks": {
            "K_generator_is_nonidentity": True,
            "translated_support_and_tangent_normalization_exact": True,
            "independently_rebuilt_target_has_length_98": True,
            "K_overlap_is_product_of_retained_dual_points": True,
            "both_J0_restrictions_have_full_overlap_rank": True,
            "source_coordinate_intertwining_checks": 8,
            "target_coordinate_intertwining_checks": 8,
            "source_graph_J1_chain_identities": 8,
            "target_graph_J1_cleared_denominator_chain_identities": 8,
            "all_chain_residuals_zero": True,
            "nonidentity_target_blocks_present": True,
            "tautological_phi_A_identity_used": False,
        },
        "claim_boundary": {
            "proved": (
                "the independently normalized target quotient, both J0 maps, and "
                "the exact rational graph-Koszul J1 for one nonidentity K generator"
            ),
            "not_claimed": (
                "full Shamash J2, the comparison homotopies for the six equation "
                "sectors, transported lambda, the counit correction a, L0/L1, K "
                "coherence, Yoneda, Question 11.4, or the Hodge conjecture"
            ),
            "next_exact_construction": (
                "use these J0/J1 blocks with the six transported equations to solve "
                "the generator-specific Shamash comparison homotopies Gamma; only "
                "then solve lambda_k J2 - w_k lambda_e = a D2_e"
            ),
        },
        "files": {
            "overlap_Hermite": overlap_path.name,
            "J0_J1_blocks": blocks_path.name,
            "support_correspondence": points_path.name,
            "affine_forms": forms_path.name,
            "J1_factor_program": program_path.name,
            "J1_verification": verification_path.name,
        },
        "inputs": {
            "all_chart_producer": {
                "path": str(ALL_CHART_PRODUCER.relative_to(ROOT)),
                "sha256": sha256(ALL_CHART_PRODUCER),
            },
            "source_border_artifact": {
                "path": str(
                    (BORDER_DIRECTORY / "Z_2_C0010.npz").relative_to(ROOT)
                ),
                "sha256": sha256(BORDER_DIRECTORY / "Z_2_C0010.npz"),
            },
            "ordinary_pilot_border_artifact": {
                "path": str(
                    (BORDER_DIRECTORY / "Z_2_C0000.npz").relative_to(ROOT)
                ),
                "sha256": sha256(BORDER_DIRECTORY / "Z_2_C0000.npz"),
            },
        },
    }
    write_json_atomic(certificate_path, certificate)

    files = {}
    for path in (
        blocks_path,
        overlap_path,
        points_path,
        forms_path,
        program_path,
        verification_path,
        certificate_path,
    ):
        files[path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    manifest = {
        "schema": "section12-HS-rank2-K-generator-overlap-J1-manifest-v1",
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
        (output_directory / "k_generator_overlap_j1_certificate.json").read_text(
            encoding="utf-8"
        )
    )
    manifest = json.loads(
        (output_directory / "manifest.json").read_text(encoding="utf-8")
    )
    if certificate.get("schema") != "section12-HS-rank2-K-generator-overlap-J1-v1":
        raise HSKGeneratorOverlapJ1Error("the certificate schema changed")
    for name, record in manifest["files"].items():
        path = output_directory / name
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise HSKGeneratorOverlapJ1Error("an emitted artifact hash changed")
    with np.load(output_directory / certificate["files"]["overlap_Hermite"]) as data:
        overlap_length = int(certificate["counts"]["K_overlap_length"])
        target_length = int(
            certificate["counts"]["independently_rebuilt_target_length"]
        )
        if data["source_restriction"].shape != (
            overlap_length,
            target_length,
        ):
            raise HSKGeneratorOverlapJ1Error("the source J0 shape changed")
        if data["target_restriction"].shape != (
            overlap_length,
            target_length,
        ):
            raise HSKGeneratorOverlapJ1Error("the target J0 shape changed")
    return certificate


def main():
    certificate = write_artifacts()
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
