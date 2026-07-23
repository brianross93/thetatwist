"""Resolve the exact global-section part of the six-surface HS chart problem.

This producer does two things which were previously only stated as missing.

First, it fixes the common product atlas on ``E^4``.  The elliptic curve is
covered by ``Z != 0`` and ``Y != 0``; their fourfold product has exactly 16
charts and 32 oriented adjacent edges.  The producer emits both cube
incidence matrices and checks their exact ranks over ``F_31``.

Second, it constructs honest global sections which generate the two required
jets of every K-orbit double point.  In the complete degree-two tensor band

    H^0(E^4, O_E(6O) external-product ... external-product O_E(6O)),

there are ``6^4=1296`` coefficients.  For each surface the 81 lifted length-
two points impose 162 value/tangent equations.  The matrix has rank 161.
Adding the complementary surface derivative and the second tangent jet adds
two independent rows.  The resulting rank-163 system is solved twice to
produce sections ``X_s`` and ``Q_s`` with

    X_s|W=Q_s|W=0,
    dX_s(u)=1, D_v^2 X_s=0,
    dQ_s(u)=0, D_v^2 Q_s=1.

Thus ``X_s,Q_s`` generate the completed ideal ``(x_s,y_s^2)`` at the chosen
lift.  Their K translates give the corresponding generators at all 81 lifts.

This is an actual global-section/jet result.  It is not the final HS Cech
cocycle.  To obtain that cocycle one must still choose chart-local syzygy
bases for the full K-orbit ideal, restrict them on the 32 edges, and solve the
resulting sheaf-valued cube equations.  The scalar cube matrices emitted here
are the exact topological skeleton of that solve; they are not substituted
for the absent rational-function coefficient blocks.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (
    cm_bridge_six_hs_local_seed_emitter as seed_emitter,
)


P = 31
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_global_hs_chart_resolution"
)

SEED_CERTIFICATE = ROOT / (
    "results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/"
    "six_hs_local_seed_certificate.json"
)
SEEDS = ROOT / (
    "results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/"
    "six_hs_local_seeds.csv"
)
SUPPORT_ORBITS = ROOT / (
    "results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/"
    "six_hs_support_K_orbits.csv"
)
SECTION_COEFFICIENTS = ROOT / (
    "results/section12_instantiation/independent_incidence_fiber/"
    "candidate_section_coefficients.csv"
)
INCIDENCE_CERTIFICATE = ROOT / (
    "results/section12_instantiation/independent_incidence_fiber/"
    "independent_incidence_certificate.json"
)

FACTOR_QUADRATIC_BASIS = tuple(itertools.combinations_with_replacement(range(3), 2))
TENSOR_BASIS = tuple(itertools.product(range(6), repeat=4))


class GlobalHSChartResolutionError(RuntimeError):
    """Raised when an exact upstream contract or matrix invariant changes."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def dot(left, right):
    return sum(int(a) * int(b) for a, b in zip(left, right)) % P


def rank_mod(matrix):
    work = [[int(value) % P for value in row] for row in matrix]
    rank = 0
    column_count = len(work[0]) if work else 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, P)
        work[rank] = [value * inverse % P for value in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                (left - multiple * right) % P
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
    return rank


def sparse_rank_solve(rows, rhs_columns, variable_count):
    """Echelonize an underdetermined system and set all free variables to zero."""

    if len(rows) != len(rhs_columns):
        raise ValueError("the coefficient and right-hand-side row counts differ")
    rhs_count = len(rhs_columns[0]) if rhs_columns else 0
    pivots = []
    pivot_rows = {}
    for coefficients, rhs in zip(rows, rhs_columns):
        row = [int(value) % P for value in coefficients]
        out = [int(value) % P for value in rhs]
        for pivot in pivots:
            multiple = row[pivot]
            if not multiple:
                continue
            known, known_rhs = pivot_rows[pivot]
            row = [
                (left - multiple * right) % P
                for left, right in zip(row, known)
            ]
            out = [
                (left - multiple * right) % P
                for left, right in zip(out, known_rhs)
            ]
        pivot = next((index for index, value in enumerate(row) if value), None)
        if pivot is None:
            if any(out):
                raise GlobalHSChartResolutionError("the normalized jet system is inconsistent")
            continue
        inverse = pow(row[pivot], -1, P)
        row = [value * inverse % P for value in row]
        out = [value * inverse % P for value in out]
        pivots.append(pivot)
        pivot_rows[pivot] = (row, out)

    solutions = [[0] * variable_count for _ in range(rhs_count)]
    for pivot in reversed(pivots):
        row, out = pivot_rows[pivot]
        for rhs_index in range(rhs_count):
            correction = sum(
                row[column] * solutions[rhs_index][column]
                for column in range(pivot + 1, variable_count)
            )
            solutions[rhs_index][pivot] = (out[rhs_index] - correction) % P

    for coefficients, rhs in zip(rows, rhs_columns):
        for rhs_index, solution in enumerate(solutions):
            if dot(coefficients, solution) != int(rhs[rhs_index]) % P:
                raise GlobalHSChartResolutionError("a canonical jet solution failed substitution")
    return pivots, solutions


def load_inputs():
    seed_certificate = load_json(SEED_CERTIFICATE)
    incidence = load_json(INCIDENCE_CERTIFICATE)
    if seed_certificate.get("status") != (
        "cyclic-seed-transport-closed-for-independent-fibre;"
        "six-rational-local-HS-seeds-emitted;Z2-seed-preserved;"
        "six-completed-local-cocycles-emitted;global-Yoneda-numerators-open"
    ):
        raise GlobalHSChartResolutionError("the six-surface seed contract changed")
    if incidence.get("status") != "independent-incidence-fiber-certified":
        raise GlobalHSChartResolutionError("the independent incidence fibre changed")
    if not incidence.get("projective_chart_ledger", {}).get(
        "all_six_fivefold_intersections_empty"
    ):
        # Older certificates place the same exact result in exact_geometry.
        fivefold = incidence.get("exact_geometry", {}).get("fivefold_intersections", [])
        if len(fivefold) != 6 or not all(row.get("all_16_charts_unit") for row in fivefold):
            raise GlobalHSChartResolutionError("the six theta opens are not certified to cover")

    seeds = read_csv(SEEDS)
    support_rows = read_csv(SUPPORT_ORBITS)
    if len(seeds) != 6 or len(support_rows) != 486:
        raise GlobalHSChartResolutionError("the six K-orbit point schemes are incomplete")

    section_coefficients = {label: [0] * 81 for label in range(6)}
    for row in read_csv(SECTION_COEFFICIENTS):
        section_coefficients[int(row["label"])][int(row["theta_index"], 3)] = (
            int(row["coefficient_mod_31"]) % P
        )
    if any(not any(row) for row in section_coefficients.values()):
        raise GlobalHSChartResolutionError("a theta section coefficient row is empty")
    return seed_certificate, incidence, seeds, support_rows, section_coefficients


def factor_second_derivatives(points):
    return [
        (0, 0, 0)
        if point is None
        else (0, 6 * point[0] * point[0] % P, 12 * point[0] * point[1] % P)
        for point in points
    ]


def factor_quadratic_jets(point_index, direction, values, first, second):
    value_rows = []
    first_rows = []
    second_rows = []
    for left, right in FACTOR_QUADRATIC_BASIS:
        left_value = values[point_index][left]
        right_value = values[point_index][right]
        left_first = direction * first[point_index][left] % P
        right_first = direction * first[point_index][right] % P
        left_second = direction * direction * second[point_index][left] % P
        right_second = direction * direction * second[point_index][right] % P
        value_rows.append(left_value * right_value % P)
        first_rows.append(
            (left_first * right_value + left_value * right_first) % P
        )
        second_rows.append(
            (
                left_second * right_value
                + 2 * left_first * right_first
                + left_value * right_second
            )
            % P
        )
    return value_rows, first_rows, second_rows


def tensor_degree_two_jets(point, direction, values, first, second):
    factors = [
        factor_quadratic_jets(point[index], direction[index], values, first, second)
        for index in range(4)
    ]
    value_row = []
    first_row = []
    second_row = []
    for indices in TENSOR_BASIS:
        factor_values = [factors[index][0][indices[index]] for index in range(4)]
        factor_first = [factors[index][1][indices[index]] for index in range(4)]
        factor_second = [factors[index][2][indices[index]] for index in range(4)]

        value = 1
        for entry in factor_values:
            value = value * entry % P

        derivative = 0
        second_derivative = 0
        for index in range(4):
            term = factor_first[index]
            second_term = factor_second[index]
            for other in range(4):
                if other != index:
                    term = term * factor_values[other] % P
                    second_term = second_term * factor_values[other] % P
            derivative += term
            second_derivative += second_term
        for left in range(4):
            for right in range(left + 1, 4):
                term = 2 * factor_first[left] * factor_first[right]
                for other in range(4):
                    if other not in (left, right):
                        term = term * factor_values[other] % P
                second_derivative += term
        value_row.append(value)
        first_row.append(derivative % P)
        second_row.append(second_derivative % P)
    return value_row, first_row, second_row


def section_gradient(coefficients, point, values, first):
    return seed_emitter.section_gradient(coefficients, point, values, first)


def complementary_surface_vector(seed, point, coefficients, values, first):
    labels = json.loads(seed["surface_section_labels"])
    x_row = json.loads(seed["x_fixed_frame_coefficients"])
    y_row = json.loads(seed["y_fixed_frame_coefficients"])
    normals = [
        section_gradient(coefficients[label], point, values, first) for label in labels
    ]
    vector = seed_emitter.solve4([*normals, x_row, y_row], [0, 0, 1, 0])
    tangent = json.loads(seed["tangent_fixed_frame_coefficients"])
    if any(dot(normal, vector) for normal in normals):
        raise GlobalHSChartResolutionError("the complementary vector left the surface")
    if dot(x_row, vector) != 1 or dot(y_row, vector):
        raise GlobalHSChartResolutionError("the complementary vector lost normalization")
    if any(dot(normal, tangent) for normal in normals):
        raise GlobalHSChartResolutionError("the recorded tangent left the surface")
    return vector


def cube_atlas():
    charts = list(itertools.product((0, 1), repeat=4))
    chart_ids = {bits: "C" + "".join(map(str, bits)) for bits in charts}
    chart_rows = []
    for bits in charts:
        chart_rows.append(
            {
                "chart": chart_ids[bits],
                "bits_Z0_Y1": "".join(map(str, bits)),
                "nonvanishing_coordinates": ",".join(
                    "Y" if bit else "Z" for bit in bits
                ),
                "H_frame": "*".join(
                    f"{'Y' if bit else 'Z'}_{index}" for index, bit in enumerate(bits)
                ),
            }
        )

    edge_rows = []
    edge_lookup = {}
    for bits in charts:
        for factor in range(4):
            if bits[factor]:
                continue
            target = list(bits)
            target[factor] = 1
            target = tuple(target)
            edge_id = len(edge_rows)
            edge_lookup[(bits, factor)] = edge_id
            edge_rows.append(
                {
                    "edge_id": edge_id,
                    "source_chart": chart_ids[bits],
                    "target_chart": chart_ids[target],
                    "changed_factor": factor,
                    "H_frame_transition_source_to_target": f"Y_{factor}/Z_{factor}",
                    "2H_frame_transition_source_to_target": f"(Y_{factor}/Z_{factor})^2",
                    "4H_frame_transition_source_to_target": f"(Y_{factor}/Z_{factor})^4",
                }
            )

    d0 = [[0] * len(charts) for _ in edge_rows]
    for row in edge_rows:
        source = tuple(map(int, row["source_chart"][1:]))
        target = tuple(map(int, row["target_chart"][1:]))
        d0[int(row["edge_id"])][charts.index(source)] = P - 1
        d0[int(row["edge_id"])][charts.index(target)] = 1

    face_rows = []
    d1 = []
    for left, right in itertools.combinations(range(4), 2):
        other = [index for index in range(4) if index not in (left, right)]
        for other_bits in itertools.product((0, 1), repeat=2):
            base = [0] * 4
            for index, bit in zip(other, other_bits):
                base[index] = bit
            base = tuple(base)
            left_vertex = list(base)
            left_vertex[left] = 1
            left_vertex = tuple(left_vertex)
            right_vertex = list(base)
            right_vertex[right] = 1
            right_vertex = tuple(right_vertex)
            row = [0] * len(edge_rows)
            row[edge_lookup[(base, left)]] = 1
            row[edge_lookup[(left_vertex, right)]] = 1
            row[edge_lookup[(right_vertex, left)]] = P - 1
            row[edge_lookup[(base, right)]] = P - 1
            d1.append(row)
            face_rows.append(
                {
                    "face_id": len(face_rows),
                    "varying_factors": f"{left},{right}",
                    "fixed_bits": ",".join(
                        f"{index}:{base[index]}" for index in other
                    ),
                }
            )

    if len(chart_rows) != 16 or len(edge_rows) != 32 or len(face_rows) != 24:
        raise ArithmeticError("the four-cube atlas count changed")
    if rank_mod(d0) != 15 or rank_mod(d1) != 17:
        raise ArithmeticError("a four-cube Cech rank changed")
    for face in d1:
        for column in range(16):
            if sum(face[edge] * d0[edge][column] for edge in range(32)) % P:
                raise ArithmeticError("the cube boundary no longer squares to zero")
    return chart_rows, edge_rows, face_rows, d0, d1


def matrix_coo(matrix, row_label, column_label):
    rows = []
    for row_index, row in enumerate(matrix):
        for column_index, value in enumerate(row):
            if value % P:
                rows.append(
                    {
                        row_label: row_index,
                        column_label: column_index,
                        "value_mod_31": value % P,
                    }
                )
    return rows


def canonical_chart(projective_points):
    bits = []
    for x, y, z in projective_points:
        if z:
            bits.append(0)
        elif y:
            bits.append(1)
        else:
            raise GlobalHSChartResolutionError("a projective point misses both atlas charts")
    return "C" + "".join(map(str, bits))


def write_artifacts(output_directory: Path = OUTPUT_DIRECTORY):
    seed_certificate, incidence, seeds, support_rows, section_coefficients = load_inputs()
    points, point_values, point_first = seed_emitter.curve_data()
    point_second = factor_second_derivatives(points)
    chart_rows, edge_rows, face_rows, d0, d1 = cube_atlas()

    supports_by_surface = {f"Z_{index}": [] for index in range(6)}
    assignment_rows = []
    for row in support_rows:
        supports_by_surface[row["surface"]].append(row)
        projective = json.loads(row["source_projective_points_X_Y_Z"])
        assignment_rows.append(
            {
                "surface": row["surface"],
                "K_orbit_lift_index": row["K_orbit_lift_index"],
                "source_point_indices": row["source_point_indices"],
                "canonical_product_chart": canonical_chart(projective),
            }
        )

    basis_rows = []
    for basis_id, indices in enumerate(TENSOR_BASIS):
        basis_rows.append(
            {
                "basis_id": basis_id,
                "factor_quadratic_basis_indices": json.dumps(list(indices)),
                "factor_monomials": "*".join(
                    f"b{FACTOR_QUADRATIC_BASIS[index][0]}b{FACTOR_QUADRATIC_BASIS[index][1]}"
                    for index in indices
                ),
            }
        )

    summary_rows = []
    generator_rows = []
    normalization_rows = []
    system_fingerprints = []
    for surface_index, seed in enumerate(seeds):
        surface = f"Z_{surface_index}"
        tangent = json.loads(seed["tangent_fixed_frame_coefficients"])
        representative = json.loads(seed["source_point_indices"])
        complement = complementary_surface_vector(
            seed, representative, section_coefficients, point_values, point_first
        )

        equations = []
        for support in supports_by_surface[surface]:
            point = json.loads(support["source_point_indices"])
            value_row, tangent_row, _ = tensor_degree_two_jets(
                point, tangent, point_values, point_first, point_second
            )
            equations.extend((value_row, tangent_row))
        if len(equations) != 162:
            raise GlobalHSChartResolutionError("a K-orbit jet system is incomplete")

        _, transverse_row, _ = tensor_degree_two_jets(
            representative, complement, point_values, point_first, point_second
        )
        _, _, tangent_second_row = tensor_degree_two_jets(
            representative, tangent, point_values, point_first, point_second
        )
        rhs = [[0, 0] for _ in equations] + [[1, 0], [0, 1]]
        full_system = [*equations, transverse_row, tangent_second_row]
        pivots, solutions = sparse_rank_solve(full_system, rhs, len(TENSOR_BASIS))

        base_rank = rank_mod(equations)
        full_rank = len(pivots)
        if base_rank != 161 or full_rank != 163:
            raise GlobalHSChartResolutionError("the degree-two global jet ranks changed")
        if rank_mod([*equations, transverse_row]) != 162:
            raise GlobalHSChartResolutionError("the transverse jet is no longer independent")
        if rank_mod([*equations, tangent_second_row]) != 162:
            raise GlobalHSChartResolutionError("the quadratic tangent jet is no longer independent")

        for generator, solution in zip(("X", "Q"), solutions):
            for basis_id, coefficient in enumerate(solution):
                if coefficient:
                    generator_rows.append(
                        {
                            "surface": surface,
                            "generator": generator,
                            "basis_id": basis_id,
                            "coefficient_mod_31": coefficient,
                        }
                    )
        normalization_rows.extend(
            [
                {
                    "surface": surface,
                    "generator": "X",
                    "vanishes_on_all_81_length_two_lifts": True,
                    "representative_transverse_derivative_mod_31": 1,
                    "representative_second_tangent_derivative_mod_31": 0,
                },
                {
                    "surface": surface,
                    "generator": "Q",
                    "vanishes_on_all_81_length_two_lifts": True,
                    "representative_transverse_derivative_mod_31": 0,
                    "representative_second_tangent_derivative_mod_31": 1,
                },
            ]
        )
        summary_rows.append(
            {
                "surface": surface,
                "unknown_coefficients": len(TENSOR_BASIS),
                "K_orbit_lifts": 81,
                "value_plus_tangent_equations": len(equations),
                "value_plus_tangent_rank_mod_31": base_rank,
                "value_plus_tangent_kernel_dimension": len(TENSOR_BASIS) - base_rank,
                "normalization_equations": 2,
                "full_equation_count": len(full_system),
                "full_rank_mod_31": full_rank,
                "full_solution_dimension": len(TENSOR_BASIS) - full_rank,
                "global_generators_solved": True,
                "final_Cech_cocycle_solved": False,
            }
        )
        encoded = bytearray()
        for row in equations:
            encoded.extend(bytes(row))
        system_fingerprints.append(
            {
                "surface": surface,
                "jet_matrix_sha256": hashlib.sha256(encoded).hexdigest(),
                "complementary_surface_vector": complement,
                "representative_tangent": tangent,
            }
        )

    output_directory.mkdir(parents=True, exist_ok=True)
    paths = {}
    artifacts = [
        ("product_theta_atlas.csv", chart_rows),
        ("product_theta_atlas_edges.csv", edge_rows),
        ("product_theta_atlas_faces.csv", face_rows),
        ("cech_d0_vertices_to_edges.coo.csv", matrix_coo(d0, "edge_id", "chart_index")),
        ("cech_d1_edges_to_faces.coo.csv", matrix_coo(d1, "face_id", "edge_id")),
        ("support_chart_assignment.csv", assignment_rows),
        ("degree2_tensor_basis.csv", basis_rows),
        ("degree2_global_generator_coefficients.coo.csv", generator_rows),
        ("degree2_global_generator_normalization.csv", normalization_rows),
        ("degree2_global_jet_system_summary.csv", summary_rows),
    ]
    for name, rows in artifacts:
        path = output_directory / name
        write_csv(path, rows, list(rows[0]))
        paths[name] = path

    edge_contract = {
        "schema": "section12-global-hs-chart-edge-system-v1",
        "target_output": "global_hs_theta_chart_cocycle.coo.csv",
        "common_atlas": {
            "charts": 16,
            "adjacent_edges": 32,
            "square_faces": 24,
            "rank_d0_mod_31": 15,
            "rank_d1_mod_31": 17,
            "scalar_ker_d1_dimension": 15,
            "scalar_cube_H1_dimension": 0,
        },
        "actual_global_section_input": {
            "surface_sectors": 6,
            "degree2_unknowns_per_generator": 1296,
            "jet_equations_per_surface": 162,
            "jet_rank_per_surface": 161,
            "normalized_equations_per_surface": 164,
            "normalized_rank_per_surface": 163,
            "generators": ["X_s", "Q_s"],
            "K_transport": "use all 81 pullbacks of each emitted generator",
        },
        "smallest_unsolved_sheaf_valued_system": {
            "known_topological_rows": "32 edge restrictions and 24 square-face equations",
            "known_topological_ranks": "rank(d0)=15 and rank(d1)=17 per scalar coefficient channel",
            "unknown_coefficient_blocks": [
                "chart-local syzygy bases of the K-invariant ideal generated by the 81 translates of X_s and Q_s",
                "their restriction matrices on the 32 adjacent product-chart edges",
                "the K-linearization matrices on those chart-local syzygy bases",
            ],
            "why_no_total_unknown_count_is_printed": (
                "the rank of each rational-function coefficient block is itself an output "
                "of the missing chart-local syzygy computation; multiplying the scalar cube "
                "counts by an invented block size would fabricate a finite system"
            ),
            "solve_after_blocks_exist": (
                "assemble the sheaf-valued d0 and d1 by replacing each nonzero scalar cube "
                "entry with the corresponding chart restriction block, impose d1*g=0 and "
                "K-equivariance, quotient by d0 gauge, and emit the surviving edge matrices"
            ),
        },
        "scope": (
            "This file is the exact finite atlas skeleton and global-section input. "
            "It is not the final HS Cech or Yoneda cocycle."
        ),
    }
    edge_contract_path = output_directory / "global_hs_theta_chart_cocycle_linear_system.json"
    edge_contract_path.write_text(
        json.dumps(edge_contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    paths[edge_contract_path.name] = edge_contract_path

    fingerprints_path = output_directory / "global_jet_matrix_fingerprints.json"
    fingerprints_path.write_text(
        json.dumps(system_fingerprints, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    paths[fingerprints_path.name] = fingerprints_path

    certificate = {
        "schema": "section12-global-hs-chart-resolution-v1",
        "status": (
            "common-16-chart-32-edge-atlas-certified;"
            "six-degree2-global-HS-jet-pairs-solved;"
            "chart-local-syzygy-edge-blocks-open;final-Cech-and-Yoneda-cocycles-open"
        ),
        "field": "F_31",
        "common_product_atlas": edge_contract["common_atlas"],
        "degree2_global_jet_result": {
            "line_bundle_upstairs": "O_E(6O) external-product^4",
            "basis_dimension": 1296,
            "surface_count": 6,
            "all_value_tangent_ranks": [int(row["value_plus_tangent_rank_mod_31"]) for row in summary_rows],
            "all_normalized_ranks": [int(row["full_rank_mod_31"]) for row in summary_rows],
            "actual_generators_emitted": True,
            "interpretation": (
                "X_s and Q_s generate the x and y^2 jets at the chosen lift; "
                "their K pullbacks do so at all 81 lifts"
            ),
        },
        "global_output": {
            "global_sections_and_jets_computed": True,
            "global_hs_theta_chart_cocycle_computed": False,
            "global_Yoneda_numerators_computed": False,
            "T_computed": False,
        },
        "exact_remaining_input": edge_contract["smallest_unsolved_sheaf_valued_system"],
        "scope": {
            "proved": [
                "the common product atlas has 16 charts, 32 adjacent edges, and 24 square faces",
                "the scalar atlas incidence ranks are 15 and 17",
                "each K-orbit double point has a 162 by 1296 degree-two jet matrix of rank 161",
                "two independent normalization jets raise the rank to 163",
                "explicit global degree-two sections X_s and Q_s with the required local jets exist for all six surfaces",
            ],
            "not_claimed": [
                "that two global sections cut only W_s",
                "a chart-local syzygy presentation of the full K-invariant ideal",
                "a solved sheaf-valued 32-edge Cech system",
                "a global HS or off-diagonal Yoneda cocycle",
                "a value for T or any Hodge-conjecture claim",
            ],
        },
        "provenance": {
            "six_seed_certificate": str(SEED_CERTIFICATE.relative_to(ROOT)),
            "six_seed_certificate_sha256": sha256(SEED_CERTIFICATE),
            "independent_incidence_certificate": str(INCIDENCE_CERTIFICATE.relative_to(ROOT)),
            "independent_incidence_certificate_sha256": sha256(INCIDENCE_CERTIFICATE),
        },
    }
    certificate_path = output_directory / "global_hs_chart_resolution_certificate.json"
    certificate_path.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    paths[certificate_path.name] = certificate_path

    manifest = {
        "schema": "section12-global-hs-chart-resolution-manifest-v1",
        "generator": str(Path(__file__).relative_to(ROOT)),
        "generator_sha256": sha256(Path(__file__)),
        "files": {
            name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for name, path in sorted(paths.items())
        },
    }
    manifest_path = output_directory / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return certificate


if __name__ == "__main__":
    result = write_artifacts()
    print(json.dumps(result, indent=2, sort_keys=True))
