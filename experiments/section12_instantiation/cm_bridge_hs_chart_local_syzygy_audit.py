"""Compute the exact chart-local *symbol* syzygies of the six HS ideals.

The global HS chart producer emits two normalized sections ``X_s,Q_s`` in
``H0(E^4,M^2)`` for each of the six surfaces.  The cached theta-group data
also contain the projective action of every three-torsion translation on
``H0(E,O_E(3O))``.  This producer induces that action on
``H0(E,O_E(6O))`` and hence reconstructs the 81 K-translates of each of
``X_s,Q_s`` without inventing a translation matrix.

For every support lift, reduce those 162 generators modulo the maximal ideal
times ``(x_s,y_s^2)``.  The result is the two-dimensional conormal-symbol
map.  On each of the 16 product charts and each of the 32 adjacent overlaps,
this is an honest finite matrix over F_31.  The module emits canonical
nullspace bases and exact restriction matrices between them.

This is deliberately not called the full sheaf syzygy module.  The latter is
the kernel over the affine surface coordinate ring, not the kernel after
passing to the residue symbols at the finite support.  It additionally needs
an ideal saturation/base-locus calculation and polynomial-module Groebner
bases.  A second exact obstruction is recorded here: the 16-chart product
atlas is not K-stable.  Its K-orbit consists of 1296 distinct charts, so a
K-linearized Cech solve cannot use only the original 32 edges.
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
    cm_bridge_global_hs_chart_resolution as global_hs,
)
from experiments.section12_instantiation import (
    cm_bridge_six_hs_local_seed_emitter as seed_emitter,
)


P = 31
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit"
)
GLOBAL_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_global_hs_chart_resolution"
)
GLOBAL_CERTIFICATE = GLOBAL_DIRECTORY / "global_hs_chart_resolution_certificate.json"
GLOBAL_GENERATORS = GLOBAL_DIRECTORY / "degree2_global_generator_coefficients.coo.csv"
ATLAS = GLOBAL_DIRECTORY / "product_theta_atlas.csv"
ATLAS_EDGES = GLOBAL_DIRECTORY / "product_theta_atlas_edges.csv"
SEEDS = ROOT / (
    "results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/"
    "six_hs_local_seeds.csv"
)
SUPPORTS = ROOT / (
    "results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/"
    "six_hs_support_K_orbits.csv"
)
DESCENT = ROOT / (
    "results/section12_instantiation/descent_section/descent_section_certificate.json"
)
THETA_COEFFICIENTS = ROOT / (
    "results/section12_instantiation/independent_incidence_fiber/"
    "candidate_section_coefficients.csv"
)


class ChartLocalSyzygyAuditError(RuntimeError):
    """A bound exact input or a finite-field invariant changed."""


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


def matrix_multiply(left, right):
    rows = len(left)
    middle = len(right)
    columns = len(right[0])
    return [
        [
            sum(left[row][index] * right[index][column] for index in range(middle))
            % P
            for column in range(columns)
        ]
        for row in range(rows)
    ]


def matrix_power(matrix, exponent):
    result = [[int(row == column) for column in range(len(matrix))] for row in range(len(matrix))]
    factor = matrix
    while exponent:
        if exponent & 1:
            result = matrix_multiply(result, factor)
        factor = matrix_multiply(factor, factor)
        exponent //= 2
    return result


def matrix_vector(matrix, vector):
    return [dot(row, vector) for row in matrix]


def canonical_projective(vector):
    vector = [int(value) % P for value in vector]
    pivot = next((value for value in vector if value), None)
    if pivot is None:
        raise ArithmeticError("the zero vector has no projective class")
    inverse = pow(pivot, -1, P)
    return tuple(value * inverse % P for value in vector)


def projective_scalar(left, right):
    """Return the unique nonzero scalar with left=scalar*right."""

    if right and isinstance(right[0], list):
        left = [value for row in left for value in row]
        right = [value for row in right for value in row]
    pivot = next(index for index, value in enumerate(right) if value % P)
    scalar = left[pivot] * pow(right[pivot], -1, P) % P
    if not scalar or any((a - scalar * b) % P for a, b in zip(left, right)):
        raise ChartLocalSyzygyAuditError("two asserted projective vectors disagree")
    return scalar


def rref_nullspace(matrix, column_count):
    """Return RREF pivots, free columns, and the canonical nullspace basis."""

    work = [[int(value) % P for value in row] for row in matrix if any(value % P for value in row)]
    pivot_columns = []
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, P)
        work[pivot_row] = [value * inverse % P for value in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                (left - multiple * right) % P
                for left, right in zip(work[row], work[pivot_row])
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(work):
            break
    work = work[:pivot_row]
    free_columns = [column for column in range(column_count) if column not in set(pivot_columns)]
    basis = []
    for free in free_columns:
        vector = [0] * column_count
        vector[free] = 1
        for row, pivot in enumerate(pivot_columns):
            vector[pivot] = -work[row][free] % P
        if any(dot(source, vector) for source in matrix):
            raise ArithmeticError("a canonical nullspace vector failed substitution")
        basis.append(vector)
    return pivot_columns, free_columns, basis


def coordinates_in_nullspace(vector, free_columns, basis):
    coefficients = [vector[column] % P for column in free_columns]
    reconstructed = [0] * len(vector)
    for coefficient, source in zip(coefficients, basis):
        if not coefficient:
            continue
        reconstructed = [
            (left + coefficient * right) % P
            for left, right in zip(reconstructed, source)
        ]
    if reconstructed != [value % P for value in vector]:
        raise ArithmeticError("a restriction vector left the target nullspace basis")
    return coefficients


def load_inputs():
    global_certificate = load_json(GLOBAL_CERTIFICATE)
    descent = load_json(DESCENT)
    if global_certificate.get("status") != (
        "common-16-chart-32-edge-atlas-certified;"
        "six-degree2-global-HS-jet-pairs-solved;"
        "chart-local-syzygy-edge-blocks-open;final-Cech-and-Yoneda-cocycles-open"
    ):
        raise ChartLocalSyzygyAuditError("the global HS chart contract changed")
    if descent.get("status") != "certified":
        raise ChartLocalSyzygyAuditError("the exact K translation contract changed")

    seeds = read_csv(SEEDS)
    supports = read_csv(SUPPORTS)
    atlas = read_csv(ATLAS)
    edges = read_csv(ATLAS_EDGES)
    if len(seeds) != 6 or len(supports) != 486:
        raise ChartLocalSyzygyAuditError("the six support orbits are incomplete")
    if len(atlas) != 16 or len(edges) != 32:
        raise ChartLocalSyzygyAuditError("the product atlas changed")

    theta = {label: [0] * 81 for label in range(6)}
    for row in read_csv(THETA_COEFFICIENTS):
        theta[int(row["label"])][int(row["theta_index"], 3)] = int(
            row["coefficient_mod_31"]
        ) % P

    generators = {(surface, generator): [0] * 1296 for surface in range(6) for generator in ("X", "Q")}
    for row in read_csv(GLOBAL_GENERATORS):
        generators[(int(row["surface"].split("_")[1]), row["generator"])][
            int(row["basis_id"])
        ] = int(row["coefficient_mod_31"]) % P
    if any(not any(vector) for vector in generators.values()):
        raise ChartLocalSyzygyAuditError("a normalized global generator is empty")
    return global_certificate, descent, seeds, supports, atlas, edges, theta, generators


def group_data(descent, points):
    translation_table = {
        tuple(map(int, key.split(","))): [[int(value) % P for value in row] for row in matrix]
        for key, matrix in descent["translation_matrices"].items()
    }
    torsion_p = (0, 1)
    torsion_q = (3, 11)
    torsion = {
        (a, b): seed_emitter.add_points(
            seed_emitter.multiply_point(a, torsion_p),
            seed_emitter.multiply_point(b, torsion_q),
        )
        for a in range(3)
        for b in range(3)
    }
    rows = descent["kernel_rows"]
    eigenvalues = [int(value) % P for value in descent["eigenline"]["eigenvalues"]]
    elements = list(itertools.product(range(3), repeat=4))
    element_index = {element: index for index, element in enumerate(elements)}
    records = []
    for exponents in elements:
        flat = [
            sum(exponents[index] * rows[index][column] for index in range(4)) % 3
            for column in range(8)
        ]
        keys = tuple((flat[2 * factor], flat[2 * factor + 1]) for factor in range(4))
        shifts = tuple(torsion[key] for key in keys)
        factor_matrices = []
        for factor in range(4):
            matrix = [[int(row == column) for column in range(3)] for row in range(3)]
            for generator, exponent in enumerate(exponents):
                key = (rows[generator][2 * factor], rows[generator][2 * factor + 1])
                matrix = matrix_multiply(matrix, matrix_power(translation_table[key], exponent))
            factor_matrices.append(matrix)
        normalization = 1
        for exponent, eigenvalue in zip(exponents, eigenvalues):
            normalization = normalization * pow(eigenvalue, -2 * exponent, P) % P
        records.append(
            {
                "exponents": exponents,
                "keys": keys,
                "shifts": shifts,
                "factor_matrices": factor_matrices,
                "normalization": normalization,
            }
        )

    # The normalized degree-two operators form an honest K representation.
    group_checks = 0
    for left in records:
        for right in records:
            target_exp = tuple((a + b) % 3 for a, b in zip(left["exponents"], right["exponents"]))
            target = records[element_index[target_exp]]
            scalar = left["normalization"] * right["normalization"] % P
            for factor in range(4):
                product = matrix_multiply(left["factor_matrices"][factor], right["factor_matrices"][factor])
                factor_scalar = projective_scalar(product, target["factor_matrices"][factor])
                scalar = scalar * factor_scalar * factor_scalar % P
            if scalar != target["normalization"]:
                raise ChartLocalSyzygyAuditError("the normalized degree-two K action is projective")
            group_checks += 1
    return records, element_index, group_checks


def action_scale(record, point, point_values, point_index):
    scale = record["normalization"]
    translated = []
    for factor in range(4):
        source = list(point_values[point[factor]])
        image = matrix_vector(record["factor_matrices"][factor], source)
        target_point = seed_emitter.add_points(point_index[point[factor]], record["shifts"][factor])
        target = list(point_values[point_index.index(target_point)])
        factor_scale = projective_scalar(image, target)
        scale = scale * factor_scale * factor_scale % P
        translated.append(point_index.index(target_point))
    return scale, tuple(translated)


def chart_membership(point, bits, point_values):
    return all(point_values[index][2 if bit else 0] for index, bit in zip(point, bits))


def edge_membership(point, source_bits, changed_factor, point_values):
    for factor, bit in enumerate(source_bits):
        if factor == changed_factor:
            if not point_values[point[factor]][0] or not point_values[point[factor]][2]:
                return False
        elif not point_values[point[factor]][2 if bit else 0]:
            return False
    return True


def line_pullback(matrix, line):
    return canonical_projective(
        [sum(line[row] * matrix[row][column] for row in range(3)) % P for column in range(3)]
    )


def quadratic_operator(matrix):
    """Induce coordinate substitution on the six quadratic monomials."""

    pairs = global_hs.FACTOR_QUADRATIC_BASIS
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    operator = [[0] * 6 for _ in range(6)]
    for source, (left, right) in enumerate(pairs):
        for first in range(3):
            for second in range(3):
                target = pair_index[tuple(sorted((first, second)))]
                operator[target][source] = (
                    operator[target][source]
                    + matrix[left][first] * matrix[right][second]
                ) % P
    return operator


def apply_factor_operator(vector, operator, factor):
    result = [0] * len(vector)
    stride = 6 ** (3 - factor)
    block = stride * 6
    for start in range(0, len(vector), block):
        for offset in range(stride):
            source_values = [vector[start + source * stride + offset] for source in range(6)]
            for target in range(6):
                result[start + target * stride + offset] = dot(
                    operator[target], source_values
                )
    return result


def apply_degree2_action(vector, record):
    result = list(vector)
    for factor, matrix in enumerate(record["factor_matrices"]):
        result = apply_factor_operator(result, quadratic_operator(matrix), factor)
    return [record["normalization"] * value % P for value in result]


def stable_atlas_rows(records):
    rows = []
    keys = set()
    for record in records:
        for bits in itertools.product((0, 1), repeat=4):
            forms = tuple(
                line_pullback(
                    record["factor_matrices"][factor],
                    (0, 0, 1) if bit else (1, 0, 0),
                )
                for factor, bit in enumerate(bits)
            )
            key = json.dumps(forms)
            keys.add(key)
            rows.append(
                {
                    "K_exponents": json.dumps(list(record["exponents"])),
                    "base_chart": "C" + "".join(map(str, bits)),
                    "pulled_back_factor_linear_forms_Z_X_Y": key,
                }
            )
    if len(rows) != 1296 or len(keys) != 1296:
        raise ChartLocalSyzygyAuditError("the K-stable chart refinement count changed")
    return rows


def surface_symbol_matrix(
    surface,
    seed,
    support_rows,
    theta,
    generators,
    records,
    element_index,
    points,
    point_values,
    point_first,
    point_second,
):
    point_lookup = {point: index for index, point in enumerate(points)}
    representative = tuple(json.loads(seed["source_point_indices"]))
    tangent = json.loads(seed["tangent_fixed_frame_coefficients"])
    ordered_points = []
    for record in records:
        ordered_points.append(
            tuple(
                point_lookup[seed_emitter.add_points(points[representative[factor]], record["shifts"][factor])]
                for factor in range(4)
            )
        )
    expected = {tuple(json.loads(row["source_point_indices"])) for row in support_rows}
    if set(ordered_points) != expected or len(set(ordered_points)) != 81:
        raise ChartLocalSyzygyAuditError("the ordered K action does not reproduce a support orbit")

    base_vectors = [generators[(surface, "X")], generators[(surface, "Q")]]
    base_germs = []
    germ_rows = []
    complements = []
    for point in ordered_points:
        complement = global_hs.complementary_surface_vector(
            seed, point, theta, point_values, point_first
        )
        complements.append(complement)
        value_row, tangent_row, _ = global_hs.tensor_degree_two_jets(
            point, tangent, point_values, point_first, point_second
        )
        _, transverse_row, _ = global_hs.tensor_degree_two_jets(
            point, complement, point_values, point_first, point_second
        )
        _, _, tangent_second_row = global_hs.tensor_degree_two_jets(
            point, tangent, point_values, point_first, point_second
        )
        for vector in base_vectors:
            if dot(value_row, vector) or dot(tangent_row, vector):
                raise ChartLocalSyzygyAuditError("a global generator does not vanish on W")
        base_germs.append(
            [
                [dot(transverse_row, vector) for vector in base_vectors],
                [dot(tangent_second_row, vector) for vector in base_vectors],
            ]
        )
        germ_rows.append((transverse_row, tangent_second_row))

    matrix = [[0] * 162 for _ in range(162)]
    automorphy_checks = 0
    for support_index, point in enumerate(ordered_points):
        for group_index, record in enumerate(records):
            scale, translated = action_scale(
                record, point, point_values, points
            )
            translated_index = ordered_points.index(translated)
            expected_exp = tuple(
                (left + right) % 3
                for left, right in zip(records[support_index]["exponents"], record["exponents"])
            )
            if translated_index != element_index[expected_exp]:
                raise ChartLocalSyzygyAuditError("support translation lost the K group law")
            for component in range(2):
                for generator in range(2):
                    matrix[2 * support_index + component][2 * group_index + generator] = (
                        scale * base_germs[translated_index][component][generator]
                    ) % P
            automorphy_checks += 1
    pivots, _, _ = rref_nullspace(matrix, 162)
    for support_index in range(81):
        local_pivots, _, _ = rref_nullspace(
            matrix[2 * support_index : 2 * support_index + 2], 162
        )
        if len(local_pivots) != 2:
            raise ChartLocalSyzygyAuditError(
                "the translate-generated ideal lost a local conormal direction"
            )

    # Compare the automorphy formula with direct coefficient transport for
    # each of the four K generators at every support.  This checks both the
    # quadratic action convention and the invariant-derivative scaling.
    direct_checks = 0
    unit_indices = [
        element_index[tuple(int(index == generator) for index in range(4))]
        for generator in range(4)
    ]
    for group_index in unit_indices:
        for generator, base_vector in enumerate(base_vectors):
            translated_vector = apply_degree2_action(base_vector, records[group_index])
            for support_index, (transverse_row, tangent_second_row) in enumerate(germ_rows):
                for component, row in enumerate((transverse_row, tangent_second_row)):
                    if dot(row, translated_vector) != matrix[2 * support_index + component][
                        2 * group_index + generator
                    ]:
                        raise ChartLocalSyzygyAuditError(
                            "direct degree-two transport disagrees with the automorphy formula"
                        )
                    direct_checks += 1
    return ordered_points, matrix, len(pivots), automorphy_checks, direct_checks


def emit_basis_rows(surface, location_type, location, basis):
    rows = []
    for basis_index, vector in enumerate(basis):
        for generator_column, value in enumerate(vector):
            if value:
                rows.append(
                    {
                        "surface": f"Z_{surface}",
                        "location_type": location_type,
                        "location": location,
                        "basis_index": basis_index,
                        "generator_column": generator_column,
                        "K_element_index": generator_column // 2,
                        "generator": "X" if generator_column % 2 == 0 else "Q",
                        "value_mod_31": value,
                    }
                )
    return rows


def write_artifacts(output_directory: Path = OUTPUT_DIRECTORY):
    (
        global_certificate,
        descent,
        seeds,
        support_rows,
        atlas,
        edges,
        theta,
        generators,
    ) = load_inputs()
    points, point_values, point_first = seed_emitter.curve_data()
    point_second = global_hs.factor_second_derivatives(points)
    records, element_index, group_law_checks = group_data(descent, points)
    stable_rows = stable_atlas_rows(records)

    group_rows = []
    for index, record in enumerate(records):
        group_rows.append(
            {
                "K_element_index": index,
                "K_generator_exponents": json.dumps(list(record["exponents"])),
                "factor_E3_keys": json.dumps(record["keys"]),
                "degree2_linearization_scalar_mod_31": record["normalization"],
            }
        )

    chart_summary = []
    edge_summary = []
    chart_basis_rows = []
    edge_basis_rows = []
    restriction_rows = []
    surface_summary = []
    support_by_surface = {
        surface: [row for row in support_rows if row["surface"] == f"Z_{surface}"]
        for surface in range(6)
    }
    chart_bits = {
        row["chart"]: tuple(map(int, row["bits_Z0_Y1"])) for row in atlas
    }

    for surface, seed in enumerate(seeds):
        ordered_points, symbol, full_rank, automorphy_checks, direct_checks = surface_symbol_matrix(
            surface,
            seed,
            support_by_surface[surface],
            theta,
            generators,
            records,
            element_index,
            points,
            point_values,
            point_first,
            point_second,
        )
        surface_summary.append(
            {
                "surface": f"Z_{surface}",
                "translated_generators": 162,
                "support_lifts": 81,
                "conormal_symbol_shape": "162x162",
                "conormal_symbol_rank_mod_31": full_rank,
                "every_2x162_support_block_rank": 2,
                "automorphy_block_checks": automorphy_checks,
                "direct_quadratic_transport_checks": direct_checks,
            }
        )

        chart_data = {}
        for chart, bits in chart_bits.items():
            support_indices = [
                index
                for index, point in enumerate(ordered_points)
                if chart_membership(point, bits, point_values)
            ]
            rows = [symbol[2 * index + component] for index in support_indices for component in range(2)]
            pivots, free, basis = rref_nullspace(rows, 162)
            chart_data[chart] = (free, basis)
            chart_summary.append(
                {
                    "surface": f"Z_{surface}",
                    "chart": chart,
                    "support_lifts_in_chart": len(support_indices),
                    "symbol_rows": len(rows),
                    "symbol_rank_mod_31": len(pivots),
                    "constant_symbol_syzygy_dimension": len(basis),
                }
            )
            chart_basis_rows.extend(emit_basis_rows(surface, "chart", chart, basis))

        for edge in edges:
            edge_id = int(edge["edge_id"])
            source = edge["source_chart"]
            target = edge["target_chart"]
            changed = int(edge["changed_factor"])
            source_bits = chart_bits[source]
            support_indices = [
                index
                for index, point in enumerate(ordered_points)
                if edge_membership(point, source_bits, changed, point_values)
            ]
            rows = [symbol[2 * index + component] for index in support_indices for component in range(2)]
            pivots, free, basis = rref_nullspace(rows, 162)
            location = f"E{edge_id}"
            edge_summary.append(
                {
                    "surface": f"Z_{surface}",
                    "edge_id": edge_id,
                    "source_chart": source,
                    "target_chart": target,
                    "changed_factor": changed,
                    "support_lifts_on_overlap": len(support_indices),
                    "symbol_rows": len(rows),
                    "symbol_rank_mod_31": len(pivots),
                    "constant_symbol_syzygy_dimension": len(basis),
                }
            )
            edge_basis_rows.extend(emit_basis_rows(surface, "edge", location, basis))
            for endpoint, chart in (("source", source), ("target", target)):
                _, chart_basis = chart_data[chart]
                for source_basis_index, vector in enumerate(chart_basis):
                    coefficients = coordinates_in_nullspace(vector, free, basis)
                    for target_basis_index, value in enumerate(coefficients):
                        if value:
                            restriction_rows.append(
                                {
                                    "surface": f"Z_{surface}",
                                    "edge_id": edge_id,
                                    "endpoint": endpoint,
                                    "chart": chart,
                                    "edge_basis_index": target_basis_index,
                                    "chart_basis_index": source_basis_index,
                                    "value_mod_31": value,
                                }
                            )

    output_directory.mkdir(parents=True, exist_ok=True)
    artifact_rows = [
        ("K_degree2_action_elements.csv", group_rows),
        ("K_stable_product_chart_refinement.csv", stable_rows),
        ("surface_conormal_symbol_summary.csv", surface_summary),
        ("chart_conormal_symbol_summary.csv", chart_summary),
        ("edge_conormal_symbol_summary.csv", edge_summary),
        ("chart_constant_symbol_syzygy_basis.coo.csv", chart_basis_rows),
        ("edge_constant_symbol_syzygy_basis.coo.csv", edge_basis_rows),
        ("chart_to_edge_constant_symbol_restrictions.coo.csv", restriction_rows),
    ]
    paths = {}
    for name, rows in artifact_rows:
        path = output_directory / name
        write_csv(path, rows, list(rows[0]))
        paths[name] = path

    full_contract = {
        "schema": "section12-full-chart-local-HS-syzygy-next-contract-v1",
        "completed_finite_quotient": (
            "the residue conormal-symbol kernels for all six surfaces on all 16 charts "
            "and their 32 adjacent overlaps, including exact restriction matrices"
        ),
        "full_sheaf_map": "R_(s,c)^162 -> J_(s,c)",
        "affine_surface_ring": (
            "R_(s,c)=F_31[eight affine factor coordinates]/(four plane-cubic equations,"
            " dehomogenized s_s, dehomogenized s_(s+1))"
        ),
        "generator_ideal": (
            "J_(s,c) is generated by the 81 normalized K pullbacks of each emitted "
            "X_s and Q_s; its coefficients are exactly reconstructible from the emitted "
            "degree-two action table and the upstream coefficient vectors"
        ),
        "smallest_remaining_exact_operations": [
            "saturate J_(s,c) and compare it with the 81-lift curvilinear ideal I_W on every chart",
            "compute a module Groebner basis for ker(R_(s,c)^162 -> J_(s,c)) in a fixed term order",
            "localize those module bases on each chart overlap and solve the change-of-basis matrices",
            "repeat the restriction calculation on the 1296-chart K-stable refinement before imposing K-equivariance",
        ],
        "why_the_32_edges_are_not_the_K_descent_system": (
            "the exact three-torsion matrices send Z and Y to general linear forms; "
            "the K-orbit of the 16 product charts has 1296 distinct members"
        ),
        "no_guessed_block_size": True,
    }
    contract_path = output_directory / "full_sheaf_syzygy_next_contract.json"
    contract_path.write_text(json.dumps(full_contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    paths[contract_path.name] = contract_path

    certificate = {
        "schema": "section12-HS-chart-local-syzygy-audit-v1",
        "status": (
            "degree2-K-action-reconstructed;all-six-supportwise-conormal-symbols-surjective;"
            "16-chart-and-32-edge-constant-symbol-syzygies-and-restrictions-emitted;"
            "16-chart-atlas-not-K-stable;full-sheaf-syzygies-and-global-Cech-open"
        ),
        "field": "F_31",
        "K_action": {
            "group_order": 81,
            "normalized_group_law_checks": group_law_checks,
            "source": "induced quadratic action of the exact cached 3x3 plane-cubic translation matrices",
            "translation_coefficients_fabricated": False,
        },
        "finite_symbol_result": {
            "surfaces": 6,
            "generators_per_surface": 162,
            "support_lifts_per_surface": 81,
            "conormal_symbol_dimension_per_support": 2,
            "joint_constant_symbol_ranks": [int(row["conormal_symbol_rank_mod_31"]) for row in surface_summary],
            "every_support_block_rank": 2,
            "automorphy_block_checks": sum(
                int(row["automorphy_block_checks"]) for row in surface_summary
            ),
            "direct_quadratic_transport_checks": sum(
                int(row["direct_quadratic_transport_checks"]) for row in surface_summary
            ),
            "chart_system_count": len(chart_summary),
            "edge_system_count": len(edge_summary),
            "restriction_matrices_emitted": True,
        },
        "atlas_result": {
            "original_product_charts": 16,
            "original_adjacent_edges": 32,
            "K_stable": False,
            "distinct_K_translates_of_product_charts": len(stable_rows),
            "consequence": "the original 32 edges cannot carry the global K-linearization by themselves",
        },
        "claim_boundary": {
            "proved": [
                "the exact normalized degree-two K action reconstructs all 162 generator germs",
                "the 162 generators span both residue conormal directions at each of the 81 supports",
                "by Nakayama, the translate-generated ideal J and I_W have the same stalk at every support lift",
                "the joint constant-coefficient conormal-symbol ranks are computed rather than assumed",
                "the constant-coefficient conormal-symbol syzygy bases and chart-to-edge restrictions are exact",
                "the 16 product charts are not K-stable and have 1296 distinct K translates",
            ],
            "not_claimed": [
                "that the translate-generated ideal J equals the full ideal I_W after saturation",
                "a polynomial-module syzygy basis over any affine surface chart ring",
                "K-equivariant descent on only the original 16 charts",
                "the global HS Cech cocycle, Yoneda numerator, transgression T, or a Hodge result",
            ],
        },
        "next_contract": full_contract,
        "provenance": {
            "global_HS_certificate": str(GLOBAL_CERTIFICATE.relative_to(ROOT)),
            "global_HS_certificate_sha256": sha256(GLOBAL_CERTIFICATE),
            "descent_certificate": str(DESCENT.relative_to(ROOT)),
            "descent_certificate_sha256": sha256(DESCENT),
        },
    }
    certificate_path = output_directory / "hs_chart_local_syzygy_audit_certificate.json"
    certificate_path.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    paths[certificate_path.name] = certificate_path

    manifest = {
        "schema": "section12-HS-chart-local-syzygy-audit-manifest-v1",
        "generator": str(Path(__file__).relative_to(ROOT)),
        "generator_sha256": sha256(Path(__file__)),
        "files": {
            name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for name, path in sorted(paths.items())
        },
    }
    manifest_path = output_directory / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return certificate


if __name__ == "__main__":
    print(json.dumps(write_artifacts(), indent=2, sort_keys=True))
