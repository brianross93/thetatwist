"""Emit six exact local Hartshorne--Serre seeds for the main Markman lane.

The independent F31 incidence fibre is not the earlier cyclic six-divisor
fibre.  Its six theta divisors are six independent translates in
``E[3]^4/K``.  Consequently, one cannot transport the recorded seed on
``Z_2=D_2 intersect D_3`` to the other five surfaces.  This module checks
that statement from the exact translation classes before constructing any
new data.

The cached theta coefficients nevertheless determine enough geometry to
construct the missing *local* inputs directly.  For every

    Z_s = D_s intersect D_(s+1),

the producer selects an F31-rational smooth support point away from the two
adjacent conductor divisors, chooses fixed-frame surface coordinates
``(x_s,y_s)``, and defines the curvilinear double point

    W_s = (x_s, y_s^2).

It then evaluates

    delta_s = s_(s-1) s_(s+2)

on ``W_s`` and emits the actual completed-local Cech representative

    [[1, delta_s/(x_s*y_s^2)], [0, 1]].

The original Z2 support, tangent, coordinates, and residue orientation are
preserved exactly.  No old node-wedge matrix is identified with the global
Yoneda transgression ``T``.  The emitted cocycles are completed-local
cocycles, not a fabricated global theta-chart resolution.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
from array import array
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
P = 31
POINT_COUNT = 36
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter"
)

COEFFICIENTS = ROOT / (
    "results/section12_instantiation/independent_incidence_fiber/"
    "candidate_section_coefficients.csv"
)
INCIDENCE_CERTIFICATE = ROOT / (
    "results/section12_instantiation/independent_incidence_fiber/"
    "independent_incidence_certificate.json"
)
DESCENT_CERTIFICATE = ROOT / (
    "results/section12_instantiation/descent_section/"
    "descent_section_certificate.json"
)
MAIN_PAIR_CERTIFICATE = ROOT / (
    "results/section12_instantiation/cm_bridge_rank2_chow_search/"
    "rank2_chow_solubility_certificate.json"
)
SEED_CERTIFICATE = ROOT / (
    "results/section12_instantiation/cm_bridge_surface_hartshorne_serre_seed/"
    "f31_surface_double_point_seed.json"
)
SEED_JETS = ROOT / (
    "results/section12_instantiation/f31_rkeep_point_sector/"
    "rkeep_local_section_jets.csv"
)


class SixSeedContractError(RuntimeError):
    """Raised when a bound exact input changes."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def add_points(left, right):
    """Add points on y^2=x^3+1 over F31; None is the origin."""

    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and (y1 + y2) % P == 0:
        return None
    if left == right:
        if y1 == 0:
            return None
        slope = 3 * x1 * x1 * pow(2 * y1, -1, P) % P
    else:
        slope = (y2 - y1) * pow((x2 - x1) % P, -1, P) % P
    x3 = (slope * slope - x1 - x2) % P
    return x3, (slope * (x1 - x3) - y1) % P


def multiply_point(scalar, point):
    result = None
    for _ in range(scalar):
        result = add_points(result, point)
    return result


def rank_mod(matrix, prime=P):
    work = [[int(value) % prime for value in row] for row in matrix]
    rank = 0
    columns = len(work[0]) if work else 0
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [value * inverse % prime for value in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                (left - multiple * right) % prime
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
    return rank


def determinant4(matrix):
    result = 0
    for permutation in itertools.permutations(range(4)):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(4)
            for right in range(left + 1, 4)
        )
        term = -1 if inversions % 2 else 1
        for row in range(4):
            term *= matrix[row][permutation[row]]
        result += term
    return result % P


def solve4(matrix, target):
    augmented = [
        [int(value) % P for value in row] + [int(target[index]) % P]
        for index, row in enumerate(matrix)
    ]
    for column in range(4):
        pivot = next(
            row for row in range(column, 4) if augmented[row][column]
        )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = pow(augmented[column][column], -1, P)
        augmented[column] = [value * inverse % P for value in augmented[column]]
        for row in range(4):
            if row == column or not augmented[row][column]:
                continue
            multiple = augmented[row][column]
            augmented[row] = [
                (left - multiple * right) % P
                for left, right in zip(augmented[row], augmented[column])
            ]
    return [augmented[index][4] for index in range(4)]


def dot(left, right):
    return sum(a * b for a, b in zip(left, right)) % P


def dual_mul(left, right):
    return (
        left[0] * right[0] % P,
        (left[0] * right[1] + left[1] * right[0]) % P,
    )


def dual_inverse(value):
    inverse = pow(value[0], -1, P)
    return inverse, -value[1] * inverse * inverse % P


def load_inputs():
    incidence = load_json(INCIDENCE_CERTIFICATE)
    descent = load_json(DESCENT_CERTIFICATE)
    main_pair = load_json(MAIN_PAIR_CERTIFICATE)
    seed = load_json(SEED_CERTIFICATE)

    if incidence.get("status") != "independent-incidence-fiber-certified":
        raise SixSeedContractError("the independent incidence fibre changed")
    if descent.get("status") != "certified":
        raise SixSeedContractError("the descended theta section changed")
    if not main_pair.get("status", "").startswith(
        "determinant-matched-HS-pair-certified"
    ):
        raise SixSeedContractError("the quotient-valid determinant-2H HS pair changed")
    hs = main_pair["determinant_matched_hartshorne_serre"]
    if hs.get("split_bundle") != "V_split=O_S(3H) plus O_S(-H)":
        raise SixSeedContractError("the main-lane split bundle changed")
    if hs.get("exact_difference") != "c2(V_HS)-c2(V_split)=[W]=2[p] in CH_0(S)":
        raise SixSeedContractError("the main-lane point correction changed")
    if seed.get("surface", {}).get("name") != "Z_2":
        raise SixSeedContractError("the preserved Z2 seed changed")

    coefficients = {label: [] for label in range(6)}
    with COEFFICIENTS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            coefficients[int(row["label"])].append(
                int(row["coefficient_mod_31"]) % P
            )
    if any(len(coefficients[label]) != 81 for label in range(6)):
        raise SixSeedContractError("a translated theta coefficient row is incomplete")

    seed_jets = {}
    with SEED_JETS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["quotient_orbit"]) == 0 and int(row["orbit_lift_index"]) == 0:
                seed_jets[int(row["label"])] = {
                    "value": int(row["value_mod_31"]) % P,
                    "gradient": [
                        int(row[f"direction_{direction}_derivative_mod_31"]) % P
                        for direction in range(4)
                    ],
                }
    if set(seed_jets) != set(range(6)):
        raise SixSeedContractError("the preserved seed lost a section jet")
    return incidence, descent, main_pair, seed, coefficients, seed_jets


def curve_data():
    points = [None] + [
        (x, y)
        for x in range(P)
        for y in range(P)
        if (y * y - x**3 - 1) % P == 0
    ]
    if len(points) != POINT_COUNT:
        raise ArithmeticError("the F31 elliptic point count changed")
    values = [(0, 0, 1) if point is None else (1, point[0], point[1]) for point in points]
    # The fixed invariant vector field used by f31_principal_parts.py.
    derivatives = [
        (0, P - 1, 0)
        if point is None
        else (0, 2 * point[1] % P, 3 * point[0] * point[0] % P)
        for point in points
    ]
    return points, values, derivatives


def section_value_table(coefficients, point_values):
    """Evaluate one four-factor tensor on all 36^4 rational product points."""

    first = [
        [
            sum(
                value[a] * coefficients[((a * 3 + b) * 3 + c) * 3 + d]
                for a in range(3)
            )
            % P
            for b in range(3)
            for c in range(3)
            for d in range(3)
        ]
        for value in point_values
    ]
    second = []
    for first_row in first:
        for value in point_values:
            second.append(
                [
                    sum(
                        value[b] * first_row[(b * 3 + c) * 3 + d]
                        for b in range(3)
                    )
                    % P
                    for c in range(3)
                    for d in range(3)
                ]
            )
    third = []
    for second_row in second:
        for value in point_values:
            third.append(
                tuple(
                    sum(value[c] * second_row[c * 3 + d] for c in range(3)) % P
                    for d in range(3)
                )
            )
    result = array("B")
    append = result.append
    for third_row in third:
        for value in point_values:
            append(sum(value[d] * third_row[d] for d in range(3)) % P)
    if len(result) != POINT_COUNT**4:
        raise ArithmeticError("a section value table has the wrong size")
    return result


def decode_index(index):
    result = []
    for power in (POINT_COUNT**3, POINT_COUNT**2, POINT_COUNT, 1):
        result.append(index // power)
        index %= power
    return tuple(result)


def flatten_point(point):
    return ((point[0] * POINT_COUNT + point[1]) * POINT_COUNT + point[2]) * POINT_COUNT + point[3]


def section_gradient(coefficients, point, point_values, point_derivatives):
    result = []
    for direction in range(4):
        total = 0
        for flat, coefficient in enumerate(coefficients):
            if not coefficient:
                continue
            digits = tuple((flat // 3 ** (3 - factor)) % 3 for factor in range(4))
            term = coefficient
            for factor, digit in enumerate(digits):
                vector = (
                    point_derivatives[point[factor]]
                    if factor == direction
                    else point_values[point[factor]]
                )
                term *= vector[digit]
            total += term
        result.append(total % P)
    return result


def kernel_orbit_data(descent, points):
    point_index = {point: index for index, point in enumerate(points)}
    torsion_p = (0, 1)
    torsion_q = (3, 11)
    e3 = {
        (a, b): add_points(multiply_point(a, torsion_p), multiply_point(b, torsion_q))
        for a in range(3)
        for b in range(3)
    }
    shifts = {
        key: [point_index[add_points(point, translation)] for point in points]
        for key, translation in e3.items()
    }
    kernel_rows = descent["kernel_rows"]
    kernel = []
    for scalars in itertools.product(range(3), repeat=4):
        row = [
            sum(scalars[index] * kernel_rows[index][column] for index in range(4)) % 3
            for column in range(8)
        ]
        kernel.append(tuple((row[2 * factor], row[2 * factor + 1]) for factor in range(4)))
    if len(set(kernel)) != 81:
        raise ArithmeticError("the quotient kernel no longer has 81 elements")

    def orbit(point):
        return sorted(
            {
                tuple(shifts[key[factor]][point[factor]] for factor in range(4))
                for key in kernel
            }
        )

    return orbit


def quotient_vector(translation):
    return tuple(value for pair in translation for value in pair)


def subtract_vectors(left, right):
    return tuple((a - b) % 3 for a, b in zip(left, right))


def adjacent_transport_audit(incidence, descent):
    translations = [
        quotient_vector(translation)
        for translation in incidence["translations_e3_keys"]
    ]
    differences = [
        subtract_vectors(translations[(index + 1) % 6], translations[index])
        for index in range(6)
    ]
    kernel_rows = descent["kernel_rows"]
    base_rank = rank_mod(kernel_rows, prime=3)

    def same_class(left, right):
        difference = subtract_vectors(left, right)
        return rank_mod([*kernel_rows, difference], prime=3) == base_rank

    rows = []
    for source in range(6):
        equivalent = [
            target for target in range(6) if same_class(differences[source], differences[target])
        ]
        rows.append(
            {
                "surface": f"Z_{source}=D_{source} intersection D_{(source + 1) % 6}",
                "adjacent_difference_E3_quotient": json.dumps(list(differences[source])),
                "translation_equivalent_surface_indices": json.dumps(equivalent),
                "transport_from_Z2_available": same_class(differences[source], differences[2]),
            }
        )
    if any(json.loads(row["translation_equivalent_surface_indices"]) != [index] for index, row in enumerate(rows)):
        raise ArithmeticError("two independent surface pairs unexpectedly became translates")
    return rows


def coordinate_choice(normals):
    for first, second in itertools.combinations(range(4), 2):
        x_row = [int(index == first) for index in range(4)]
        y_row = [int(index == second) for index in range(4)]
        matrix = [*normals, x_row, y_row]
        if determinant4(matrix):
            tangent = solve4(matrix, [0, 0, 0, 1])
            return x_row, y_row, tangent, determinant4(matrix)
    raise ArithmeticError("a smooth surface point has no fixed-frame coordinate pair")


def select_six_seeds(seed, seed_jets, coefficients, value_tables, point_values, point_derivatives, orbit):
    seed_point = tuple(seed["point"]["lift_product_point_indices"])
    seed_orbit = orbit(seed_point)
    if len(seed_orbit) != 81:
        raise ArithmeticError("the preserved seed does not have a free K orbit")
    used = {seed_orbit[0]}
    chosen = {}

    for surface in range(6):
        if surface == 2:
            point = seed_point
            gradients = {
                label: seed_jets[label]["gradient"] for label in range(6)
            }
            tangent = list(seed["curvilinear_double_point"]["tangent_vector"])
            x_row = [
                (15 * gradients[0][index] - 20 * gradients[5][index]) % P
                for index in range(4)
            ]
            y_row = [14 * gradients[0][index] % P for index in range(4)]
            determinant = determinant4([gradients[2], gradients[3], x_row, y_row])
            point_values_row = [seed_jets[label]["value"] for label in range(6)]
            if determinant != 23 or dot(x_row, tangent) or dot(y_row, tangent) != 1:
                raise ArithmeticError("the preserved Z2 adapted coordinates changed")
            chosen[surface] = {
                "point": point,
                "quotient_orbit_representative": seed_orbit[0],
                "quotient_orbit": seed_orbit,
                "gradients": gradients,
                "values": point_values_row,
                "x_row": x_row,
                "y_row": y_row,
                "tangent": tangent,
                "coordinate_determinant": determinant,
                "selection": "preserved-certified-Z2-seed",
            }
            continue

        surface_next = (surface + 1) % 6
        boundary_left = (surface - 1) % 6
        boundary_right = (surface + 2) % 6
        selected = None
        for flat in range(POINT_COUNT**4):
            if value_tables[surface][flat] or value_tables[surface_next][flat]:
                continue
            if not value_tables[boundary_left][flat] or not value_tables[boundary_right][flat]:
                continue
            point = decode_index(flat)
            point_orbit = orbit(point)
            representative = point_orbit[0]
            if representative in used:
                continue
            gradients = {
                label: section_gradient(
                    coefficients[label], point, point_values, point_derivatives
                )
                for label in range(6)
            }
            normals = [gradients[surface], gradients[surface_next]]
            if rank_mod(normals) != 2:
                continue
            x_row, y_row, tangent, determinant = coordinate_choice(normals)
            selected = {
                "point": point,
                "quotient_orbit_representative": representative,
                "quotient_orbit": point_orbit,
                "gradients": gradients,
                "values": [value_tables[label][flat] for label in range(6)],
                "x_row": x_row,
                "y_row": y_row,
                "tangent": tangent,
                "coordinate_determinant": determinant,
                "selection": "first-smooth-boundary-avoiding-distinct-quotient-orbit",
            }
            break
        if selected is None:
            raise ArithmeticError(f"no rational local HS support was found on Z_{surface}")
        chosen[surface] = selected
        used.add(selected["quotient_orbit_representative"])

    if len(chosen) != 6 or len({chosen[index]["quotient_orbit_representative"] for index in chosen}) != 6:
        raise ArithmeticError("the six selected quotient supports are not distinct")
    return chosen


def point_projective(points, point):
    result = []
    for index in point:
        value = points[index]
        result.append([0, 1, 0] if value is None else [value[0], value[1], 1])
    return result


def write_artifacts(output_directory: Path = OUTPUT_DIRECTORY):
    incidence, descent, main_pair, seed, coefficients, seed_jets = load_inputs()
    points, point_values, point_derivatives = curve_data()
    value_tables = [
        section_value_table(coefficients[label], point_values) for label in range(6)
    ]

    # Bind the independently generated fixed-frame values and derivatives at
    # the preserved seed before extending the construction to five new points.
    seed_point = tuple(seed["point"]["lift_product_point_indices"])
    seed_flat = flatten_point(seed_point)
    for label in range(6):
        if value_tables[label][seed_flat] != seed_jets[label]["value"]:
            raise SixSeedContractError("the seed value disagrees with the theta tensor")
        if section_gradient(
            coefficients[label], seed_point, point_values, point_derivatives
        ) != seed_jets[label]["gradient"]:
            raise SixSeedContractError("the seed derivative disagrees with the fixed-frame jet")

    orbit = kernel_orbit_data(descent, points)
    transport_rows = adjacent_transport_audit(incidence, descent)
    chosen = select_six_seeds(
        seed,
        seed_jets,
        coefficients,
        value_tables,
        point_values,
        point_derivatives,
        orbit,
    )

    seed_rows = []
    restriction_rows = []
    support_orbit_rows = []
    cocycles = []
    for surface in range(6):
        data = chosen[surface]
        surface_labels = (surface, (surface + 1) % 6)
        boundary_labels = ((surface - 1) % 6, (surface + 2) % 6)
        tangent_derivatives = [
            dot(data["gradients"][label], data["tangent"]) for label in range(6)
        ]
        left = (
            data["values"][boundary_labels[0]],
            tangent_derivatives[boundary_labels[0]],
        )
        right = (
            data["values"][boundary_labels[1]],
            tangent_derivatives[boundary_labels[1]],
        )
        delta = dual_mul(left, right)
        if not delta[0]:
            raise ArithmeticError("a selected W meets its conductor boundary")
        delta_inverse = dual_inverse(delta)
        if dual_mul(delta, delta_inverse) != (1, 0):
            raise ArithmeticError("a conductor unit failed inversion")
        residue_coefficients = (delta[1], delta[0])
        if (
            residue_coefficients[0] * delta_inverse[0]
            + residue_coefficients[1] * delta_inverse[1]
        ) % P:
            raise ArithmeticError("the synchronized boundary direction has nonzero residue")

        seed_rows.append(
            {
                "surface": f"Z_{surface}",
                "surface_section_labels": json.dumps(list(surface_labels)),
                "boundary_section_labels": json.dumps(list(boundary_labels)),
                "source_point_indices": json.dumps(list(data["point"])),
                "source_projective_points_X_Y_Z": json.dumps(
                    point_projective(points, data["point"])
                ),
                "quotient_orbit_representative": json.dumps(
                    list(data["quotient_orbit_representative"])
                ),
                "quotient_orbit_size": len(data["quotient_orbit"]),
                "selection": data["selection"],
                "x_fixed_frame_coefficients": json.dumps(data["x_row"]),
                "y_fixed_frame_coefficients": json.dumps(data["y_row"]),
                "tangent_fixed_frame_coefficients": json.dumps(data["tangent"]),
                "adapted_coordinate_determinant_mod_31": data[
                    "coordinate_determinant"
                ],
                "delta0_mod_31": delta[0],
                "delta1_mod_31": delta[1],
                "delta_inverse0_mod_31": delta_inverse[0],
                "delta_inverse1_mod_31": delta_inverse[1],
                "residue_f0_coefficient_mod_31": residue_coefficients[0],
                "residue_f1_coefficient_mod_31": residue_coefficients[1],
            }
        )
        for lift_index, lift_point in enumerate(data["quotient_orbit"]):
            support_orbit_rows.append(
                {
                    "surface": f"Z_{surface}",
                    "K_orbit_lift_index": lift_index,
                    "source_point_indices": json.dumps(list(lift_point)),
                    "source_projective_points_X_Y_Z": json.dumps(
                        point_projective(points, lift_point)
                    ),
                    "transported_tangent_fixed_frame_coefficients": json.dumps(
                        data["tangent"]
                    ),
                    "transport_reason": (
                        "the four factor directions use the invariant elliptic vector field; "
                        "translation preserves their coefficients"
                    ),
                }
            )
        for label in range(6):
            restriction_rows.append(
                {
                    "surface": f"Z_{surface}",
                    "section": f"s_{label}",
                    "value_f0_mod_31": data["values"][label],
                    "tangent_f1_mod_31": tangent_derivatives[label],
                    "is_surface_equation": label in surface_labels,
                    "is_boundary_equation": label in boundary_labels,
                }
            )

        cocycles.append(
            {
                "surface": f"Z_{surface}",
                "completed_surface_ring": "F_31[[x_s,y_s]]",
                "double_point_ideal": "(x_s,y_s^2)",
                "charts": ["D(x_s)", "D(y_s)"],
                "overlap": "D(x_s*y_s)",
                "delta": {
                    "global_formula": f"s_{boundary_labels[0]}*s_{boundary_labels[1]}",
                    "restriction_to_W": list(delta),
                },
                "transition_Dx_to_Dy": [
                    ["1", "delta_s/(x_s*y_s^2)"],
                    ["0", "1"],
                ],
                "transition_Dy_to_Dx": [
                    ["1", "-delta_s/(x_s*y_s^2)"],
                    ["0", "1"],
                ],
                "local_resolution": {
                    "kernel_column": ["-y_s^2", "x_s"],
                    "quotient_row": ["x_s", "y_s^2"],
                },
                "residue_oracle": (
                    f"phi_{surface}(f0+f1*epsilon)="
                    f"{residue_coefficients[0]}*f0+{residue_coefficients[1]}*f1 mod 31"
                ),
                "scope": "actual completed-local HS cocycle; not a global theta-chart cocycle",
            }
        )

    if seed_rows[2]["delta0_mod_31"] != 28 or seed_rows[2]["delta1_mod_31"] != 17:
        raise ArithmeticError("the preserved Z2 conductor orientation changed")

    output_directory.mkdir(parents=True, exist_ok=True)
    seed_path = output_directory / "six_hs_local_seeds.csv"
    write_csv(seed_path, seed_rows, list(seed_rows[0]))
    restrictions_path = output_directory / "six_hs_section_restrictions.csv"
    write_csv(restrictions_path, restriction_rows, list(restriction_rows[0]))
    support_orbits_path = output_directory / "six_hs_support_K_orbits.csv"
    write_csv(support_orbits_path, support_orbit_rows, list(support_orbit_rows[0]))
    cocycle_path = output_directory / "six_hs_completed_local_cocycles.json"
    cocycle_path.write_text(
        json.dumps(cocycles, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    transport_path = output_directory / "surface_pair_transport_audit.csv"
    write_csv(transport_path, transport_rows, list(transport_rows[0]))

    next_contract = {
        "schema": "section12-global-hs-chart-producer-input-v1",
        "smallest_next_producer": "global_HS_chart_resolution",
        "available_exact_inputs": [
            "six rational quotient supports through their full K-orbit representatives",
            "six fixed-frame tangent double points W_s=(x_s,y_s^2)",
            "six boundary units delta_s and exact local residue functionals",
            "six completed-local Hartshorne--Serre Cech matrices",
            "all six translated theta coefficient tensors",
        ],
        "next_operations": [
            "form the 81-lift K-invariant ideal of each W_s on the product cover",
            "solve the connecting-homomorphism lift of delta_s/(x_s*y_s^2) in a common theta-chart cover",
            "emit the global rank-two HS transition matrices with their K linearizations",
            "embed the six presentations into a common ambient first-syzygy complex",
        ],
        "still_not_determined": [
            "the global off-diagonal Yoneda numerator nu_s,r,c|W_s",
            "the six surface-residue columns",
            "the transgression matrix T",
            "equality in CH_0 with the original twelve retained-node cycle",
        ],
        "important_scope": (
            "the six new W_s give the exact degree-12 cohomological point correction; "
            "their sum has not been identified with the original retained-node cycle in CH_0"
        ),
    }
    contract_path = output_directory / "next_global_hs_chart_input.json"
    contract_path.write_text(
        json.dumps(next_contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    certificate = {
        "schema": "section12-six-hs-local-seed-emitter-v1",
        "status": (
            "cyclic-seed-transport-closed-for-independent-fibre;"
            "six-rational-local-HS-seeds-emitted;Z2-seed-preserved;"
            "six-completed-local-cocycles-emitted;global-Yoneda-numerators-open"
        ),
        "field": "F_31",
        "main_lane_pair": main_pair["determinant_matched_hartshorne_serre"],
        "surface_family": "Z_s=D_s intersection D_(s+1), indices modulo 6",
        "transport_verdict": {
            "seed_transport_to_other_surfaces": False,
            "reason": (
                "the six adjacent divisor-difference classes in E[3]^4/K are pairwise "
                "distinct; a translation carrying an ordered theta pair to another must "
                "preserve this difference class"
            ),
            "replacement": "construct each local seed directly from its exact theta coefficients",
        },
        "local_seed_count": len(seed_rows),
        "distinct_quotient_support_count": len(
            {row["quotient_orbit_representative"] for row in seed_rows}
        ),
        "total_double_point_length": 2 * len(seed_rows),
        "all_boundary_constants_units": all(row["delta0_mod_31"] for row in seed_rows),
        "preserved_Z2_orientation": {
            "point_indices": seed_rows[2]["source_point_indices"],
            "tangent": seed_rows[2]["tangent_fixed_frame_coefficients"],
            "delta": [seed_rows[2]["delta0_mod_31"], seed_rows[2]["delta1_mod_31"]],
            "matches_existing_oracle": True,
        },
        "cohomological_ledger": (
            "sum_s length(W_s)=12, hence the top-degree cohomology correction is 12[pt]"
        ),
        "CH0_ledger": (
            "sum_s[W_s]=2*sum_s[p_s] exactly for the emitted supports; no equality with "
            "the original twelve retained points is claimed"
        ),
        "global_output": {
            "global_HS_theta_chart_cocycles_computed": False,
            "global_Yoneda_numerators_computed": False,
            "T_computed": False,
            "reason": (
                "completed-local Cech data determine the residue oracles but not the "
                "global connecting-homomorphism representatives or off-diagonal chain maps"
            ),
        },
        "artifacts": {
            "seeds": seed_path.name,
            "section_restrictions": restrictions_path.name,
            "support_K_orbits": support_orbits_path.name,
            "completed_local_cocycles": cocycle_path.name,
            "transport_audit": transport_path.name,
            "next_global_input": contract_path.name,
        },
        "scope": {
            "proved": [
                "the actual independent fibre does not admit cyclic transport of the Z2 ordered surface pair",
                "six distinct F31-rational smooth surface supports exist",
                "each support defines a boundary-avoiding curvilinear length-two W_s",
                "each delta_s is a unit on W_s and fixes an exact local residue oracle",
                "the six completed-local Hartshorne--Serre transition cocycles are explicit",
            ],
            "not_claimed": [
                "a global theta-chart HS transition cocycle",
                "an off-diagonal global Yoneda class",
                "a value or rank for T",
                "vanishing of the six corrected witnesses",
                "an answer to Question 11.4 or the Hodge conjecture",
            ],
        },
    }
    certificate_path = output_directory / "six_hs_local_seed_certificate.json"
    certificate_path.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    files = [
        seed_path,
        restrictions_path,
        support_orbits_path,
        cocycle_path,
        transport_path,
        contract_path,
        certificate_path,
    ]
    manifest = {
        "schema": "section12-six-hs-local-seed-emitter-manifest-v1",
        "generator": str(Path(__file__).relative_to(ROOT)),
        "generator_sha256": sha256(Path(__file__)),
        "sources": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                COEFFICIENTS,
                INCIDENCE_CERTIFICATE,
                DESCENT_CERTIFICATE,
                MAIN_PAIR_CERTIFICATE,
                SEED_CERTIFICATE,
                SEED_JETS,
            )
        },
        "files": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in files
        },
    }
    manifest_path = output_directory / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return certificate


if __name__ == "__main__":
    print(json.dumps(write_artifacts(), indent=2, sort_keys=True))
