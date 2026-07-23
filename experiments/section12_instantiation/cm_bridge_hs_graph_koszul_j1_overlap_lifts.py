"""Emit exact graph--Koszul ``J1`` lifts on all original chart edges.

The edge carrier uses the source-chart affine coordinates.  On the changed
factor the target coordinates are

    u' = u / v,    v' = 1 / v.

For a target quotient restriction ``R`` and target multiplication matrices
``M'_u,M'_v``, the changed part of the degree-one chain map is stored as the
factored Laurent block

                 target u'                  target v'
    edge u       v^-1 R                     0
    edge v      -v^-1 R M'_u              -v^-1 R M'_v.

Every unchanged generator maps diagonally by ``R``.  A source endpoint is
the eight-block diagonal map with coefficient ``R``.  The producer expands
these programs in memory and verifies, for all 384 endpoint maps and all
eight graph generators,

    d_edge J1 = J0 d_chart

as an identity of Laurent-polynomial matrices over ``F_31``.  It writes the
coefficient blocks and the factored programs, not the enormous expanded
``8*n_edge`` by ``8*n_chart`` matrices.

This is only the graph--Koszul degree-one lift.  The full Shamash ``J2`` and
``J3`` maps additionally need an exact change-of-equations matrix for the
six chart regular sequences and the resulting comparison homotopies.  Those
maps are not inferred from quotient restrictions and are not zero-filled.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import os
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
P = 31
SURFACE_COUNT = 6
FACTOR_COUNT = 4
VARIABLE_COUNT = 8
EDGE_COUNT = 32

INPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_iw_hermite_border_basis"
)
OVERLAP_DIRECTORY = INPUT_DIRECTORY / "overlaps"
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_graph_koszul_j1_overlap_lifts"
)


class HSGraphKoszulJ1Error(RuntimeError):
    """An input, rational chart transition, or chain identity failed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def array_sha256(value: np.ndarray) -> str:
    return sha256_bytes(
        np.ascontiguousarray(np.asarray(value, dtype=np.uint8)).tobytes()
    )


def canonical_sha256(value) -> str:
    return sha256_bytes(
        json.dumps(
            value, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    )


def write_csv(path: Path, rows, fieldnames=None) -> None:
    rows = list(rows)
    if not rows and fieldnames is None:
        raise ValueError("fieldnames are required for an empty CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames or list(rows[0].keys())
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json_atomic(path: Path, value) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def chart_name(bits) -> str:
    return "C" + "".join(map(str, bits))


def coordinate_name(variable: int) -> str:
    return f"{'u' if variable % 2 == 0 else 'v'}{variable // 2}"


def load_inputs():
    charts = {}
    overlaps = {}
    input_paths = []

    for surface in range(SURFACE_COUNT):
        for bits in itertools.product((0, 1), repeat=FACTOR_COUNT):
            chart = chart_name(bits)
            binary = INPUT_DIRECTORY / f"Z_{surface}_{chart}.npz"
            summary = INPUT_DIRECTORY / f"Z_{surface}_{chart}.json"
            if not binary.exists() or not summary.exists():
                raise HSGraphKoszulJ1Error(
                    f"missing chart input Z_{surface}/{chart}"
                )
            record = json.loads(summary.read_text(encoding="utf-8"))
            if record.get("schema") != "section12-IW-Hermite-border-basis-v1":
                raise HSGraphKoszulJ1Error("a chart schema changed")
            with np.load(binary) as artifact:
                field = int(artifact["field"][0])
                multiplication = np.asarray(
                    artifact["multiplication_matrices"], dtype=np.int64
                ) % P
                dimension = int(artifact["standard_monomials"].shape[0])
            if field != P or multiplication.shape != (
                VARIABLE_COUNT,
                dimension,
                dimension,
            ):
                raise HSGraphKoszulJ1Error(
                    "a chart multiplication artifact changed shape or field"
                )
            charts[(surface, chart)] = {
                "bits": tuple(map(int, bits)),
                "dimension": dimension,
                "multiplication": multiplication,
                "binary": binary,
                "summary": summary,
            }
            input_paths.extend((binary, summary))

    for surface in range(SURFACE_COUNT):
        for edge in range(EDGE_COUNT):
            binary = OVERLAP_DIRECTORY / f"Z_{surface}_E{edge:02d}.npz"
            summary = OVERLAP_DIRECTORY / f"Z_{surface}_E{edge:02d}.json"
            if not binary.exists() or not summary.exists():
                raise HSGraphKoszulJ1Error(
                    f"missing overlap input Z_{surface}/E{edge:02d}"
                )
            record = json.loads(summary.read_text(encoding="utf-8"))
            if (
                record.get("schema")
                != "section12-IW-Hermite-overlap-restriction-v1"
            ):
                raise HSGraphKoszulJ1Error("an overlap schema changed")
            if not all(record["exact_checks"].values()):
                raise HSGraphKoszulJ1Error(
                    "an upstream overlap exact check is false"
                )
            source = record["source_chart"]
            target = record["target_chart"]
            changed = int(record["changed_factor"])
            source_bits = charts[(surface, source)]["bits"]
            target_bits = charts[(surface, target)]["bits"]
            differences = [
                index
                for index, pair in enumerate(zip(source_bits, target_bits))
                if pair[0] != pair[1]
            ]
            if differences != [changed]:
                raise HSGraphKoszulJ1Error(
                    "an edge no longer changes exactly its recorded factor"
                )
            if source_bits[changed] != 0 or target_bits[changed] != 1:
                raise HSGraphKoszulJ1Error(
                    "the stored edge orientation is not source Z to target Y"
                )
            with np.load(binary) as artifact:
                field = int(artifact["field"][0])
                multiplication = np.asarray(
                    artifact["multiplication_matrices"], dtype=np.int64
                ) % P
                source_map = np.asarray(
                    artifact["source_restriction"], dtype=np.int64
                ) % P
                target_map = np.asarray(
                    artifact["target_restriction"], dtype=np.int64
                ) % P
            edge_dimension = source_map.shape[0]
            expected = (
                (VARIABLE_COUNT, edge_dimension, edge_dimension),
                (
                    edge_dimension,
                    charts[(surface, source)]["dimension"],
                ),
                (
                    edge_dimension,
                    charts[(surface, target)]["dimension"],
                ),
            )
            if field != P or (
                multiplication.shape,
                source_map.shape,
                target_map.shape,
            ) != expected:
                raise HSGraphKoszulJ1Error(
                    "an overlap multiplication or restriction shape changed"
                )
            overlaps[(surface, edge)] = {
                "source": source,
                "target": target,
                "changed": changed,
                "dimension": edge_dimension,
                "multiplication": multiplication,
                "source_map": source_map,
                "target_map": target_map,
                "binary": binary,
                "summary": summary,
            }
            input_paths.extend((binary, summary))

    if len(charts) != 96 or len(overlaps) != 192:
        raise HSGraphKoszulJ1Error("the chart or overlap census changed")
    return charts, overlaps, tuple(sorted(set(input_paths)))


def _add_polynomial_matrix(polynomial, exponent, matrix) -> None:
    exponent = tuple(map(int, exponent))
    value = np.asarray(matrix, dtype=np.int64) % P
    if exponent in polynomial:
        value = (polynomial[exponent] + value) % P
    if np.any(value):
        polynomial[exponent] = value
    else:
        polynomial.pop(exponent, None)


def _shift(exponent, variable, amount=1):
    result = list(map(int, exponent))
    result[int(variable)] += int(amount)
    return tuple(result)


def _pullback_coordinate(endpoint, variable, changed):
    """Return Laurent monomials for one chart coordinate on the edge."""

    zero = (0,) * VARIABLE_COUNT
    if endpoint == "source" or variable not in (2 * changed, 2 * changed + 1):
        return [(_shift(zero, variable), 1)]
    u = 2 * changed
    v = u + 1
    if variable == u:
        return [(_shift(_shift(zero, u), v, -1), 1)]
    return [(_shift(zero, v, -1), 1)]


def endpoint_program(
    surface,
    edge,
    endpoint,
    chart,
    changed,
    restriction,
    chart_multiplication,
    coefficient_blocks,
):
    """Return the exact factored J1 rows for one endpoint."""

    map_id = f"Z{surface}_E{edge:02d}_{endpoint}"
    base_key = f"{map_id}_J0"
    coefficient_blocks[base_key] = np.asarray(restriction, dtype=np.uint8)
    block_hashes = {base_key: array_sha256(restriction)}
    u = 2 * changed
    v = u + 1
    derived = {}
    if endpoint == "target":
        for variable in (u, v):
            key = f"{map_id}_J0M{variable}"
            matrix = restriction @ chart_multiplication[variable] % P
            coefficient_blocks[key] = np.asarray(matrix, dtype=np.uint8)
            block_hashes[key] = array_sha256(matrix)
            derived[variable] = key

    rows = []

    def add(domain, codomain, exponent, scalar, key, kind, multiplier=-1):
        rows.append(
            {
                "map_id": map_id,
                "surface": f"Z_{surface}",
                "edge_id": edge,
                "endpoint": endpoint,
                "chart": chart,
                "changed_factor": changed,
                "domain_generator": domain,
                "domain_coordinate": coordinate_name(domain),
                "codomain_generator": codomain,
                "codomain_coordinate": coordinate_name(codomain),
                "laurent_exponents_source_coordinates": json.dumps(
                    list(map(int, exponent)), separators=(",", ":")
                ),
                "scalar_mod_31": int(scalar) % P,
                "coefficient_block_key": key,
                "coefficient_block_kind": kind,
                "chart_multiplication_index": multiplier,
                "coefficient_block_shape": (
                    f"{restriction.shape[0]}x{restriction.shape[1]}"
                ),
                "coefficient_block_sha256": block_hashes[key],
            }
        )

    zero = (0,) * VARIABLE_COUNT
    if endpoint == "source":
        for variable in range(VARIABLE_COUNT):
            add(variable, variable, zero, 1, base_key, "J0")
        return rows

    denominator = _shift(zero, v, -1)
    for variable in range(VARIABLE_COUNT):
        if variable not in (u, v):
            add(variable, variable, zero, 1, base_key, "J0")
    add(u, u, denominator, 1, base_key, "J0")
    add(u, v, denominator, -1, derived[u], "J0_TIMES_M_CHART", u)
    add(v, v, denominator, -1, derived[v], "J0_TIMES_M_CHART", v)
    return rows


def verify_endpoint_chain_identity(
    rows,
    coefficient_blocks,
    endpoint,
    changed,
    restriction,
    chart_multiplication,
    edge_multiplication,
):
    """Verify ``d_edge J1 = J0 d_chart`` in the Laurent polynomial ring."""

    by_domain = {variable: [] for variable in range(VARIABLE_COUNT)}
    for row in rows:
        by_domain[int(row["domain_generator"])].append(row)

    generator_hashes = []
    zero = (0,) * VARIABLE_COUNT
    for domain in range(VARIABLE_COUNT):
        left = {}
        for row in by_domain[domain]:
            codomain = int(row["codomain_generator"])
            exponent = tuple(
                map(
                    int,
                    json.loads(row["laurent_exponents_source_coordinates"]),
                )
            )
            coefficient = (
                int(row["scalar_mod_31"])
                * np.asarray(
                    coefficient_blocks[row["coefficient_block_key"]],
                    dtype=np.int64,
                )
            ) % P
            _add_polynomial_matrix(
                left, _shift(exponent, codomain), coefficient
            )
            _add_polynomial_matrix(
                left,
                exponent,
                -edge_multiplication[codomain] @ coefficient,
            )

        right = {}
        for exponent, scalar in _pullback_coordinate(
            endpoint, domain, changed
        ):
            _add_polynomial_matrix(right, exponent, scalar * restriction)
        chart_product = restriction @ chart_multiplication[domain] % P
        _add_polynomial_matrix(right, zero, -chart_product)

        if set(left) != set(right):
            raise HSGraphKoszulJ1Error(
                "a J1 graph identity has different Laurent support"
            )
        for exponent in left:
            if not np.array_equal(left[exponent] % P, right[exponent] % P):
                raise HSGraphKoszulJ1Error(
                    "d_edge J1 != J0 d_chart for a graph generator"
                )
        generator_hashes.append(
            canonical_sha256(
                {
                    "domain": domain,
                    "laurent_support": [list(value) for value in sorted(left)],
                    "coefficient_hashes": [
                        array_sha256(left[value]) for value in sorted(left)
                    ],
                }
            )
        )
    return generator_hashes


def build():
    charts, overlaps, input_paths = load_inputs()
    program_rows = []
    verification_rows = []
    coefficient_blocks = {}
    generator_hashes = []

    for surface in range(SURFACE_COUNT):
        for edge in range(EDGE_COUNT):
            overlap = overlaps[(surface, edge)]
            for endpoint in ("source", "target"):
                chart = overlap[endpoint]
                chart_data = charts[(surface, chart)]
                restriction = overlap[f"{endpoint}_map"]
                rows = endpoint_program(
                    surface,
                    edge,
                    endpoint,
                    chart,
                    overlap["changed"],
                    restriction,
                    chart_data["multiplication"],
                    coefficient_blocks,
                )
                hashes = verify_endpoint_chain_identity(
                    rows,
                    coefficient_blocks,
                    endpoint,
                    overlap["changed"],
                    restriction,
                    chart_data["multiplication"],
                    overlap["multiplication"],
                )
                program_rows.extend(rows)
                generator_hashes.extend(hashes)
                verification_rows.append(
                    {
                        "map_id": rows[0]["map_id"],
                        "surface": f"Z_{surface}",
                        "edge_id": edge,
                        "endpoint": endpoint,
                        "chart": chart,
                        "changed_factor": overlap["changed"],
                        "J0_shape": (
                            f"{restriction.shape[0]}x{restriction.shape[1]}"
                        ),
                        "J1_factored_shape": (
                            f"{VARIABLE_COUNT * restriction.shape[0]}x"
                            f"{VARIABLE_COUNT * restriction.shape[1]}"
                        ),
                        "J1_term_count": len(rows),
                        "graph_generators_checked": VARIABLE_COUNT,
                        "residual_laurent_terms": 0,
                        "identity": "d_edge*J1=J0*d_chart",
                        "status": "EXACT_CHAIN_IDENTITY_VERIFIED",
                        "program_sha256": canonical_sha256(rows),
                    }
                )

    if len(verification_rows) != 384 or len(program_rows) != 3264:
        raise HSGraphKoszulJ1Error("the J1 endpoint or factor-term census changed")
    if len(coefficient_blocks) != 768 or len(generator_hashes) != 3072:
        raise HSGraphKoszulJ1Error(
            "the J1 coefficient-block or generator-check census changed"
        )

    input_rows = [
        {
            "path": str(path.relative_to(ROOT)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in input_paths
    ]
    return {
        "program_rows": program_rows,
        "verification_rows": verification_rows,
        "coefficient_blocks": coefficient_blocks,
        "generator_hashes": generator_hashes,
        "input_rows": input_rows,
    }


def write_artifacts(output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    result = build()

    program_path = output_directory / "j1_factor_program.csv"
    blocks_path = output_directory / "j1_coefficient_blocks.npz"
    verification_path = output_directory / "j1_endpoint_verification.csv"
    input_path = output_directory / "input_hashes.csv"
    pending_path = output_directory / "pending_j2_j3_shamash_contract.json"
    certificate_path = output_directory / "j1_overlap_lift_certificate.json"

    write_csv(program_path, result["program_rows"])
    np.savez_compressed(blocks_path, **result["coefficient_blocks"])
    write_csv(verification_path, result["verification_rows"])
    write_csv(input_path, result["input_rows"])

    pending = {
        "schema": "section12-HS-J2-J3-Shamash-comparison-contract-v1",
        "field": "F_31",
        "closed_here": {
            "J0": "all 384 original chart-to-edge quotient restrictions",
            "J1": (
                "all 384 exact graph-Koszul Laurent lifts; "
                "d_edge J1=J0 d_chart for all 3072 graph generators"
            ),
        },
        "graph_exterior_observation": {
            "wedge2_and_wedge3_are_formally_determined_by_J1": True,
            "formulas": [
                "J2_K=(Lambda^2 D)(I_28 tensor J0)",
                "J3_K=(Lambda^3 D)(I_56 tensor J0)",
            ],
            "why_not_promoted_to_full_J2_or_J3": (
                "P2=K2 plus six K0 equation sectors and P3=K3 plus six K1 "
                "equation sectors; exterior powers cover only K2 and K3"
            ),
        },
        "known_regular_sequence_change": {
            "notation": "t=v_source^-1 in the Laurent overlap ring",
            "B": (
                "diag(1,1,1,1,t,t), with the changed curve entry among the "
                "first four replaced by t^3"
            ),
            "ordering": [
                "curve_0",
                "curve_1",
                "curve_2",
                "curve_3",
                "theta_0",
                "theta_1",
            ],
            "identities": [
                "f_target_changed_curve(u/v,1/v)=t^3*f_source_changed_curve",
                "f_target_theta_j(pulled back)=t*f_source_theta_j for j=0,1",
            ],
            "warning": (
                "Laurent t is not the finite quotient multiplication operator "
                "T=M_v^-1"
            ),
        },
        "next_exact_operation": {
            "bind_prequotient_data": (
                "emit and SHA-bind the six source equations, six pulled target "
                "equations, diagonal B identity, and source/target/edge H0,H1 "
                "divided-difference blocks on every edge"
            ),
            "degree2_comparison": (
                "Delta_j=J1_K*H0_target_j-sum_l(H0_edge_l*B_lj*J0); "
                "solve d2*Gamma_j=Delta_j and emit "
                "J2=[[J2_K,Gamma],[0,B tensor J0]]"
            ),
            "degree3_coherence": (
                "Xi_j=J2_K*H1_target_j+Gamma_j*d1_target-"
                "sum_l(H1_edge_l*B_lj*J1_K); solve d3*Omega_j=Xi_j and "
                "emit J3=[[J3_K,Omega],[0,B tensor J1_K]]"
            ),
        },
        "not_available_from_the_current_quotient_artifacts_alone": [
            "the six pre-quotient surface polynomials and their Laurent supports",
            "the 1152 degree-two comparison homotopies Gamma_j",
            "the 1152 degree-three comparison homotopies Omega_j",
            "the square/cube coherence of those homotopies on higher overlaps",
        ],
        "not_emitted": [
            "full J2",
            "full J3",
            "transported HS lambda",
            "global extension cochain",
            "V1",
            "C",
            "rho_U",
            "rho_F",
            "Yoneda product",
            "T",
        ],
        "missing_maps_are_not_zero": True,
    }
    write_json_atomic(pending_path, pending)

    endpoint_counts = Counter(
        row["endpoint"] for row in result["verification_rows"]
    )
    term_counts = Counter(row["endpoint"] for row in result["program_rows"])
    certificate = {
        "schema": "section12-HS-graph-Koszul-J1-overlap-lifts-v1",
        "field": "F_31",
        "status": "all-384-exact-J1-lifts-emitted-and-verified;full-J2-J3-pending",
        "coordinate_transition": {
            "edge_coordinates": "source chart coordinates",
            "changed_factor": ["u_target=u_source/v_source", "v_target=1/v_source"],
            "unchanged_factors": "identity",
            "denominator": "v_source, invertible on the overlap",
        },
        "counts": {
            "surfaces": 6,
            "edges": 192,
            "endpoint_maps": len(result["verification_rows"]),
            "source_endpoint_maps": endpoint_counts["source"],
            "target_endpoint_maps": endpoint_counts["target"],
            "graph_generators_per_map": 8,
            "graph_generator_chain_identities": len(result["generator_hashes"]),
            "factored_J1_terms": len(result["program_rows"]),
            "source_factored_terms": term_counts["source"],
            "target_factored_terms": term_counts["target"],
            "stored_coefficient_blocks": len(result["coefficient_blocks"]),
            "input_files_bound": len(result["input_rows"]),
        },
        "exact_checks": {
            "all_edges_change_exactly_one_0_to_1_factor": True,
            "all_J0_shapes_match_chart_and_edge_quotient_dimensions": True,
            "all_3072_Laurent_graph_generator_residuals_are_zero": True,
            "all_changed_u_blocks_use_the_exact_divided_difference": True,
            "all_changed_v_blocks_use_the_exact_inverse_difference": True,
            "all_576_inputs_are_bound_by_SHA256": True,
        },
        "storage": {
            "representation": "factored Laurent block program",
            "base_blocks": "one stored J0 matrix per endpoint",
            "derived_blocks": (
                "only J0*M_target_u and J0*M_target_v on target endpoints"
            ),
            "dense_8n_by_8n_matrices_emitted": False,
            "program_is_fully_executable": True,
        },
        "claim_boundary": {
            "proved": (
                "the exact resolution-level graph-Koszul J1 restriction on all "
                "384 endpoints of the original 192-edge product atlas"
            ),
            "not_claimed": (
                "full Shamash J2/J3, transported lambda, a global HS extension, "
                "V1, C, an endpoint comparison, Yoneda, T, Question 11.4, or Hodge"
            ),
        },
        "files": {
            "J1_program": program_path.name,
            "coefficient_blocks": blocks_path.name,
            "endpoint_verification": verification_path.name,
            "input_hashes": input_path.name,
            "pending_J2_J3_contract": pending_path.name,
        },
    }
    write_json_atomic(certificate_path, certificate)

    manifest_files = {}
    for path in (
        program_path,
        blocks_path,
        verification_path,
        input_path,
        pending_path,
        certificate_path,
    ):
        manifest_files[path.name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    manifest = {
        "schema": "section12-HS-graph-Koszul-J1-overlap-lifts-manifest-v1",
        "generator": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "files": manifest_files,
    }
    write_json_atomic(output_directory / "manifest.json", manifest)
    return certificate, result


def main():
    certificate, _result = write_artifacts()
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
