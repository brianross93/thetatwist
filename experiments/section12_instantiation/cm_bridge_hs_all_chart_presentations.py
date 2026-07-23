"""Apply one Hartshorne--Serre presentation law to every affine chart.

The 96 persisted Hermite charts have different finite lengths, so their
matrix sizes are not constant.  This producer derives every size from the
selected Hermite artifact and applies the same graph--Koszul, Shamash,
tree-contraction, and local-duality construction used by the representative
``Z_2/C0010`` computation.

The default run materializes the independent second-chart pilot
``Z_2/C0000``.  ``--all`` applies the same function sequentially to the full
six-by-sixteen chart census.  The inventory artifact binds all 96 available
inputs even when only the pilot is materialized.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import sys
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (
    cm_bridge_hs_extension_from_surface_presentation as extension,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_chart_local_syzygy_audit as chart_context,
)
from experiments.section12_instantiation import (
    cm_bridge_iw_presentation_solver as iw_solver,
)
from experiments.section12_instantiation import (
    cm_bridge_iw_surface_ring_presentation as surface_presentation,
)


P = 31
VARIABLE_COUNT = 8
PILOT_SURFACE = 2
PILOT_CHART = "C0000"
BORDER_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_iw_hermite_border_basis"
)
ORIENTATION_PATH = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_border_to_local_frames/"
    "hs_completed_lci_supports.csv"
)
REGULAR_PATH = ROOT / (
    "results/section12_instantiation/cm_bridge_o8_section_ring_extension/"
    "o8_section_ring_extension_certificate.json"
)
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_all_chart_presentations"
)
SURFACE_EQUATION_INPUTS = {
    "global_chart_certificate": chart_context.GLOBAL_CERTIFICATE,
    "global_degree2_generators": chart_context.GLOBAL_GENERATORS,
    "product_theta_atlas": chart_context.ATLAS,
    "product_theta_atlas_edges": chart_context.ATLAS_EDGES,
    "six_HS_seeds": chart_context.SEEDS,
    "six_HS_support_orbits": chart_context.SUPPORTS,
    "K_descent_certificate": chart_context.DESCENT,
    "theta_coefficients": chart_context.THETA_COEFFICIENTS,
}


class HSAllChartPresentationError(RuntimeError):
    """A bound chart, orientation, or common presentation identity changed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def chart_names():
    return [f"C{value:04b}" for value in range(16)]


def chart_paths(surface: int, chart: str):
    stem = f"Z_{int(surface)}_{chart}"
    return BORDER_DIRECTORY / f"{stem}.npz", BORDER_DIRECTORY / f"{stem}.json"


def surface_equations(surface: int, chart: str):
    context = iw_solver.load_problem_context()
    if chart not in context["chart_bits"]:
        raise HSAllChartPresentationError(f"unknown affine chart {chart}")
    if not 0 <= int(surface) < len(context["seeds"]):
        raise HSAllChartPresentationError(f"unknown surface Z_{surface}")
    bits = context["chart_bits"][chart]
    labels = list(
        map(
            int,
            json.loads(
                context["seeds"][int(surface)]["surface_section_labels"]
            ),
        )
    )
    equations = [
        *iw_solver.curve_polynomials(bits),
        *[
            iw_solver.theta_polynomial(context["theta"][label], bits)
            for label in labels
        ],
    ]
    if len(equations) != 6:
        raise HSAllChartPresentationError("a surface lost its six equations")
    return equations, labels


def transition_quotient(data):
    """Return the graph-to-border quotient and its triangular tree."""

    standard = data["standard"]
    border = data["border"]
    standard_index = {value: index for index, value in enumerate(standard)}
    border_index = {value: index for index, value in enumerate(border)}
    quotient = -np.ones((VARIABLE_COUNT, len(standard)), dtype=np.int64)
    internal = 0
    border_transitions = 0
    for variable in range(VARIABLE_COUNT):
        for source, exponent in enumerate(standard):
            product = list(exponent)
            product[variable] += 1
            product = tuple(product)
            if product in standard_index:
                internal += 1
            elif product in border_index:
                quotient[variable, source] = border_index[product]
                border_transitions += 1
            else:
                raise HSAllChartPresentationError(
                    "a graph transition is neither standard nor border"
                )
    if internal + border_transitions != VARIABLE_COUNT * len(standard):
        raise HSAllChartPresentationError("the graph-transition census is incomplete")
    if set(map(int, quotient[quotient >= 0])) != set(range(len(border))):
        raise HSAllChartPresentationError("the graph-to-border map is not onto")

    tree = []
    for target, exponent in enumerate(standard):
        if exponent == (0,) * VARIABLE_COUNT:
            continue
        for variable, power in enumerate(exponent):
            if not power:
                continue
            parent = list(exponent)
            parent[variable] -= 1
            parent = tuple(parent)
            if parent in standard_index:
                tree.append((variable, standard_index[parent], target))
                break
        else:
            raise HSAllChartPresentationError("the order ideal lost a parent")
    if len(tree) != len(standard) - 1:
        raise HSAllChartPresentationError("the triangular tree has the wrong size")
    return quotient, tree, internal, border_transitions


def orientation_map(surface: int):
    rows = {}
    with ORIENTATION_PATH.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["surface"] != f"Z_{int(surface)}":
                continue
            key = tuple(map(int, json.loads(row["support_point_indices"])))
            rows[key] = (
                int(row["delta0_mod_31"]) % P,
                int(row["delta1_mod_31"]) % P,
            )
    if len(rows) != 81:
        raise HSAllChartPresentationError(
            f"Z_{surface} has {len(rows)} rather than 81 orientation rows"
        )
    return rows


def load_chart_inputs(surface: int, chart: str):
    binary_path, summary_path = chart_paths(surface, chart)
    for path in (
        binary_path,
        summary_path,
        ORIENTATION_PATH,
        REGULAR_PATH,
        *SURFACE_EQUATION_INPUTS.values(),
    ):
        if not path.is_file():
            raise HSAllChartPresentationError(f"missing exact input: {path}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    with np.load(binary_path) as artifact:
        field = int(artifact["field"][0])
        standard = np.asarray(artifact["standard_monomials"], dtype=np.int64)
        multiplication = np.asarray(
            artifact["multiplication_matrices"], dtype=np.int64
        ) % P
        border = np.asarray(artifact["border_monomials"], dtype=np.int64)
        normal_forms = np.asarray(
            artifact["border_normal_form_coefficients"], dtype=np.int64
        ) % P
        support_indices = np.asarray(
            artifact["support_point_indices"], dtype=np.int64
        )
        affine_values = np.asarray(artifact["affine_values"], dtype=np.int64) % P
        affine_tangents = (
            np.asarray(artifact["affine_tangents"], dtype=np.int64) % P
        )
        hermite = np.asarray(artifact["hermite_matrix"], dtype=np.int64) % P
        hermite_inverse = (
            np.asarray(artifact["hermite_inverse"], dtype=np.int64) % P
        )
    dimension = int(standard.shape[0])
    support_count = int(support_indices.shape[0])
    border_rank = int(border.shape[0])
    summary_counts = summary["summary"]
    if (
        int(summary_counts["quotient_length"]) != dimension
        or int(summary_counts["support_points"]) != support_count
        or int(summary_counts["border_monomials"]) != border_rank
    ):
        raise HSAllChartPresentationError(
            "the Hermite binary and summary dimensions disagree"
        )
    if field != P or dimension != 2 * support_count:
        raise HSAllChartPresentationError("the dual-number chart length changed")
    if standard.shape != (dimension, VARIABLE_COUNT):
        raise HSAllChartPresentationError("the standard basis shape changed")
    if multiplication.shape != (VARIABLE_COUNT, dimension, dimension):
        raise HSAllChartPresentationError("the multiplication shape changed")
    if support_indices.shape != (support_count, 4):
        raise HSAllChartPresentationError("the support-index shape changed")
    if affine_values.shape != (support_count, VARIABLE_COUNT):
        raise HSAllChartPresentationError("the affine-value shape changed")
    if affine_tangents.shape != (support_count, VARIABLE_COUNT):
        raise HSAllChartPresentationError("the affine-tangent shape changed")
    if border.shape != (border_rank, VARIABLE_COUNT):
        raise HSAllChartPresentationError("the border basis shape changed")
    if normal_forms.shape != (border_rank, dimension):
        raise HSAllChartPresentationError("the border normal forms changed")
    if hermite.shape != (dimension, dimension) or hermite_inverse.shape != (
        dimension,
        dimension,
    ):
        raise HSAllChartPresentationError("the Hermite matrices changed shape")
    identity = np.eye(dimension, dtype=np.int64)
    if np.any((hermite @ hermite_inverse - identity) % P) or np.any(
        (hermite_inverse @ hermite - identity) % P
    ):
        raise HSAllChartPresentationError("the Hermite basis is not invertible")

    equations, labels = surface_equations(surface, chart)
    for equation in equations:
        if np.any(
            surface_presentation._evaluate_matrix_polynomial(
                equation, multiplication
            )
        ):
            raise HSAllChartPresentationError(
                "a surface equation does not annihilate the Hermite algebra"
            )
    data = {
        "path": binary_path,
        "summary_path": summary_path,
        "standard": [tuple(map(int, row)) for row in standard],
        "multiplication": multiplication,
        "border": [tuple(map(int, row)) for row in border],
        "normal_forms": normal_forms,
        "surface_equations": equations,
        "surface_labels": labels,
        "regular_path": REGULAR_PATH,
    }
    quotient, tree, internal, border_transitions = transition_quotient(data)
    data["_transition_quotient"] = quotient

    orientations = orientation_map(surface)
    tau_on_jets = np.zeros(dimension, dtype=np.int64)
    local_rows = []
    for position, support in enumerate(support_indices):
        key = tuple(map(int, support))
        if key not in orientations:
            raise HSAllChartPresentationError(
                f"support {key} has no Z_{surface} orientation"
            )
        delta0, delta1 = orientations[key]
        if not delta0:
            raise HSAllChartPresentationError("a local orientation is not a unit")
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
    return {
        "surface": int(surface),
        "chart": chart,
        "summary": summary,
        "data": data,
        "quotient": quotient,
        "tree": tree,
        "internal": internal,
        "border_transitions": border_transitions,
        "support_indices": support_indices,
        "hermite": hermite,
        "hermite_inverse": hermite_inverse,
        "affine_values": affine_values,
        "affine_tangents": affine_tangents,
        "tau_on_jets": tau_on_jets,
        "tau": tau_on_jets @ hermite % P,
        "local_rows": local_rows,
    }


@contextmanager
def extension_dimension(dimension: int):
    old = extension.DIMENSION
    extension.DIMENSION = int(dimension)
    try:
        yield
    finally:
        extension.DIMENSION = old


def build_chart(surface=PILOT_SURFACE, chart=PILOT_CHART):
    started = time.monotonic()
    inputs = load_chart_inputs(surface, chart)
    dimension = len(inputs["data"]["standard"])
    with extension_dimension(dimension):
        contraction = extension.build_tree_contraction(
            inputs["data"], inputs["tree"]
        )
        d2_tree, unit_cycles = extension.build_d2_tree(
            inputs["data"], contraction
        )
        koszul_d3 = extension.build_koszul_d3(inputs["data"])
        h_operators = extension.build_h_operators(inputs["data"])
        h1_raw = extension.build_h1_raw(h_operators)
        d3_factor = extension.build_tree_d3_factor(
            koszul_d3, h1_raw, contraction, contraction["generators"]
        )
        cocycle = extension.build_change_of_rings_cocycle(inputs)
        verification = extension.verify_shamash_and_cocycle(
            inputs["data"],
            d2_tree,
            d3_factor,
            cocycle,
            contraction,
            h_operators,
        )
        relation_bottom = extension.build_module_presentation(d2_tree)

        merged = {}
        for (row, column), polynomial in d2_tree.items():
            for target in np.flatnonzero(contraction["q"][:, row]):
                extension._matrix_entry_add(
                    merged,
                    int(target),
                    int(column),
                    polynomial,
                    int(contraction["q"][int(target), row]),
                )
        expected_columns = [
            *surface_presentation.koszul_columns(
                inputs["data"], inputs["quotient"]
            ),
            *surface_presentation.surface_columns(
                inputs["data"], inputs["quotient"]
            )[0],
        ]
        expected_merged = {
            (int(row), int(column)): dict(polynomial)
            for column, entries in enumerate(expected_columns)
            for row, polynomial in entries.items()
        }
        if merged != expected_merged:
            raise HSAllChartPresentationError(
                "the tree model does not merge to the chart's border model"
            )

    retained_rank = 7 * dimension + 1
    relation_rank = 28 * dimension + 6
    d3_source_rank = 56 * dimension + 6 * retained_rank
    return {
        "inputs": inputs,
        "contraction": contraction,
        "d2_tree": d2_tree,
        "unit_cycles": unit_cycles,
        "d3_factor": d3_factor,
        "h_operators": h_operators,
        "cocycle": cocycle,
        "verification": verification,
        "module_relation_bottom": relation_bottom,
        "merged_d2": merged,
        "dimensions": {
            "finite_algebra": dimension,
            "support_count": dimension // 2,
            "raw_K1": 8 * dimension,
            "tree_F1": retained_rank,
            "K2": 28 * dimension,
            "F2": relation_rank,
            "K3": 56 * dimension,
            "F3_factored": d3_source_rank,
            "module_target": retained_rank + 1,
            "border_carrier": len(inputs["data"]["border"]),
            "merge_kernel": retained_rank - len(inputs["data"]["border"]),
        },
        "elapsed_seconds": time.monotonic() - started,
    }


def inventory_rows():
    rows = []
    for surface, chart in itertools.product(range(6), chart_names()):
        binary, summary_path = chart_paths(surface, chart)
        if not binary.is_file() or not summary_path.is_file():
            raise HSAllChartPresentationError(
                f"the 96-chart input census is missing Z_{surface}/{chart}"
            )
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        dimension = int(summary["summary"]["quotient_length"])
        support_count = int(summary["summary"]["support_points"])
        border_rank = int(summary["summary"]["border_monomials"])
        if dimension != 2 * support_count:
            raise HSAllChartPresentationError("an inventory length is not curvilinear")
        retained_rank = 7 * dimension + 1
        rows.append(
            {
                "surface": surface,
                "chart": chart,
                "finite_algebra_dimension": dimension,
                "support_count": support_count,
                "border_rank": border_rank,
                "tree_F1_rank": retained_rank,
                "F2_rank": 28 * dimension + 6,
                "F3_factored_rank": 56 * dimension + 6 * retained_rank,
                "module_target_rank": retained_rank + 1,
                "merge_kernel_rank": retained_rank - border_rank,
                "binary_path": str(binary.relative_to(ROOT)),
                "binary_sha256": sha256(binary),
                "summary_path": str(summary_path.relative_to(ROOT)),
                "summary_sha256": sha256(summary_path),
            }
        )
    return rows


def write_chart(result, output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    inputs = result["inputs"]
    surface = inputs["surface"]
    chart = inputs["chart"]
    dimension = result["dimensions"]["finite_algebra"]
    contraction = result["contraction"]
    d3 = result["d3_factor"]
    cocycle = result["cocycle"]
    payload = {
        "field": np.asarray([P], dtype=np.uint8),
        "surface": np.asarray([surface], dtype=np.uint8),
        "chart_bits": np.asarray([int(value) for value in chart[1:]], dtype=np.uint8),
        "tree_complex_ranks": np.asarray(
            [
                1,
                result["dimensions"]["tree_F1"],
                result["dimensions"]["F2"],
                result["dimensions"]["F3_factored"],
            ],
            dtype=np.uint32,
        ),
        "module_presentation_shape": np.asarray(
            [
                result["dimensions"]["module_target"],
                result["dimensions"]["F2"],
            ],
            dtype=np.uint32,
        ),
        "support_point_indices": np.asarray(
            inputs["support_indices"], dtype=np.uint8
        ),
        "affine_values": np.asarray(inputs["affine_values"], dtype=np.uint8),
        "affine_tangents": np.asarray(inputs["affine_tangents"], dtype=np.uint8),
        "tau_on_jets": np.asarray(inputs["tau_on_jets"], dtype=np.uint8),
        "tau_on_standard_basis": np.asarray(inputs["tau"], dtype=np.uint8),
        "tree_edges": np.asarray(contraction["tree_edges"], dtype=np.uint32),
        "non_tree_edges": np.asarray(contraction["non_tree_edges"], dtype=np.uint32),
        "merge_q": np.asarray(contraction["q"], dtype=np.uint8),
        "merge_section": np.asarray(contraction["section"], dtype=np.uint8),
        "merge_kernel": np.asarray(contraction["kernel"], dtype=np.uint8),
        "merge_kernel_left_inverse": np.asarray(
            contraction["kernel_left_inverse"], dtype=np.uint8
        ),
        "tree_projection_p1": np.asarray(contraction["p1"], dtype=np.uint8),
        "hermite_standard_to_jets": np.asarray(cocycle["hermite"], dtype=np.uint8),
        "hermite_jets_to_standard": np.asarray(
            cocycle["hermite_inverse"], dtype=np.uint8
        ),
    }
    with extension_dimension(dimension):
        extension._store_matrix(payload, "D2_tree", result["d2_tree"])
        extension._store_matrix(
            payload, "D3_K3_to_K2", d3["koszul_K3_to_K2"]
        )
        extension._store_matrix(
            payload,
            "tree_K1_inclusion_i1",
            extension.i1_matrix(contraction["i1"]),
        )
        for equation, matrix in enumerate(d3["H1_raw_K1_to_K2"]):
            extension._store_matrix(
                payload, f"D3_H1_raw_equation_{equation}", matrix
            )
        lambda_matrix = extension.lambda_jet_matrix(
            cocycle["local_lambda_K2"]
        )
        extension._store_matrix(payload, "lambda_K2_jet", lambda_matrix)
        extension._store_matrix(
            payload,
            "Vprime_relation_bottom_minus_D2",
            result["module_relation_bottom"],
        )
        extension._store_matrix(payload, "merged_D2_border", result["merged_d2"])

    artifact_path = output_directory / f"Z_{surface}_{chart}_hs_extension.npz"
    np.savez_compressed(artifact_path, **payload)
    certificate = {
        "schema": "section12-all-chart-HS-presentation-v1",
        "field": "F_31",
        "surface": f"Z_{surface}",
        "chart": chart,
        "common_law": (
            "graph Koszul plus six divided-difference Shamash homotopies, "
            "triangular order-tree contraction, and supportwise local-duality lambda"
        ),
        "dimensions": result["dimensions"],
        "orientation_binding": {
            "rows": len(inputs["local_rows"]),
            "all_support_indices_found_in_surface_ledger": True,
            "all_delta0_nonzero": all(row["delta0"] for row in inputs["local_rows"]),
        },
        "change_of_rings_cocycle": {
            "supportwise_dual_number_factors": len(cocycle["local_lambda_K2"]),
            "mu_factorization": cocycle["mu_factorization"],
        },
        "exact_checks": {
            "Hermite_basis_invertible": True,
            "six_surface_equations_annihilate_A_W": True,
            "tree_inclusion_and_projection_split": True,
            "D2_D3_zero_over_surface_ring": True,
            "lambda_D3_zero_over_surface_ring": True,
            "A_D3_zero_over_surface_ring": True,
            "lambda_is_closed_Ext2_cochain": True,
            "tree_model_merges_to_independent_border_model": True,
            "C0010_curve_reduction_helpers_not_used_by_this_build": True,
        },
        "verification_ledger": result["verification"],
        "elapsed_seconds": result["elapsed_seconds"],
        "artifact": {
            "path": artifact_path.name,
            "bytes": artifact_path.stat().st_size,
            "sha256": sha256(artifact_path),
        },
        "upstream": {
            "Hermite_binary": {
                "path": str(inputs["data"]["path"].relative_to(ROOT)),
                "sha256": sha256(inputs["data"]["path"]),
            },
            "Hermite_summary": {
                "path": str(inputs["data"]["summary_path"].relative_to(ROOT)),
                "sha256": sha256(inputs["data"]["summary_path"]),
            },
            "orientation_ledger": {
                "path": str(ORIENTATION_PATH.relative_to(ROOT)),
                "sha256": sha256(ORIENTATION_PATH),
            },
            "regular_sequence_certificate": {
                "path": str(REGULAR_PATH.relative_to(ROOT)),
                "sha256": sha256(REGULAR_PATH),
            },
            **{
                name: {
                    "path": str(path.relative_to(ROOT)),
                    "sha256": sha256(path),
                }
                for name, path in SURFACE_EQUATION_INPUTS.items()
            },
        },
        "claim_boundary": {
            "proved": (
                "one exact parameterized chart presentation with factored D3 and "
                "lambda, obtained without a chart-specific rule"
            ),
            "not_claimed": (
                "materialization of every other chart, inter-chart transition maps, "
                "the independently normalized K-target presentations required by "
                "the K-chain comparison contract, a global projective idempotent, "
                "or a Hodge-theoretic verdict"
            ),
        },
    }
    certificate_path = output_directory / f"Z_{surface}_{chart}_certificate.json"
    certificate_path.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return artifact_path, certificate_path


def write_inventory(output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    rows = inventory_rows()
    path = output_directory / "all_chart_input_inventory.json"
    payload = {
        "schema": "section12-all-chart-HS-input-inventory-v1",
        "chart_count": len(rows),
        "surfaces": 6,
        "charts_per_surface": 16,
        "dimension_law": {
            "tree_F1": "7*N+1",
            "F2": "28*N+6",
            "F3_factored": "56*N+6*(7*N+1)",
            "module_target": "7*N+2",
        },
        "rows": rows,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def reproduce(
    surface=PILOT_SURFACE,
    chart=PILOT_CHART,
    output_directory=OUTPUT_DIRECTORY,
):
    inventory = write_inventory(output_directory)
    result = build_chart(surface, chart)
    artifact, certificate = write_chart(result, output_directory)
    manifest_path = Path(output_directory) / "manifest.json"
    paths = [inventory, artifact, certificate]
    manifest_path.write_text(
        json.dumps(
            {
                path.name: {
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in paths
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return result, {
        "inventory": inventory,
        "artifact": artifact,
        "certificate": certificate,
        "manifest": manifest_path,
    }


def reproduce_all(output_directory=OUTPUT_DIRECTORY):
    inventory = write_inventory(output_directory)
    outputs = []
    for surface, chart in itertools.product(range(6), chart_names()):
        result = build_chart(surface, chart)
        outputs.append(write_chart(result, output_directory))
    return inventory, outputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface", type=int, default=PILOT_SURFACE)
    parser.add_argument("--chart", default=PILOT_CHART)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--output-directory", type=Path, default=OUTPUT_DIRECTORY)
    args = parser.parse_args(argv)
    if args.all:
        inventory, outputs = reproduce_all(args.output_directory)
        print(json.dumps({"inventory": str(inventory), "charts": len(outputs)}, indent=2))
        return 0
    result, paths = reproduce(args.surface, args.chart, args.output_directory)
    print(
        json.dumps(
            {
                "surface": result["inputs"]["surface"],
                "chart": result["inputs"]["chart"],
                "dimensions": result["dimensions"],
                "elapsed_seconds": result["elapsed_seconds"],
                "certificate": str(paths["certificate"]),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
