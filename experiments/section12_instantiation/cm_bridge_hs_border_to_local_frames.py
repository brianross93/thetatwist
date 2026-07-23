"""Convert the Hermite border package into completed-local HS frames.

The border-basis producer gives the exact finite quotient algebras
``O(W_s cap U)`` on all 96 surface charts and their restrictions on all 192
cube edges.  This producer performs the homological conversion that is
already intrinsic at the support points.

For every one of the 6*81 lifted curvilinear points it writes the lci model

    I=(x,y^2),             R --(-y^2,x)^t--> R^2 --> I,

and the oriented Hartshorne--Serre pushout

    R --(Delta,y^2,-x)^t--> R + R^2 --> V1,

where ``Delta=delta_0+delta_1*y`` is a unit.  Consequently ``V1`` is free
with basis inherited from ``R^2``; the line inclusion and quotient are

    (-y^2/Delta,x/Delta)^t,       (x,y^2).

On the punctured overlap the two quotient splittings differ by
``Delta/(x*y^2)``, which is the recorded local Ext orientation.

The producer also verifies that every Hermite chart-to-overlap map preserves
the direct product of point dual-number factors, including both idempotents
and tangent nilpotents.  What it deliberately does not invent is a lift of
those quotient-algebra maps to chain maps over the localized surface rings.
That final lift is emitted as one exact algebraic contract.

No Markman, Question 11.4, or Hodge-conjecture verdict is claimed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_iw_hermite_border_basis as border,
    cm_bridge_iw_presentation_solver as iw_solver,
)


P = 31
BORDER_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_iw_hermite_border_basis"
)
SEED_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter"
)
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_border_to_local_frames"
)


class HSBorderFrameError(RuntimeError):
    """An exact upstream invariant or local-frame identity changed."""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows, fieldnames=None):
    rows = list(rows)
    if not rows and fieldnames is None:
        raise ValueError("fieldnames are required for an empty CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_inputs():
    manifest_path = BORDER_DIRECTORY / "hermite_border_basis_manifest.json"
    manifest = border.verify_manifest(manifest_path)
    seed_rows = read_csv(SEED_DIRECTORY / "six_hs_local_seeds.csv")
    support_rows = read_csv(SEED_DIRECTORY / "six_hs_support_K_orbits.csv")
    cocycles = json.loads(
        (SEED_DIRECTORY / "six_hs_completed_local_cocycles.json").read_text(
            encoding="utf-8"
        )
    )
    if len(seed_rows) != 6 or len(cocycles) != 6 or len(support_rows) != 486:
        raise HSBorderFrameError("the six-surface HS seed package changed")
    if manifest.get("chart_count") != 96 or manifest.get("overlap_count") != 192:
        raise HSBorderFrameError("the complete Hermite border package is absent")
    return manifest_path, manifest, seed_rows, support_rows, cocycles


def support_key(value):
    return tuple(map(int, json.loads(value)))


def support_index(seed_rows, support_rows):
    """Index the 486 K lifts and bind the invariant adapted-coordinate rows."""

    seed_by_surface = {row["surface"]: row for row in seed_rows}
    result = {}
    for row in support_rows:
        surface = row["surface"]
        point = support_key(row["source_point_indices"])
        key = (surface, point)
        if key in result:
            raise HSBorderFrameError("a lifted support appears twice")
        seed = seed_by_surface[surface]
        tangent = json.loads(row["transported_tangent_fixed_frame_coefficients"])
        if tangent != json.loads(seed["tangent_fixed_frame_coefficients"]):
            raise HSBorderFrameError("a transported tangent changed fixed-frame coefficients")
        result[key] = {
            "K_lift_index": int(row["K_orbit_lift_index"]),
            "x_row": tuple(map(int, json.loads(seed["x_fixed_frame_coefficients"]))),
            "y_row": tuple(map(int, json.loads(seed["y_fixed_frame_coefficients"]))),
            "tangent": tuple(map(int, tangent)),
        }
    if len(result) != 486:
        raise HSBorderFrameError("the lifted support index is incomplete")
    return result


def chart_rows_by_key(manifest):
    return {(row["surface"], row["chart"]): row for row in manifest["charts"]}


def chart_artifact(surface: str, chart: str):
    return BORDER_DIRECTORY / f"{surface}_{chart}.npz"


def canonical_support_charts(manifest, support_lookup):
    """Choose the first exact affine chart containing each support lift."""

    choices = {}
    for row in sorted(manifest["charts"], key=lambda item: (item["surface"], item["chart"])):
        with np.load(BORDER_DIRECTORY / row["artifact"]) as stored:
            for raw in stored["support_point_indices"]:
                point = tuple(map(int, raw))
                key = (row["surface"], point)
                if key not in support_lookup:
                    raise HSBorderFrameError("a border support is absent from the K-orbit ledger")
                choices.setdefault(key, row["chart"])
    if set(choices) != set(support_lookup):
        raise HSBorderFrameError("not every K lift occurs in the affine chart cover")
    return choices


def boundary_delta_on_chart(context, surface_index, chart, values, tangent):
    bits = context["chart_bits"][chart]
    labels = ((surface_index - 1) % 6, (surface_index + 2) % 6)
    jets = []
    for label in labels:
        polynomial = iw_solver.theta_polynomial(context["theta"][label], bits)
        jets.append(border.polynomial_jet(polynomial, values, tangent))
    delta0 = jets[0][0] * jets[1][0] % P
    delta1 = (jets[0][0] * jets[1][1] + jets[0][1] * jets[1][0]) % P
    if not delta0:
        raise HSBorderFrameError("a lifted W support met its boundary divisor")
    inverse0 = pow(delta0, -1, P)
    inverse1 = -delta1 * inverse0 * inverse0 % P
    if (delta0 * inverse0) % P != 1 or (
        delta0 * inverse1 + delta1 * inverse0
    ) % P:
        raise ArithmeticError("a boundary dual-number inverse failed")
    return labels, (delta0, delta1), (inverse0, inverse1)


def completed_support_records(manifest, support_lookup, canonical_charts):
    """Compute one chart-trivialized boundary unit for each lifted support."""

    context = iw_solver.load_problem_context()
    records = []
    by_chart = chart_rows_by_key(manifest)
    for key in sorted(support_lookup, key=lambda item: (int(item[0].split("_")[1]), item[1])):
        surface, point = key
        surface_index = int(surface.split("_")[1])
        chart = canonical_charts[key]
        row = by_chart[(surface, chart)]
        with np.load(BORDER_DIRECTORY / row["artifact"]) as stored:
            supports = [tuple(map(int, value)) for value in stored["support_point_indices"]]
            local_index = supports.index(point)
            values = tuple(map(int, stored["affine_values"][local_index]))
            tangent = tuple(map(int, stored["affine_tangents"][local_index]))
        labels, delta, delta_inverse = boundary_delta_on_chart(
            context, surface_index, chart, values, tangent
        )
        support = support_lookup[key]
        records.append(
            {
                "surface": surface,
                "surface_index": surface_index,
                "support_point_indices": json.dumps(list(point)),
                "K_lift_index": support["K_lift_index"],
                "canonical_chart": chart,
                "canonical_chart_support_position": local_index,
                "x_fixed_frame_coefficients": json.dumps(list(support["x_row"])),
                "y_fixed_frame_coefficients": json.dumps(list(support["y_row"])),
                "tangent_fixed_frame_coefficients": json.dumps(list(support["tangent"])),
                "boundary_section_labels": json.dumps(list(labels)),
                "delta0_mod_31": delta[0],
                "delta1_mod_31": delta[1],
                "delta_inverse0_mod_31": delta_inverse[0],
                "delta_inverse1_mod_31": delta_inverse[1],
                "residue_row_on_1_epsilon": json.dumps([delta[1], delta[0]]),
            }
        )
    if len(records) != 486:
        raise HSBorderFrameError("the completed support record count changed")
    return records


FRAME_LAYERS = (
    ("I_resolution_degree_1", ("r",)),
    ("I_resolution_degree_0", ("q_x", "q_y2")),
    ("HS_pushout_source", ("rho",)),
    ("HS_pushout_target", ("a", "q_x", "q_y2")),
    ("V1_free_frame", ("v_x", "v_y2")),
    ("V1_frame_Dx", ("a", "sigma_x")),
    ("V1_frame_Dy", ("a", "sigma_y")),
)


def frame_basis_rows(records):
    rows = []
    for record in records:
        support_id = f"{record['surface']}:K{record['K_lift_index']:02d}"
        for layer, labels in FRAME_LAYERS:
            for index, label in enumerate(labels):
                rows.append(
                    {
                        "support_id": support_id,
                        "surface": record["surface"],
                        "K_lift_index": record["K_lift_index"],
                        "support_point_indices": record["support_point_indices"],
                        "canonical_chart": record["canonical_chart"],
                        "ring": "F_31[[x,y]][Delta^-1]",
                        "layer": layer,
                        "basis_index": index,
                        "basis_label": label,
                    }
                )
    return rows


def frame_map_rows(records):
    """Emit sparse exact matrices in the intrinsic completed-local basis."""

    rows = []

    def add(record, name, kind, source, target, row, column, expression):
        rows.append(
            {
                "support_id": f"{record['surface']}:K{record['K_lift_index']:02d}",
                "map_id": name,
                "map_kind": kind,
                "source_layer": source,
                "target_layer": target,
                "row": row,
                "column": column,
                "coefficient_expression": expression,
                "field": "F_31",
            }
        )

    for record in records:
        d0 = int(record["delta0_mod_31"])
        d1 = int(record["delta1_mod_31"])
        delta = f"({d0}+{d1}*y)"
        prefix = f"{record['surface']}_K{record['K_lift_index']:02d}"

        # R -> R^2 -> I_W.
        add(record, f"{prefix}:I_d1", "lci_resolution", "I_resolution_degree_1", "I_resolution_degree_0", 0, 0, "-y^2")
        add(record, f"{prefix}:I_d1", "lci_resolution", "I_resolution_degree_1", "I_resolution_degree_0", 1, 0, "x")

        # Pushout presentation R -> R + R^2.
        add(record, f"{prefix}:HS_pushout", "HS_pushout_presentation", "HS_pushout_source", "HS_pushout_target", 0, 0, delta)
        add(record, f"{prefix}:HS_pushout", "HS_pushout_presentation", "HS_pushout_source", "HS_pushout_target", 1, 0, "y^2")
        add(record, f"{prefix}:HS_pushout", "HS_pushout_presentation", "HS_pushout_source", "HS_pushout_target", 2, 0, "-x")

        # After Delta is inverted, the quotient is free on q_x,q_y2.
        add(record, f"{prefix}:A_inclusion", "HS_exact_sequence", "A_line", "V1_free_frame", 0, 0, f"-y^2/{delta}")
        add(record, f"{prefix}:A_inclusion", "HS_exact_sequence", "A_line", "V1_free_frame", 1, 0, f"x/{delta}")
        add(record, f"{prefix}:I_projection", "HS_exact_sequence", "V1_free_frame", "I_W", 0, 0, "x")
        add(record, f"{prefix}:I_projection", "HS_exact_sequence", "V1_free_frame", "I_W", 0, 1, "y^2")

        # Quotient splittings on D(x), D(y), expressed in the free V1 frame.
        add(record, f"{prefix}:sigma_x", "punctured_splitting", "I_W_unit", "V1_free_frame", 0, 0, "x^-1")
        add(record, f"{prefix}:sigma_y", "punctured_splitting", "I_W_unit", "V1_free_frame", 1, 0, "y^-2")

        # Cech transition from (a,sigma_x) to (a,sigma_y) and its inverse.
        add(record, f"{prefix}:Dx_to_Dy", "completed_local_transition", "V1_frame_Dx", "V1_frame_Dy", 0, 0, "1")
        add(record, f"{prefix}:Dx_to_Dy", "completed_local_transition", "V1_frame_Dx", "V1_frame_Dy", 0, 1, f"{delta}/(x*y^2)")
        add(record, f"{prefix}:Dx_to_Dy", "completed_local_transition", "V1_frame_Dx", "V1_frame_Dy", 1, 1, "1")
        add(record, f"{prefix}:Dy_to_Dx", "completed_local_transition_inverse", "V1_frame_Dy", "V1_frame_Dx", 0, 0, "1")
        add(record, f"{prefix}:Dy_to_Dx", "completed_local_transition_inverse", "V1_frame_Dy", "V1_frame_Dx", 0, 1, f"-{delta}/(x*y^2)")
        add(record, f"{prefix}:Dy_to_Dx", "completed_local_transition_inverse", "V1_frame_Dy", "V1_frame_Dx", 1, 1, "1")
    return rows


def verify_chart_point_decomposition(manifest):
    rows = []
    for raw in manifest["charts"]:
        path = BORDER_DIRECTORY / raw["artifact"]
        with np.load(path) as stored:
            hermite = stored["hermite_matrix"].astype(np.int64) % P
            inverse = stored["hermite_inverse"].astype(np.int64) % P
            idempotents = stored["point_idempotent_coefficients"].astype(np.int64) % P
            nilpotents = stored["point_nilpotent_coefficients"].astype(np.int64) % P
            dimension = hermite.shape[0]
            support_count = idempotents.shape[1]
            identity = np.eye(dimension, dtype=np.int64)
            if np.any((hermite @ inverse - identity) % P):
                raise HSBorderFrameError("a chart Hermite inverse changed")
            expected_idempotents = identity[:, 0::2]
            expected_nilpotents = identity[:, 1::2]
            if np.any((hermite @ idempotents - expected_idempotents) % P):
                raise HSBorderFrameError("a point idempotent lost its value/jet coordinates")
            if np.any((hermite @ nilpotents - expected_nilpotents) % P):
                raise HSBorderFrameError("a tangent nilpotent lost its value/jet coordinates")
            rows.append(
                {
                    "surface": raw["surface"],
                    "chart": raw["chart"],
                    "support_points": support_count,
                    "quotient_dimension": dimension,
                    "standard_monomial_to_point_dual_map": "hermite_matrix",
                    "point_dual_to_standard_monomial_map": "hermite_inverse",
                    "point_idempotent_columns": "hermite_inverse[:,0::2]",
                    "tangent_nilpotent_columns": "hermite_inverse[:,1::2]",
                    "binary_artifact": str(path.relative_to(ROOT)),
                    "binary_artifact_sha256": sha256(path),
                    "two_sided_coordinate_check": True,
                }
            )
    return rows


def selected_jet_rows(hermite, positions):
    row_indices = [index for position in positions for index in (2 * position, 2 * position + 1)]
    return hermite[row_indices]


def verify_overlap_point_decomposition(manifest):
    chart_lookup = chart_rows_by_key(manifest)
    rows = []
    for raw in manifest["overlaps"]:
        surface = raw["surface"]
        source_chart = raw["source_chart"]
        target_chart = raw["target_chart"]
        source_path = BORDER_DIRECTORY / chart_lookup[(surface, source_chart)]["artifact"]
        target_path = BORDER_DIRECTORY / chart_lookup[(surface, target_chart)]["artifact"]
        overlap_path = BORDER_DIRECTORY / raw["artifact"]
        with np.load(source_path) as source, np.load(target_path) as target, np.load(overlap_path) as overlap:
            source_support = [tuple(map(int, value)) for value in source["support_point_indices"]]
            target_support = [tuple(map(int, value)) for value in target["support_point_indices"]]
            common_support = [tuple(map(int, value)) for value in overlap["support_point_indices"]]
            source_positions = [source_support.index(point) for point in common_support]
            target_positions = [target_support.index(point) for point in common_support]
            hs = source["hermite_matrix"].astype(np.int64) % P
            ht = target["hermite_matrix"].astype(np.int64) % P
            he = overlap["hermite_matrix"].astype(np.int64) % P
            rs = overlap["source_restriction"].astype(np.int64) % P
            rt = overlap["target_restriction"].astype(np.int64) % P
            if np.any((he @ rs - selected_jet_rows(hs, source_positions)) % P):
                raise HSBorderFrameError("a source overlap map changed point-dual coordinates")
            if np.any((he @ rt - selected_jet_rows(ht, target_positions)) % P):
                raise HSBorderFrameError("a target overlap map changed point-dual coordinates")
            edge_inverse = border.inverse_mod(he)
            source_idempotents = source["point_idempotent_coefficients"].astype(np.int64)[:, source_positions] % P
            target_idempotents = target["point_idempotent_coefficients"].astype(np.int64)[:, target_positions] % P
            source_nilpotents = source["point_nilpotent_coefficients"].astype(np.int64)[:, source_positions] % P
            target_nilpotents = target["point_nilpotent_coefficients"].astype(np.int64)[:, target_positions] % P
            edge_idempotents = edge_inverse[:, 0::2]
            edge_nilpotents = edge_inverse[:, 1::2]
            if np.any((rs @ source_idempotents - edge_idempotents) % P) or np.any(
                (rt @ target_idempotents - edge_idempotents) % P
            ):
                raise HSBorderFrameError("an overlap restriction changed a point idempotent")
            if np.any((rs @ source_nilpotents - edge_nilpotents) % P) or np.any(
                (rt @ target_nilpotents - edge_nilpotents) % P
            ):
                raise HSBorderFrameError("an overlap restriction changed a tangent nilpotent")
            rows.append(
                {
                    "surface": surface,
                    "edge_id": int(raw["edge_id"]),
                    "source_chart": source_chart,
                    "target_chart": target_chart,
                    "common_support_points": len(common_support),
                    "overlap_quotient_dimension": he.shape[0],
                    "source_jet_square_commutes": True,
                    "target_jet_square_commutes": True,
                    "point_idempotents_preserved": True,
                    "tangent_nilpotents_preserved": True,
                    "binary_artifact": str(overlap_path.relative_to(ROOT)),
                    "binary_artifact_sha256": sha256(overlap_path),
                }
            )
    return rows


def verify_symbolic_frame_identities(records):
    # Every identity is coefficient-independent except invertibility of Delta,
    # which is exactly delta_0 != 0 in F_31[[x,y]].
    if any(int(record["delta0_mod_31"]) % P == 0 for record in records):
        raise HSBorderFrameError("a completed-local Delta is not a unit")
    return {
        "lci_complex_q_times_d": "x*(-y^2)+y^2*x=0",
        "pushout_elimination": "Delta*A+y^2*q_x-x*q_y2=0",
        "free_frame_inclusion": "Delta*j_A=(-y^2,x)^t",
        "exact_sequence_composition": "(x,y^2)*(-y^2/Delta,x/Delta)^t=0",
        "punctured_splitting_difference": (
            "sigma_y-sigma_x=(-x^-1,y^-2)^t="
            "Delta/(x*y^2)*(-y^2/Delta,x/Delta)^t"
        ),
        "transition_inverse": "[[1,c],[0,1]]*[[1,-c],[0,1]]=I_2",
        "support_instances_checked": len(records),
    }


def missing_chain_lift_contract():
    return {
        "schema": "section12-HS-surface-ring-edge-chain-lift-contract-v1",
        "target": "surface_ring_lci_edge_lifts.coo.csv",
        "scope": "all 192 oriented surface-chart cube edges, before K transport",
        "known_on_every_edge": [
            "the source and target Hermite border presentations of O_W",
            "the exact quotient-algebra restriction maps to the common overlap",
            "the preserved point idempotents and tangent nilpotents",
            "the intrinsic completed-local Koszul and HS pushout matrices",
            "the two line-bundle diagonal theta transitions and determinant transition",
        ],
        "single_missing_conversion": (
            "lift each quotient-algebra overlap square to a syzygy-compatible chain "
            "map between the localized surface-ring lci resolutions"
        ),
        "required_identity": {
            "rings": "R_e is the localized surface ring on the edge overlap",
            "source_quotient_row": "q_s=(x_s,y_s^2)",
            "target_quotient_row": "q_t=(x_t,y_t^2)",
            "source_syzygy_column": "d_s=(-y_s^2,x_s)^t",
            "target_syzygy_column": "d_t=(-y_t^2,x_t)^t",
            "unknowns": "U_e in GL_2(R_e), lambda_e in R_e^times",
            "equations": ["q_s*U_e=q_t", "U_e*d_t=d_s*lambda_e"],
        },
        "why_border_quotients_do_not_choose_the_lift": (
            "the Hermite maps live after quotienting by I_W; both q_s and q_t are "
            "zero there, so their higher surface-ring re-expression is not encoded"
        ),
        "completion_after_this_file": [
            "push out each chain map along the emitted Delta orientation",
            "transport the 192 base lifts through the exact 81-element K action",
            "emit the reserved hs_localized_frame_basis.csv and hs_localized_frame_maps.coo.csv",
            "assemble the existing global HS hyper-Cech solve",
        ],
        "no_new_geometric_parameter": True,
    }


def write_artifacts(output_directory: Path = OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    manifest_path, manifest, seed_rows, support_rows, _cocycles = load_inputs()
    support_lookup = support_index(seed_rows, support_rows)
    canonical_charts = canonical_support_charts(manifest, support_lookup)
    records = completed_support_records(manifest, support_lookup, canonical_charts)
    basis_rows = frame_basis_rows(records)
    map_rows = frame_map_rows(records)
    chart_rows = verify_chart_point_decomposition(manifest)
    overlap_rows = verify_overlap_point_decomposition(manifest)
    symbolic_checks = verify_symbolic_frame_identities(records)
    missing = missing_chain_lift_contract()

    paths = {}
    tables = {
        "hs_completed_lci_supports.csv": records,
        "hs_completed_lci_frame_basis.csv": basis_rows,
        "hs_completed_lci_frame_maps.coo.csv": map_rows,
        "hs_chart_point_dual_coordinate_ledger.csv": chart_rows,
        "hs_overlap_point_dual_coordinate_ledger.csv": overlap_rows,
    }
    for name, rows in tables.items():
        path = output_directory / name
        write_csv(path, rows)
        paths[name] = path

    contract_path = output_directory / "surface_ring_lci_edge_lift_contract.json"
    contract_path.write_text(
        json.dumps(missing, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    paths[contract_path.name] = contract_path

    certificate = {
        "schema": "section12-HS-border-to-completed-local-frames-v1",
        "status": (
            "486-completed-local-lci-and-HS-frames-emitted;"
            "96-Hermite-point-decompositions-bound;"
            "192-overlap-dual-number-squares-certified;"
            "surface-ring-edge-chain-lifts-open"
        ),
        "field": "F_31",
        "counts": {
            "surface_count": 6,
            "K_lifts_per_surface": 81,
            "completed_support_frames": len(records),
            "completed_frame_basis_rows": len(basis_rows),
            "completed_frame_sparse_map_entries": len(map_rows),
            "Hermite_surface_charts": len(chart_rows),
            "Hermite_surface_edges": len(overlap_rows),
        },
        "completed_local_model": {
            "ring": "F_31[[x,y]]",
            "ideal": "I_W=(x,y^2)",
            "ideal_resolution": "R --(-y^2,x)^t--> R^2 --(x,y^2)--> I_W",
            "HS_pushout": "R --(Delta,y^2,-x)^t--> R plus R^2 --> V1",
            "Delta": "delta_0+delta_1*y with delta_0 nonzero",
            "V1_free_frame": True,
            "A_inclusion": "(-y^2/Delta,x/Delta)^t",
            "I_projection": "(x,y^2)",
            "punctured_transition": "[[1,Delta/(x*y^2)],[0,1]]",
            "determinant": "the local transition has determinant 1; global determinant line remains O_S(2H)",
        },
        "exact_checks": {
            **symbolic_checks,
            "all_chart_Hermite_inverse_and_point_dual_checks": True,
            "all_overlap_jet_squares_commute": True,
            "all_overlap_idempotents_and_nilpotents_preserved": True,
        },
        "global_theta_handoff": {
            "completed_support_frame_bases_available": True,
            "completed_support_frame_maps_available": True,
            "quotient_chart_to_overlap_coordinate_maps_available": True,
            "surface_ring_chain_lifts_available": False,
            "reserved_hs_localized_frame_basis_emitted": False,
            "reserved_hs_localized_frame_maps_emitted": False,
            "reason": missing["single_missing_conversion"],
        },
        "smallest_missing_homological_conversion": missing,
        "claim_boundary": {
            "proved": [
                "the exact lci and oriented HS pushout matrices at all 486 completed support stalks",
                "local freeness of every V1 stalk from the nonzero boundary-unit constant",
                "the chart Hermite bases split functorially into point idempotents and tangent nilpotents on all 192 overlaps",
            ],
            "not_claimed": [
                "surface-ring lifts of the quotient overlap maps",
                "the reserved full localized frame pair on the 1296-chart K-stable atlas",
                "a global non-split theta cocycle, Yoneda numerator, or transgression T",
                "a result on Markman's Question 11.4 or the Hodge conjecture",
            ],
        },
        "provenance": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                manifest_path,
                SEED_DIRECTORY / "six_hs_local_seeds.csv",
                SEED_DIRECTORY / "six_hs_support_K_orbits.csv",
                SEED_DIRECTORY / "six_hs_completed_local_cocycles.json",
            )
        },
    }
    certificate_path = output_directory / "hs_border_to_local_frames_certificate.json"
    certificate_path.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    paths[certificate_path.name] = certificate_path

    manifest_out = {
        "schema": "section12-HS-border-to-completed-local-frames-manifest-v1",
        "artifacts": {
            name: {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
            }
            for name, path in sorted(paths.items())
        },
    }
    manifest_out_path = output_directory / "manifest.json"
    manifest_out_path.write_text(
        json.dumps(manifest_out, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return certificate_path


def verify_artifacts(output_directory: Path = OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    certificate = json.loads(
        (output_directory / "hs_border_to_local_frames_certificate.json").read_text(
            encoding="utf-8"
        )
    )
    manifest = json.loads((output_directory / "manifest.json").read_text(encoding="utf-8"))
    if certificate.get("counts", {}).get("completed_support_frames") != 486:
        raise HSBorderFrameError("the completed support frame count changed")
    if certificate.get("counts", {}).get("Hermite_surface_edges") != 192:
        raise HSBorderFrameError("the overlap certification count changed")
    if certificate.get("global_theta_handoff", {}).get(
        "reserved_hs_localized_frame_basis_emitted"
    ):
        raise HSBorderFrameError("a partial completed-local basis was mislabeled global")
    for artifact in manifest.get("artifacts", {}).values():
        path = ROOT / artifact["path"]
        if sha256(path) != artifact["sha256"]:
            raise HSBorderFrameError(f"artifact hash changed: {path}")
    return certificate


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--output-directory", type=Path, default=OUTPUT_DIRECTORY)
    args = parser.parse_args(argv)
    if args.verify:
        certificate = verify_artifacts(args.output_directory)
        print(json.dumps(certificate["counts"], sort_keys=True))
        return 0
    path = write_artifacts(args.output_directory)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
