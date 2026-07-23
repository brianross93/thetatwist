from __future__ import annotations

import ast
import hashlib
import inspect
import json
import unittest

import numpy as np

from experiments.section12_instantiation import (
    cm_bridge_hs_all_chart_presentations as all_charts,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_extension_from_surface_presentation as representative,
)


class HSAllChartPresentationsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.output = all_charts.OUTPUT_DIRECTORY
        cls.certificate_path = cls.output / "Z_2_C0000_certificate.json"
        cls.artifact_path = cls.output / "Z_2_C0000_hs_extension.npz"
        cls.inventory_path = cls.output / "all_chart_input_inventory.json"
        cls.certificate = json.loads(
            cls.certificate_path.read_text(encoding="utf-8")
        )
        cls.inventory = json.loads(
            cls.inventory_path.read_text(encoding="utf-8")
        )

    def test_all_96_inputs_are_bound_by_one_dimension_law(self):
        self.assertEqual(self.inventory["chart_count"], 96)
        self.assertEqual(self.inventory["surfaces"], 6)
        self.assertEqual(self.inventory["charts_per_surface"], 16)
        rows = self.inventory["rows"]
        self.assertEqual(
            {(row["surface"], row["chart"]) for row in rows},
            {
                (surface, chart)
                for surface in range(6)
                for chart in all_charts.chart_names()
            },
        )
        for row in rows:
            dimension = row["finite_algebra_dimension"]
            self.assertEqual(dimension, 2 * row["support_count"])
            self.assertEqual(row["tree_F1_rank"], 7 * dimension + 1)
            self.assertEqual(row["F2_rank"], 28 * dimension + 6)
            self.assertEqual(
                row["F3_factored_rank"],
                56 * dimension + 6 * (7 * dimension + 1),
            )

    def test_second_chart_has_the_derived_exact_dimensions(self):
        dimensions = self.certificate["dimensions"]
        self.assertEqual(self.certificate["surface"], "Z_2")
        self.assertEqual(self.certificate["chart"], "C0000")
        self.assertEqual(dimensions["finite_algebra"], 112)
        self.assertEqual(dimensions["support_count"], 56)
        self.assertEqual(dimensions["raw_K1"], 896)
        self.assertEqual(dimensions["tree_F1"], 785)
        self.assertEqual(dimensions["K2"], 3136)
        self.assertEqual(dimensions["F2"], 3142)
        self.assertEqual(dimensions["K3"], 6272)
        self.assertEqual(dimensions["F3_factored"], 10982)
        self.assertEqual(dimensions["module_target"], 786)
        self.assertEqual(dimensions["border_carrier"], 313)
        self.assertEqual(dimensions["merge_kernel"], 472)

    def test_chain_cocycle_and_orientation_checks_closed(self):
        checks = self.certificate["exact_checks"]
        for key in (
            "Hermite_basis_invertible",
            "six_surface_equations_annihilate_A_W",
            "tree_inclusion_and_projection_split",
            "D2_D3_zero_over_surface_ring",
            "lambda_D3_zero_over_surface_ring",
            "A_D3_zero_over_surface_ring",
            "lambda_is_closed_Ext2_cochain",
            "tree_model_merges_to_independent_border_model",
            "C0010_curve_reduction_helpers_not_used_by_this_build",
        ):
            self.assertTrue(checks[key], key)
        binding = self.certificate["orientation_binding"]
        self.assertEqual(binding["rows"], 56)
        self.assertTrue(binding["all_support_indices_found_in_surface_ledger"])
        self.assertTrue(binding["all_delta0_nonzero"])

    def test_parameterized_build_does_not_call_legacy_curve_reduction(self):
        for function in (
            all_charts.load_chart_inputs,
            all_charts.build_chart,
            all_charts.write_chart,
        ):
            self.assertNotIn("_curve_reduce", inspect.getsource(function))

        source = inspect.getsource(representative)
        tree = ast.parse(source)
        callers = {}

        class CurveCallVisitor(ast.NodeVisitor):
            def __init__(self):
                self.function = None

            def visit_FunctionDef(self, node):
                previous = self.function
                self.function = node.name
                self.generic_visit(node)
                self.function = previous

            def visit_Call(self, node):
                name = node.func.id if isinstance(node.func, ast.Name) else None
                if name and name.startswith("_curve_"):
                    callers.setdefault(name, set()).add(self.function)
                self.generic_visit(node)

        CurveCallVisitor().visit(tree)
        self.assertEqual(
            callers.get("_curve_reduce_monomial"), {"_curve_reduce_polynomial"}
        )
        self.assertEqual(
            callers.get("_curve_reduce_polynomial"), {"_curve_multiply"}
        )
        self.assertNotIn("_curve_multiply", callers)

    def test_cocycle_metadata_uses_the_selected_chart_support_count(self):
        cocycle = self.certificate["change_of_rings_cocycle"]
        self.assertEqual(cocycle["supportwise_dual_number_factors"], 56)
        self.assertIn("56 dual-number factors", cocycle["mu_factorization"])
        self.assertNotIn("49 dual-number factors", cocycle["mu_factorization"])

    def test_claim_boundary_does_not_promote_chart_presentations_to_k_targets(self):
        boundary = self.certificate["claim_boundary"]
        self.assertIn("chart presentation", boundary["proved"])
        self.assertIn("independently normalized K-target", boundary["not_claimed"])
        self.assertIn("K-chain comparison contract", boundary["not_claimed"])

    def test_factored_artifact_shapes_match_certificate(self):
        with np.load(self.artifact_path) as artifact:
            self.assertEqual(
                artifact["tree_complex_ranks"].tolist(),
                [1, 785, 3142, 10982],
            )
            self.assertEqual(
                artifact["module_presentation_shape"].tolist(), [786, 3142]
            )
            self.assertEqual(artifact["support_point_indices"].shape, (56, 4))
            self.assertEqual(artifact["hermite_standard_to_jets"].shape, (112, 112))
            self.assertEqual(artifact["merge_q"].shape, (313, 785))
            self.assertEqual(artifact["merge_kernel"].shape, (785, 472))
            self.assertEqual(
                int(artifact["D2_tree_column"].max()) + 1, 3142
            )
            self.assertEqual(
                int(artifact["D3_K3_to_K2_column"].max()) + 1, 6272
            )

    def test_manifest_and_upstream_custody(self):
        manifest = json.loads(
            (self.output / "manifest.json").read_text(encoding="utf-8")
        )
        for name, record in manifest.items():
            path = self.output / name
            self.assertEqual(path.stat().st_size, record["bytes"])
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(), record["sha256"]
            )
        for record in self.certificate["upstream"].values():
            path = all_charts.ROOT / record["path"]
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(), record["sha256"]
            )


if __name__ == "__main__":
    unittest.main()
