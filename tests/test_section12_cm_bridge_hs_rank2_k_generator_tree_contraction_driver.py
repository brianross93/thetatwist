import tempfile
import json
import unittest
from pathlib import Path

from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_overlap_j1 as overlap_j1,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_shamash as generator_shamash,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_tree_contraction as tree_contraction,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_tree_contraction_driver as driver,
)


class HSRank2KGeneratorTreeContractionDriverTest(unittest.TestCase):
    def test_standard_generators_have_distinct_outputs(self):
        self.assertEqual(
            driver.GENERATOR_EXPONENTS,
            {
                "e0": (1, 0, 0, 0),
                "e1": (0, 1, 0, 0),
                "e2": (0, 0, 1, 0),
                "e3": (0, 0, 0, 1),
            },
        )
        self.assertEqual(len(set(driver.OUTPUT_DIRECTORIES.values())), 4)
        self.assertEqual(
            driver.predicted_dimensions("e1"),
            {
                "overlap_factors": 31,
                "overlap_length": 62,
                "source_F0": 62,
                "source_F1": 435,
                "source_F2": 1742,
                "source_F3": 6082,
                "source_L0_rows": 436,
                "target_F0": 98,
                "target_F1": 687,
                "target_F2": 2750,
                "target_F3": 9610,
                "target_L0_columns": 688,
            },
        )

    def test_binding_is_cross_module_and_restores_previous_state(self):
        previous = (
            overlap_j1.GENERATOR_EXPONENTS,
            generator_shamash.GENERATOR_EXPONENTS,
            tree_contraction.GENERATOR_EXPONENTS,
            tree_contraction.OUTPUT_DIRECTORY,
        )
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            with driver.bind_generator("e2", output):
                expected = (0, 0, 1, 0)
                self.assertEqual(overlap_j1.GENERATOR_EXPONENTS, expected)
                self.assertEqual(generator_shamash.GENERATOR_EXPONENTS, expected)
                self.assertEqual(tree_contraction.GENERATOR_EXPONENTS, expected)
                self.assertEqual(tree_contraction.OUTPUT_DIRECTORY, output)
        self.assertEqual(
            (
                overlap_j1.GENERATOR_EXPONENTS,
                generator_shamash.GENERATOR_EXPONENTS,
                tree_contraction.GENERATOR_EXPONENTS,
                tree_contraction.OUTPUT_DIRECTORY,
            ),
            previous,
        )

    def test_dry_plan_uses_the_completed_e0_footprint_as_an_estimate(self):
        plan = driver.run_plan("e3")
        self.assertEqual(plan["exponents"], [0, 0, 0, 1])
        self.assertEqual(plan["runs_requested"], 1)
        self.assertEqual(
            plan["estimated_output_bytes"],
            plan["reference_e0_artifact_bytes"],
        )
        self.assertIn("no K coherence", plan["certificate_scope"])
        self.assertTrue(plan["execution_ready"])
        self.assertIsNone(plan["execution_blocker"])
        self.assertFalse(plan["optional_preflight"]["required"])

    def test_unpreflighted_generator_is_still_execution_ready(self):
        with tempfile.TemporaryDirectory() as temporary:
            plan = driver.run_plan("e1", temporary)
            self.assertTrue(plan["execution_ready"])
            self.assertIsNone(plan["execution_blocker"])
            self.assertFalse(
                plan["optional_preflight"]["current_source_hash_bound"]
            )

    def test_source_hash_bound_preflight_is_advisory_and_output_local(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output_with_preflight = root / "with_preflight"
            output_without_preflight = root / "without_preflight"
            output_with_preflight.mkdir()
            output_without_preflight.mkdir()
            dimensions = driver.predicted_dimensions("e2")
            certificate = {
                "schema": "section12-HS-rank2-K-generator-preflight-v1",
                "status": "PASS",
                "generator": "e2",
                "K_generator_exponents": [0, 0, 1, 0],
                "overlap_factors": dimensions["overlap_factors"],
                "overlap_length": dimensions["overlap_length"],
                "source_sha256": driver.preflight_source_hashes(),
            }
            driver.preflight_path("e2", output_with_preflight).write_text(
                json.dumps(certificate), encoding="utf-8"
            )
            present = driver.run_plan("e2", output_with_preflight)
            absent = driver.run_plan("e2", output_without_preflight)
            self.assertTrue(present["execution_ready"])
            self.assertTrue(absent["execution_ready"])
            self.assertTrue(
                present["optional_preflight"]["current_source_hash_bound"]
            )
            self.assertFalse(
                absent["optional_preflight"]["current_source_hash_bound"]
            )


if __name__ == "__main__":
    unittest.main()
