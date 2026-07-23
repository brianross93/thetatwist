import csv
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_shamash as comparison,
)


class HSRank2KGeneratorShamashTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.output = Path(cls.temporary.name)
        comparison.write_artifacts(cls.output)
        cls.certificate = comparison.verify_artifacts(cls.output)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def read_csv(self, name):
        with (self.output / name).open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    def test_one_common_law_closes_J2_and_the_D3_inputs(self):
        counts = self.certificate["counts"]
        checks = self.certificate["exact_checks"]
        self.assertEqual(counts["supports"], 35)
        self.assertEqual(counts["equations"], 6)
        self.assertEqual(counts["J2_component_identities"], 245)
        self.assertEqual(counts["D3_compatibility_inputs"], 455)
        self.assertTrue(checks["one_common_diagonal_B_law_for_all_six_equations"])
        self.assertTrue(checks["Gamma_not_six_fitted_repairs"])

    def test_adapted_orientation_and_counit_weight_are_exact(self):
        normalization = self.certificate["normalization"]
        self.assertEqual(normalization["adapted_source_orientation_scalar"], 6)
        self.assertEqual(normalization["B_determinant_scalar"], 6)
        self.assertEqual(normalization["residue_transport_scalar"], 27)
        self.assertEqual(normalization["raw_weight_scalar"], 20)
        self.assertEqual(normalization["normalized_weight_scalar"], 7)
        self.assertEqual(
            normalization["target_lambda_raw_over_normalized_scalar"], 25
        )
        self.assertEqual(normalization["top_J8_residue_transport_checks"], 35)

    def test_all_single_counit_rows_solve_all_34_columns(self):
        counts = self.certificate["counts"]
        checks = self.certificate["exact_checks"]
        self.assertEqual(counts["counit_systems_solved"], 35)
        self.assertEqual(counts["counit_systems_obstructed"], 0)
        self.assertTrue(checks["all_34_counit_columns_zero_in_local_quotient"])
        self.assertEqual(
            checks["K2_counit_class_homotopy_invariance_checks"], 35
        )
        self.assertEqual(self.certificate["counit_system"]["rank_per_support"], [12])
        self.assertEqual(
            self.certificate["counit_system"]["nullity_per_support"], [4]
        )

    def test_persisted_systems_and_maps_are_nonempty_and_exact(self):
        with np.load(
            self.output / self.certificate["files"]["counit_systems"]
        ) as data:
            systems = np.asarray(data["systems"], dtype=np.int64) % comparison.P
            rhs = np.asarray(data["right_hand_sides"], dtype=np.int64) % comparison.P
            solutions = (
                np.asarray(data["canonical_solutions"], dtype=np.int64)
                % comparison.P
            )
        self.assertEqual(systems.shape, (35, 68, 16))
        self.assertFalse(
            np.any(
                (
                    np.einsum("sij,sj->si", systems, solutions) - rhs
                )
                % comparison.P
            )
        )
        with np.load(self.output / self.certificate["files"]["maps"]) as data:
            self.assertGreater(data["graph_J2_records"].shape[0], 0)
            self.assertGreater(data["graph_J3_records"].shape[0], 0)
            self.assertGreater(data["Gamma_records"].shape[0], 0)
            self.assertGreater(data["Omega_records"].shape[0], 0)
            self.assertEqual(data["counit_a_records"].shape[0], 0)
        self.assertTrue(
            self.certificate["counit_system"][
                "canonical_solution_a_is_zero_on_all_supports"
            ]
        )

    def test_support_ledger_has_no_residual_or_fitted_exception(self):
        rows = self.read_csv(self.certificate["files"]["support_summary"])
        self.assertEqual(len(rows), 35)
        self.assertTrue(all(row["counit_status"] == "SOLVED" for row in rows))
        self.assertTrue(all(row["all_34_counit_columns_zero"] == "True" for row in rows))
        self.assertTrue(
            all(
                row["K2_counit_class_homotopy_invariant"] == "True"
                for row in rows
            )
        )
        self.assertTrue(
            all(json.loads(row["nonzero_K2_residual_coordinates"]) == [] for row in rows)
        )

    def test_claim_boundary_keeps_tree_contraction_open(self):
        boundary = self.certificate["claim_boundary"]
        self.assertIn("1966-by-2750", boundary["not_claimed"])
        self.assertIn("tree", boundary["next_exact_construction"])


if __name__ == "__main__":
    unittest.main()
