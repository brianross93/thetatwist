import json
import unittest
from pathlib import Path

import numpy as np

from experiments.section12_instantiation import (
    cm_bridge_hs_extension_from_surface_presentation as producer,
)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / (
    "results/section12_instantiation/"
    "cm_bridge_hs_extension_from_surface_presentation"
)


class Section12HSExtensionFromSurfacePresentationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads(
            (RESULTS / "hs_extension_certificate.json").read_text(
                encoding="utf-8"
            )
        )
        cls.artifact_path = RESULTS / "Z_2_C0010_hs_extension.npz"

    def test_tree_shamash_shapes_and_factorization(self):
        tree = self.certificate["tree_shamash_complex"]
        self.assertEqual(tree["ranks_F0_F1_F2_F3"], [1, 687, 2750, 9610])
        self.assertEqual(tree["D2_shape"], [687, 2750])
        factors = tree["D3_factor_shapes"]
        self.assertEqual(factors["K3_to_K2"], [2744, 5488])
        self.assertEqual(factors["six_H1_raw_blocks"], [6, 2744, 784])
        self.assertEqual(factors["tree_i1"], [784, 687])
        self.assertEqual(factors["six_H1_tree_factored"], [6, 2744, 687])
        self.assertEqual(factors["six_scalar_d1_blocks"], [6, 1, 687])

    def test_canonical_class_has_the_correct_degree_and_exact_basis_factor(self):
        cochain = self.certificate["canonical_cochain"]
        self.assertEqual(cochain["cohomological_degree"], 2)
        self.assertEqual(
            cochain["identified_HS_class"],
            "Ext^2_S(A_W,S) = Ext^1_S(I_W,S)",
        )
        self.assertEqual(cochain["lambda_jet_shape"], [1, 2744])
        self.assertEqual(cochain["basis_factor_shape"], [2744, 2744])
        self.assertIn("I_28 tensor", cochain["basis_factor_storage"])
        self.assertTrue(
            self.certificate["module"]["top_row_is_exact_factored_not_absent"]
        )

    def test_all_exact_boundary_ledgers_closed(self):
        checks = self.certificate["exact_checks"]
        for key in (
            "Hermite_and_inverse_multiply_to_I98_both_orders",
            "Hermite_conjugates_all_49_times_8_multiplication_blocks",
            "six_raw_H1_blocks_compose_to_f_j_I784_over_P",
            "therefore_D2_D3_zero_over_S",
            "lambda_d3_equals_sum_minus1_power_j_f_j_mu_j_over_P_supportwise",
            "therefore_lambda_d3_is_zero_over_S",
            "lambda_kills_all_six_repeated_wedge_H1_blocks_over_P_supportwise",
            "basis_change_transports_both_identities_to_standard_basis",
            "therefore_lambda_is_a_closed_Ext2_cochain",
        ):
            self.assertTrue(checks[key], key)
        ledger = self.certificate["verification_ledger"]
        self.assertEqual(ledger["full_H0_nullhomotopy_columns_checked"], 588)
        self.assertEqual(ledger["raw_H1_columns_checked"], 4704)
        self.assertEqual(ledger["hermite_multiplication_conjugacy_checks"], 392)
        self.assertEqual(
            ledger["lambda_supportwise_P_boundary_ledgers_checked"], 49
        )
        self.assertEqual(
            ledger["lambda_supportwise_repeated_wedge_maps_checked"], 294
        )
        self.assertEqual(len(ledger["lambda_P_boundary_hashes"]), 49)

    def test_binary_contains_executable_top_and_explicit_bottom_factors(self):
        with np.load(self.artifact_path) as artifact:
            self.assertEqual(
                list(map(int, artifact["tree_complex_ranks"])),
                [1, 687, 2750, 9610],
            )
            hermite = np.asarray(artifact["hermite_standard_to_jets"], dtype=int)
            inverse = np.asarray(artifact["hermite_jets_to_standard"], dtype=int)
            self.assertEqual(hermite.shape, (98, 98))
            self.assertEqual(inverse.shape, (98, 98))
            self.assertTrue(np.array_equal(hermite @ inverse % 31, np.eye(98, dtype=int)))
            self.assertTrue(
                np.all(artifact["lambda_K2_jet_row"] == 0)
            )
            self.assertLess(int(artifact["lambda_K2_jet_column"].max()), 2744)
            self.assertLess(
                int(artifact["Vprime_relation_bottom_minus_D2_row"].max()), 687
            )
            self.assertLess(
                int(artifact["Vprime_relation_bottom_minus_D2_column"].max()),
                2750,
            )

    def test_scope_does_not_promote_the_chart_to_a_global_object(self):
        boundary = self.certificate["claim_boundary"]
        self.assertIn("one exact chart-local HS cochain", boundary["proved"])
        for phrase in (
            "global O(3H)/O(-H) overlap descent",
            "relative O(4H) transition",
            "projective idempotent",
            "297-carrier degree3 lift B",
            "Hodge conjecture",
        ):
            self.assertIn(phrase, boundary["not_claimed"])


if __name__ == "__main__":
    unittest.main()
