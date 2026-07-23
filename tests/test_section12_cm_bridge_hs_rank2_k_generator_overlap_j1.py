import csv
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_overlap_j1 as comparison,
)


class HSRank2KGeneratorOverlapJ1Test(unittest.TestCase):
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

    def test_all_chart_pilot_is_portability_not_the_K_target(self):
        assessment = self.certificate["all_chart_pilot_assessment"]
        self.assertFalse(assessment["pilot_supplies_K_target"])
        self.assertEqual(assessment["pilot_ordinary_chart_length"], 112)
        self.assertEqual(assessment["required_C0010_target_length"], 98)

    def test_target_is_independently_normalized_on_the_translated_jets(self):
        normalization = self.certificate["independent_target_normalization"]
        self.assertFalse(normalization["target_defined_as_phi_of_source_matrix"])
        self.assertTrue(
            normalization["translated_tangent_jets_equal_canonical_target_jets"]
        )
        self.assertTrue(
            normalization[
                "independent_target_basis_and_multiplication_recover_canonical_C0010"
            ]
        )
        counts = self.certificate["counts"]
        self.assertEqual(counts["source_chart_length"], 98)
        self.assertEqual(counts["independently_rebuilt_target_length"], 98)
        self.assertEqual(counts["K_overlap_support_points"], 35)
        self.assertEqual(counts["K_overlap_length"], 70)

    def test_both_J0_maps_are_surjective_and_have_the_expected_shape(self):
        artifact = self.output / self.certificate["files"]["overlap_Hermite"]
        with np.load(artifact) as data:
            source = np.asarray(data["source_restriction"], dtype=np.int64)
            target = np.asarray(data["target_restriction"], dtype=np.int64)
        self.assertEqual(source.shape, (70, 98))
        self.assertEqual(target.shape, (70, 98))
        self.assertEqual(comparison.hermite.rank_mod(source), 70)
        self.assertEqual(comparison.hermite.rank_mod(target), 70)

    def test_nonidentity_rational_frames_and_blocks_are_present(self):
        rows = self.read_csv("k_generator_affine_forms.csv")
        self.assertEqual(len(rows), 4)
        self.assertEqual(json.loads(rows[0]["denominator_L"]), [1, 0, 1])
        self.assertEqual(json.loads(rows[3]["denominator_L"]), [1, 20, 14])
        program = self.read_csv("k_generator_j1_factor_program.csv")
        target = [row for row in program if row["endpoint"] == "target"]
        self.assertEqual(len(target), 11)
        self.assertTrue(any(row["domain_generator"] != row["codomain_generator"] for row in target))
        self.assertGreater(self.certificate["counts"]["nonidentity_target_blocks"], 0)

    def test_all_sixteen_graph_chain_identities_are_exact(self):
        rows = self.read_csv("k_generator_j1_verification.csv")
        self.assertEqual(len(rows), 16)
        self.assertTrue(all(int(row["residual_terms"]) == 0 for row in rows))
        self.assertTrue(all(row["residual_is_zero"] == "True" for row in rows))
        checks = self.certificate["exact_checks"]
        self.assertTrue(checks["all_chain_residuals_zero"])
        self.assertFalse(checks["tautological_phi_A_identity_used"])


if __name__ == "__main__":
    unittest.main()
