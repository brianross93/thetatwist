import unittest

from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_metadata_repair as repair,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_tree_contraction_driver as driver,
)


class HSRank2KGeneratorMetadataRepairTest(unittest.TestCase):
    def test_e1_sidecar_reconciles_hash_encodings_without_rewriting_maps(self):
        output = driver.OUTPUT_DIRECTORIES["e1"]
        if not (output / "tree_contracted_support_maps.npz").is_file():
            self.skipTest("the full e1 artifact is not present")
        certificate = repair.build_certificate(output, "e1")
        self.assertEqual(certificate["status"], "PASS")
        self.assertFalse(
            certificate["consistency"]["numeric_artifacts_rewritten"]
        )
        self.assertTrue(
            certificate["consistency"][
                "all_original_uint8_hashes_match"
            ]
        )
        self.assertTrue(
            certificate["consistency"][
                "all_available_canonical_preflight_hashes_match"
            ]
        )
        self.assertTrue(
            certificate["consistency"][
                "target_hash_matches_original_counts"
            ]
        )
        self.assertEqual(
            {
                tuple(ledger["shape"])
                for ledger in certificate["pair_ledgers"].values()
            },
            {(31, 2)},
        )


if __name__ == "__main__":
    unittest.main()
