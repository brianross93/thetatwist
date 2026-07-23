import csv
import json
import tempfile
import unittest
from pathlib import Path

from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_chain_comparison as audit,
)


class HSRank2KChainComparisonTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.output = Path(cls.temporary.name)
        audit.build(cls.output)
        cls.certificate = json.loads(
            (cls.output / "k_chain_comparison_certificate.json").read_text()
        )
        cls.contract = json.loads(
            (cls.output / "k_chain_comparison_input_contract.json").read_text()
        )

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def read_csv(self, name):
        with (self.output / name).open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    def test_full_orbit_and_cayley_graph_are_indexed(self):
        counts = self.certificate["counts"]
        self.assertEqual(counts["K_elements"], 81)
        self.assertEqual(counts["K_stable_chart_rows_replayed"], 1296)
        self.assertEqual(counts["factor_affine_map_rows"], 324)
        self.assertEqual(counts["typed_Cayley_edges"], 324)
        self.assertEqual(len(self.read_csv("z2_k_orbit_ring_maps.csv")), 324)
        self.assertEqual(
            len(self.read_csv("z2_k_cayley_chain_comparison_jobs.csv")), 324
        )

    def test_surface_equations_and_line_factor_are_group_coherent(self):
        counts = self.certificate["counts"]
        self.assertEqual(counts["group_product_checks"], 81**2)
        self.assertEqual(counts["projective_factor_matrix_checks"], 4 * 81**2)
        self.assertEqual(counts["normalized_H_cocycle_checks"], 81**2)
        self.assertEqual(counts["normalized_O4H_cocycle_checks"], 81**2)
        self.assertEqual(
            counts["normalized_Z2_theta_character_checks"], 2 * 81**2
        )
        rows = self.read_csv("z2_surface_equation_transport.csv")
        self.assertEqual(len(rows), 81)
        self.assertTrue(all(row["surface_preserved"] == "True" for row in rows))
        self.assertTrue(self.certificate["exact_line_factor"]["group_coherent"])

    def test_identity_element_is_the_base_affine_map(self):
        rows = [
            row
            for row in self.read_csv("z2_k_orbit_ring_maps.csv")
            if int(row["K_element_index"]) == 0
        ]
        self.assertEqual(len(rows), 4)
        for row in rows:
            self.assertEqual(
                json.loads(row["target_frame_pullback_in_source_1_u_v"]),
                [1, 0, 0],
            )
            self.assertEqual(
                json.loads(row["target_u_numerator_in_source_1_u_v"]),
                [0, 1, 0],
            )
            self.assertEqual(
                json.loads(row["target_v_numerator_in_source_1_u_v"]),
                [0, 0, 1],
            )

    def test_contract_does_not_promote_base_change_to_linearization(self):
        verdict = self.contract["verdict"]
        self.assertFalse(verdict["one_nonidentity_generator_L0_L1_currently_certified"])
        self.assertTrue(verdict["next_computation_is_well_posed"])
        boundary = self.certificate["whole_ring_boundary"]
        self.assertTrue(boundary["source_A_D3_lambda_bound"])
        self.assertFalse(boundary["nonidentity_J1_J2_bound"])
        self.assertFalse(boundary["nonidentity_L0_L1_identity_bound"])
        self.assertFalse(boundary["supportwise_evaluation_used_as_substitute"])
        self.assertFalse(boundary["tautological_phi_A_identity_used_as_linearization"])
        coherence = next(
            row
            for row in self.contract["smallest_unemitted_exact_data"]
            if row["name"] == "K_group_coherence"
        )
        self.assertEqual(coherence["generator_square_cells"], 486)
        self.assertEqual(coherence["distinct_generator_order_three_cycles"], 108)
        self.assertEqual(coherence["all_typed_pair_composites"], 81**2)

    def test_four_generator_readiness_rows_name_the_missing_matrices(self):
        rows = self.read_csv("generator_comparison_readiness.csv")
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(row["rational_ring_map"] == "exact_and_emitted" for row in rows))
        self.assertTrue(all(row["O4H_factor"] == "exact_normalized_circuit_emitted" for row in rows))
        self.assertTrue(all(row["J1_687x687"] == "not_emitted" for row in rows))
        self.assertTrue(all(row["J2_2750x2750"] == "not_emitted" for row in rows))


if __name__ == "__main__":
    unittest.main()
