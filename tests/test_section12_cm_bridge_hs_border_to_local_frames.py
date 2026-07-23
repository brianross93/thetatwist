import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_border_to_local_frames"
)


class HSBorderToLocalFramesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from experiments.section12_instantiation import (
            cm_bridge_hs_border_to_local_frames as subject,
        )

        cls.subject = subject
        cls.certificate = subject.verify_artifacts(OUTPUT)

    @staticmethod
    def rows(name):
        with (OUTPUT / name).open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    def test_all_486_completed_support_frames_are_present(self):
        counts = self.certificate["counts"]
        self.assertEqual(counts["surface_count"], 6)
        self.assertEqual(counts["K_lifts_per_surface"], 81)
        self.assertEqual(counts["completed_support_frames"], 486)
        supports = self.rows("hs_completed_lci_supports.csv")
        self.assertEqual(len(supports), 486)
        self.assertEqual(
            len({(row["surface"], row["K_lift_index"]) for row in supports}),
            486,
        )
        for row in supports:
            d0 = int(row["delta0_mod_31"])
            d1 = int(row["delta1_mod_31"])
            e0 = int(row["delta_inverse0_mod_31"])
            e1 = int(row["delta_inverse1_mod_31"])
            self.assertNotEqual(d0, 0)
            self.assertEqual(d0 * e0 % 31, 1)
            self.assertEqual((d0 * e1 + d1 * e0) % 31, 0)

    def test_completed_lci_and_pushout_matrices_are_emitted(self):
        maps = self.rows("hs_completed_lci_frame_maps.coo.csv")
        self.assertEqual(len(maps), 486 * 17)
        by_support = {}
        for row in maps:
            by_support.setdefault(row["support_id"], []).append(row)
        self.assertEqual(len(by_support), 486)
        expected_kinds = {
            "lci_resolution",
            "HS_pushout_presentation",
            "HS_exact_sequence",
            "punctured_splitting",
            "completed_local_transition",
            "completed_local_transition_inverse",
        }
        for rows in by_support.values():
            self.assertEqual({row["map_kind"] for row in rows}, expected_kinds)
            self.assertEqual(len(rows), 17)

    def test_exact_local_frame_identities_are_bound(self):
        checks = self.certificate["exact_checks"]
        self.assertEqual(checks["support_instances_checked"], 486)
        self.assertEqual(
            checks["lci_complex_q_times_d"], "x*(-y^2)+y^2*x=0"
        )
        self.assertTrue(checks["all_chart_Hermite_inverse_and_point_dual_checks"])
        self.assertTrue(checks["all_overlap_jet_squares_commute"])
        self.assertTrue(checks["all_overlap_idempotents_and_nilpotents_preserved"])

    def test_all_chart_and_overlap_coordinate_bridges_are_certified(self):
        charts = self.rows("hs_chart_point_dual_coordinate_ledger.csv")
        overlaps = self.rows("hs_overlap_point_dual_coordinate_ledger.csv")
        self.assertEqual(len(charts), 96)
        self.assertEqual(len(overlaps), 192)
        self.assertTrue(all(row["two_sided_coordinate_check"] == "True" for row in charts))
        self.assertTrue(all(row["source_jet_square_commutes"] == "True" for row in overlaps))
        self.assertTrue(all(row["target_jet_square_commutes"] == "True" for row in overlaps))
        self.assertTrue(all(row["point_idempotents_preserved"] == "True" for row in overlaps))
        self.assertTrue(all(row["tangent_nilpotents_preserved"] == "True" for row in overlaps))

    def test_only_the_surface_ring_chain_lift_remains(self):
        handoff = self.certificate["global_theta_handoff"]
        self.assertTrue(handoff["completed_support_frame_bases_available"])
        self.assertTrue(handoff["completed_support_frame_maps_available"])
        self.assertTrue(handoff["quotient_chart_to_overlap_coordinate_maps_available"])
        self.assertFalse(handoff["surface_ring_chain_lifts_available"])
        self.assertFalse(handoff["reserved_hs_localized_frame_basis_emitted"])
        self.assertFalse(handoff["reserved_hs_localized_frame_maps_emitted"])
        contract = json.loads(
            (OUTPUT / "surface_ring_lci_edge_lift_contract.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            contract["required_identity"]["equations"],
            ["q_s*U_e=q_t", "U_e*d_t=d_s*lambda_e"],
        )
        self.assertTrue(contract["no_new_geometric_parameter"])

    def test_claim_boundary_is_explicit(self):
        boundary = self.certificate["claim_boundary"]
        self.assertIn(
            "surface-ring lifts of the quotient overlap maps",
            boundary["not_claimed"],
        )
        self.assertIn(
            "a result on Markman's Question 11.4 or the Hodge conjecture",
            boundary["not_claimed"],
        )


if __name__ == "__main__":
    unittest.main()
