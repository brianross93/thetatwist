import csv
import hashlib
import json
import unittest
from pathlib import Path

from experiments.section12_instantiation import (
    cm_bridge_hs_chart_local_syzygy_audit as producer,
)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / (
    "results/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit"
)


def read_csv(name):
    with (RESULTS / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class HSChartLocalSyzygyAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads(
            (RESULTS / "hs_chart_local_syzygy_audit_certificate.json").read_text(
                encoding="utf-8"
            )
        )

    def test_exact_degree_two_K_action_closes(self):
        action = self.certificate["K_action"]
        self.assertEqual(action["group_order"], 81)
        self.assertEqual(action["normalized_group_law_checks"], 6561)
        self.assertFalse(action["translation_coefficients_fabricated"])
        rows = read_csv("K_degree2_action_elements.csv")
        self.assertEqual(len(rows), 81)
        self.assertEqual(
            {tuple(json.loads(row["K_generator_exponents"])) for row in rows},
            set(__import__("itertools").product(range(3), repeat=4)),
        )

    def test_all_support_symbols_are_locally_surjective(self):
        result = self.certificate["finite_symbol_result"]
        self.assertEqual(result["surfaces"], 6)
        self.assertEqual(result["generators_per_surface"], 162)
        self.assertEqual(result["support_lifts_per_surface"], 81)
        self.assertEqual(result["every_support_block_rank"], 2)
        self.assertEqual(
            result["joint_constant_symbol_ranks"], [158, 106, 157, 106, 107, 158]
        )
        self.assertEqual(result["automorphy_block_checks"], 6 * 81 * 81)
        self.assertEqual(result["direct_quadratic_transport_checks"], 6 * 4 * 2 * 81 * 2)
        summary = read_csv("surface_conormal_symbol_summary.csv")
        self.assertEqual(len(summary), 6)
        self.assertTrue(
            all(int(row["every_2x162_support_block_rank"]) == 2 for row in summary)
        )

    def test_chart_and_edge_nullities_are_exact(self):
        charts = read_csv("chart_conormal_symbol_summary.csv")
        edges = read_csv("edge_conormal_symbol_summary.csv")
        self.assertEqual(len(charts), 96)
        self.assertEqual(len(edges), 192)
        for row in [*charts, *edges]:
            self.assertEqual(
                int(row["constant_symbol_syzygy_dimension"]),
                162 - int(row["symbol_rank_mod_31"]),
            )
            self.assertLessEqual(
                int(row["symbol_rank_mod_31"]), int(row["symbol_rows"])
            )
        self.assertTrue(read_csv("chart_constant_symbol_syzygy_basis.coo.csv"))
        self.assertTrue(read_csv("edge_constant_symbol_syzygy_basis.coo.csv"))
        self.assertTrue(
            read_csv("chart_to_edge_constant_symbol_restrictions.coo.csv")
        )

    def test_the_original_product_atlas_is_not_K_stable(self):
        atlas = self.certificate["atlas_result"]
        self.assertFalse(atlas["K_stable"])
        self.assertEqual(atlas["original_product_charts"], 16)
        self.assertEqual(atlas["original_adjacent_edges"], 32)
        self.assertEqual(atlas["distinct_K_translates_of_product_charts"], 1296)
        rows = read_csv("K_stable_product_chart_refinement.csv")
        self.assertEqual(len(rows), 1296)
        self.assertEqual(
            len({row["pulled_back_factor_linear_forms_Z_X_Y"] for row in rows}),
            1296,
        )

    def test_full_sheaf_boundary_is_not_overclaimed(self):
        boundary = self.certificate["claim_boundary"]
        self.assertIn(
            "a polynomial-module syzygy basis over any affine surface chart ring",
            boundary["not_claimed"],
        )
        contract = json.loads(
            (RESULTS / "full_sheaf_syzygy_next_contract.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertTrue(contract["no_guessed_block_size"])
        self.assertEqual(len(contract["smallest_remaining_exact_operations"]), 4)

    def test_manifest_hashes(self):
        manifest = json.loads((RESULTS / "manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(
            manifest["generator_sha256"],
            hashlib.sha256(Path(producer.__file__).read_bytes()).hexdigest(),
        )
        for name, record in manifest["files"].items():
            path = RESULTS / name
            self.assertEqual(path.stat().st_size, record["bytes"])
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), record["sha256"])


if __name__ == "__main__":
    unittest.main()
