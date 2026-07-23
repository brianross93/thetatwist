import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from experiments.section12_instantiation import cm_bridge_global_hs_chart_resolution as producer


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/section12_instantiation/cm_bridge_global_hs_chart_resolution"


def read_csv(name):
    with (RESULTS / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class GlobalHSChartResolutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads(
            (RESULTS / "global_hs_chart_resolution_certificate.json").read_text(
                encoding="utf-8"
            )
        )

    def test_common_product_atlas_is_the_exact_four_cube(self):
        atlas = self.certificate["common_product_atlas"]
        self.assertEqual(atlas["charts"], 16)
        self.assertEqual(atlas["adjacent_edges"], 32)
        self.assertEqual(atlas["square_faces"], 24)
        self.assertEqual(atlas["rank_d0_mod_31"], 15)
        self.assertEqual(atlas["rank_d1_mod_31"], 17)
        self.assertEqual(atlas["scalar_cube_H1_dimension"], 0)
        self.assertEqual(len(read_csv("product_theta_atlas.csv")), 16)
        self.assertEqual(len(read_csv("product_theta_atlas_edges.csv")), 32)
        self.assertEqual(len(read_csv("product_theta_atlas_faces.csv")), 24)

    def test_all_six_actual_degree_two_jet_systems_have_the_certified_ranks(self):
        summary = read_csv("degree2_global_jet_system_summary.csv")
        self.assertEqual(len(summary), 6)
        for row in summary:
            self.assertEqual(int(row["unknown_coefficients"]), 1296)
            self.assertEqual(int(row["K_orbit_lifts"]), 81)
            self.assertEqual(int(row["value_plus_tangent_equations"]), 162)
            self.assertEqual(int(row["value_plus_tangent_rank_mod_31"]), 161)
            self.assertEqual(int(row["value_plus_tangent_kernel_dimension"]), 1135)
            self.assertEqual(int(row["full_equation_count"]), 164)
            self.assertEqual(int(row["full_rank_mod_31"]), 163)
            self.assertEqual(int(row["full_solution_dimension"]), 1133)
            self.assertEqual(row["global_generators_solved"], "True")
            self.assertEqual(row["final_Cech_cocycle_solved"], "False")

    def test_actual_generator_coefficients_and_normalizations_are_emitted(self):
        coefficients = read_csv("degree2_global_generator_coefficients.coo.csv")
        normalizations = read_csv("degree2_global_generator_normalization.csv")
        self.assertTrue(coefficients)
        self.assertEqual(
            {(row["surface"], row["generator"]) for row in coefficients},
            {(f"Z_{surface}", generator) for surface in range(6) for generator in ("X", "Q")},
        )
        self.assertEqual(len(normalizations), 12)
        for row in normalizations:
            self.assertEqual(row["vanishes_on_all_81_length_two_lifts"], "True")
            if row["generator"] == "X":
                self.assertEqual(row["representative_transverse_derivative_mod_31"], "1")
                self.assertEqual(row["representative_second_tangent_derivative_mod_31"], "0")
            else:
                self.assertEqual(row["representative_transverse_derivative_mod_31"], "0")
                self.assertEqual(row["representative_second_tangent_derivative_mod_31"], "1")

    def test_supports_are_assigned_to_actual_product_charts(self):
        assignments = read_csv("support_chart_assignment.csv")
        self.assertEqual(len(assignments), 486)
        self.assertTrue(all(len(row["canonical_product_chart"]) == 5 for row in assignments))
        self.assertTrue(
            all(set(row["canonical_product_chart"][1:]) <= {"0", "1"} for row in assignments)
        )

    def test_final_cocycle_is_not_fabricated(self):
        output = self.certificate["global_output"]
        self.assertTrue(output["global_sections_and_jets_computed"])
        self.assertFalse(output["global_hs_theta_chart_cocycle_computed"])
        self.assertFalse(output["global_Yoneda_numerators_computed"])
        self.assertFalse(output["T_computed"])
        contract = json.loads(
            (RESULTS / "global_hs_theta_chart_cocycle_linear_system.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(contract["target_output"], "global_hs_theta_chart_cocycle.coo.csv")
        self.assertEqual(
            contract["common_atlas"], self.certificate["common_product_atlas"]
        )
        self.assertIn(
            "chart-local syzygy bases",
            contract["smallest_unsolved_sheaf_valued_system"]["unknown_coefficient_blocks"][0],
        )

    def test_manifest_and_reproduction(self):
        manifest = json.loads((RESULTS / "manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(
            manifest["generator_sha256"],
            hashlib.sha256(Path(producer.__file__).read_bytes()).hexdigest(),
        )
        for name, record in manifest["files"].items():
            path = RESULTS / name
            self.assertEqual(path.stat().st_size, record["bytes"])
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), record["sha256"])

        with tempfile.TemporaryDirectory() as temporary:
            fresh = Path(temporary)
            producer.write_artifacts(fresh)
            self.assertEqual(
                sorted(path.name for path in fresh.iterdir()),
                sorted(path.name for path in RESULTS.iterdir()),
            )
            for path in fresh.iterdir():
                self.assertEqual(path.read_bytes(), (RESULTS / path.name).read_bytes())


if __name__ == "__main__":
    unittest.main()
