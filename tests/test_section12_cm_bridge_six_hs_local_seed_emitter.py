import csv
import json
import tempfile
import unittest
from pathlib import Path

from experiments.section12_instantiation import cm_bridge_six_hs_local_seed_emitter as producer


class SixHartshorneSerreLocalSeedEmitterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.output = Path(cls.temporary.name)
        cls.certificate = producer.write_artifacts(cls.output)
        with (cls.output / "six_hs_local_seeds.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            cls.seeds = list(csv.DictReader(handle))
        with (cls.output / "surface_pair_transport_audit.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            cls.transport = list(csv.DictReader(handle))
        with (cls.output / "six_hs_support_K_orbits.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            cls.support_orbits = list(csv.DictReader(handle))
        cls.cocycles = json.loads(
            (cls.output / "six_hs_completed_local_cocycles.json").read_text(
                encoding="utf-8"
            )
        )

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_six_distinct_boundary_avoiding_seeds_are_emitted(self):
        self.assertEqual(len(self.seeds), 6)
        self.assertEqual(
            len({row["quotient_orbit_representative"] for row in self.seeds}), 6
        )
        self.assertTrue(all(int(row["quotient_orbit_size"]) == 81 for row in self.seeds))
        self.assertTrue(all(int(row["delta0_mod_31"]) for row in self.seeds))
        self.assertEqual(self.certificate["total_double_point_length"], 12)
        self.assertEqual(len(self.support_orbits), 6 * 81)

    def test_original_Z2_seed_and_residue_oracle_are_preserved(self):
        seed = self.seeds[2]
        self.assertEqual(json.loads(seed["source_point_indices"]), [0, 0, 9, 15])
        self.assertEqual(json.loads(seed["tangent_fixed_frame_coefficients"]), [1, 1, 3, 0])
        self.assertEqual(
            [int(seed["delta0_mod_31"]), int(seed["delta1_mod_31"])],
            [28, 17],
        )
        self.assertEqual(
            [
                int(seed["residue_f0_coefficient_mod_31"]),
                int(seed["residue_f1_coefficient_mod_31"]),
            ],
            [17, 28],
        )

    def test_independent_surface_pairs_are_not_falsely_cyclicized(self):
        self.assertEqual(len(self.transport), 6)
        for index, row in enumerate(self.transport):
            self.assertEqual(
                json.loads(row["translation_equivalent_surface_indices"]), [index]
            )
            self.assertEqual(
                row["transport_from_Z2_available"], "True" if index == 2 else "False"
            )

    def test_completed_local_cocycles_are_explicit_but_not_global_T(self):
        self.assertEqual(len(self.cocycles), 6)
        for cocycle in self.cocycles:
            self.assertEqual(cocycle["double_point_ideal"], "(x_s,y_s^2)")
            self.assertEqual(
                cocycle["transition_Dx_to_Dy"][0][1],
                "delta_s/(x_s*y_s^2)",
            )
            self.assertIn("not a global theta-chart cocycle", cocycle["scope"])
        self.assertFalse(self.certificate["global_output"]["T_computed"])
        self.assertFalse(
            self.certificate["global_output"]["global_Yoneda_numerators_computed"]
        )

    def test_main_lane_certificate_is_the_determinant_two_H_pair(self):
        pair = self.certificate["main_lane_pair"]
        self.assertEqual(pair["det_V_HS"], "O_S(2H)")
        self.assertEqual(pair["split_bundle"], "V_split=O_S(3H) plus O_S(-H)")


if __name__ == "__main__":
    unittest.main()
