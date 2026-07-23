import unittest
from pathlib import Path

import numpy as np

from experiments.section12_instantiation import (
    cm_bridge_iw_surface_ring_presentation as producer,
)


class IWSurfaceRingPresentationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.output = Path(
            "results/section12_instantiation/"
            "cm_bridge_iw_surface_ring_presentation"
        )
        cls.certificate = producer.verify_artifacts(cls.output)

    def test_exact_presentation_contract(self):
        self.assertEqual(
            self.certificate["exact_presentation"],
            "S^1887 -> S^170 -> I_W -> 0",
        )
        counts = self.certificate["counts"]
        self.assertEqual(counts["original_border_generators"], 297)
        self.assertEqual(counts["corner_border_generators"], 170)
        self.assertEqual(counts["Tietze_unit_cancellations"], 127)
        self.assertEqual(counts["presentation_relation_columns"], 1887)

    def test_d_squared_ledger_is_explicit(self):
        checks = self.certificate["exact_checks"]
        self.assertTrue(checks["2744_Koszul_composites_zero_over_polynomial_ring"])
        self.assertTrue(
            checks["six_divided_difference_composites_equal_surface_equations"]
        )
        self.assertTrue(checks["therefore_d_squared_zero_over_surface_ring"])
        self.assertTrue(
            checks[
                "first_presentation_exact_by_regular_sequence_change_of_rings_and_Tietze_contraction"
            ]
        )

    def test_binary_matrix_and_cancellation_ledger(self):
        path = self.output / "Z_2_C0010_surface_ring_presentation.npz"
        with np.load(path) as artifact:
            self.assertEqual(list(map(int, artifact["matrix_shape"])), [170, 1887])
            self.assertEqual(len(artifact["original_corner_rows"]), 170)
            self.assertEqual(len(artifact["cancelled_noncorner_rows"]), 127)
            self.assertEqual(len(artifact["original_remaining_columns"]), 1887)
            self.assertEqual(len(artifact["removed_zero_original_columns"]), 554)
            self.assertEqual(
                artifact["removed_proportional_original_columns"].shape,
                (182, 3),
            )
            self.assertEqual(int(artifact["syzygy_term_row"].max()), 169)
            self.assertEqual(int(artifact["syzygy_term_column"].max()), 1886)
            self.assertEqual(len(artifact["boundary_term_equation"]), 6)
            self.assertEqual(
                sorted(map(int, artifact["boundary_term_equation"])),
                list(range(6)),
            )

    def test_bound_upstream_hashes_exist(self):
        for record in self.certificate["upstream"].values():
            path = producer.ROOT / record["path"]
            self.assertTrue(path.exists())
            self.assertEqual(producer.sha256(path), record["sha256"])


if __name__ == "__main__":
    unittest.main()
