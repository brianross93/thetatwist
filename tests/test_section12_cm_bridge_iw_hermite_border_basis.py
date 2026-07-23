import tempfile
import unittest
from pathlib import Path

import numpy as np

from experiments.section12_instantiation import (
    cm_bridge_iw_hermite_border_basis as border,
)
from experiments.section12_instantiation import (
    cm_bridge_iw_presentation_solver as iw_solver,
)


class IWHermiteBorderBasisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context = iw_solver.load_problem_context()
        border.prepare_border_context(cls.context)
        cls.z0c0000 = border.chart_border_basis(cls.context, 0, "C0000")

    def test_first_real_chart_reaches_the_exact_double_length(self):
        summary = self.z0c0000["summary"]
        self.assertEqual(summary["support_points"], 59)
        self.assertEqual(summary["quotient_length"], 118)
        self.assertEqual(summary["selected_standard_monomials"], 118)
        self.assertEqual(summary["highest_selected_degree"], 3)
        self.assertEqual(summary["translated_normalized_generators"], 162)
        self.assertEqual(summary["translated_generator_hermite_rank"], 0)

    def test_standard_monomials_are_divisor_closed(self):
        basis = set(self.z0c0000["basis"])
        self.assertIn((0,) * 8, basis)
        for exponent in basis:
            for variable, power in enumerate(exponent):
                if power:
                    divisor = list(exponent)
                    divisor[variable] -= 1
                    self.assertIn(tuple(divisor), basis)

    def test_multiplication_commutes_and_point_coordinates_are_exact(self):
        multiplication = self.z0c0000["multiplication"]
        for left in range(8):
            for right in range(left + 1, 8):
                self.assertFalse(
                    np.any(
                        (
                            multiplication[left] @ multiplication[right]
                            - multiplication[right] @ multiplication[left]
                        )
                        % 31
                    )
                )
        hermite = self.z0c0000["hermite_matrix"]
        dual = np.empty_like(self.z0c0000["hermite_inverse"])
        dual[:, 0::2] = self.z0c0000["point_idempotents"]
        dual[:, 1::2] = self.z0c0000["point_nilpotents"]
        self.assertTrue(np.array_equal(hermite @ dual % 31, np.eye(118, dtype=int)))

    def test_binary_artifact_round_trips(self):
        with tempfile.TemporaryDirectory() as temporary:
            path, summary = border.write_chart_artifact(
                self.z0c0000, Path(temporary)
            )
            self.assertTrue(path.is_file())
            self.assertTrue(summary.is_file())
            with np.load(path) as data:
                self.assertEqual(data["standard_monomials"].shape, (118, 8))
                self.assertEqual(data["multiplication_matrices"].shape, (8, 118, 118))
                self.assertEqual(data["point_idempotent_coefficients"].shape, (118, 59))

    def test_representative_overlap_restrictions_are_surjective_algebra_maps(self):
        edge = self.context["_border_edges"][0]
        target = border.chart_border_basis(
            self.context, 0, edge["target_chart"]
        )
        result = border.overlap_restriction(
            self.context,
            0,
            edge,
            {
                (0, edge["source_chart"]): self.z0c0000,
                (0, edge["target_chart"]): target,
            },
        )
        self.assertEqual(result["summary"]["overlap_quotient_length"], 118)
        self.assertEqual(result["summary"]["source_restriction_rank"], 118)
        self.assertEqual(result["summary"]["target_restriction_rank"], 118)
        self.assertEqual(result["summary"]["source_kernel_dimension"], 0)
        self.assertEqual(result["summary"]["target_kernel_dimension"], 14)


if __name__ == "__main__":
    unittest.main()
