from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from thetatwist.core import (
    certificate,
    corrected_character,
    hodge_contraction,
    natural_repair_sweep,
    normalized_character,
    presentation_character,
    printed_character,
    printed_parameters,
    secant_character,
    secant_membership,
    self_euler,
    unique_repair,
    write_results,
)


class ThetaTwistTests(unittest.TestCase):
    def test_printed_d3_character(self) -> None:
        self.assertEqual(
            printed_character(3),
            (
                Fraction(1),
                Fraction(1),
                Fraction(-11, 2),
                Fraction(37, 6),
                Fraction(-71, 24),
            ),
        )

    def test_unique_repair_closes_for_declared_examples(self) -> None:
        for d in (3, 5, 9, 15):
            with self.subTest(d=d):
                repair = unique_repair(d)
                self.assertEqual(repair["unique_positive_d_real_integral_twist"], 3)
                self.assertTrue(repair["counts_close"])
                self.assertTrue(repair["coefficients_close"])
                self.assertEqual(corrected_character(d), secant_character(d, 3))

    def test_direct_presentation_uses_printed_counts(self) -> None:
        n, r = printed_parameters(3)
        self.assertEqual(presentation_character(n, r, 1), printed_character(3))
        self.assertEqual(presentation_character(n, r, 3), secant_character(3, 3))

    def test_membership_and_euler_separation(self) -> None:
        self.assertFalse(secant_membership(printed_character(3), 3)["member"])
        self.assertTrue(secant_membership(corrected_character(3), 3)["member"])
        self.assertEqual(self_euler(printed_character(3)), 288)
        self.assertEqual(self_euler(corrected_character(3)), 288)
        self.assertEqual(self_euler(secant_character(3, 1)), 96)

    def test_normalized_character_is_preserved(self) -> None:
        for d in (3, 5, 9):
            with self.subTest(d=d):
                self.assertEqual(
                    normalized_character(printed_character(d)),
                    normalized_character(corrected_character(d)),
                )

    def test_declared_natural_repair_space_has_one_nondual_survivor(self) -> None:
        sweep = natural_repair_sweep(3)
        candidates = {row["id"]: row for row in sweep["candidates"]}
        self.assertTrue(candidates["nondual_plus_2"]["secant_membership"]["member"])
        self.assertFalse(candidates["nondual_minus_4"]["secant_membership"]["member"])
        self.assertTrue(candidates["dual_minus_2"]["secant_membership"]["member"])
        self.assertFalse(candidates["dual_plus_4"]["secant_membership"]["member"])
        self.assertEqual(
            sweep["exact_alpha_plus_beta_autoequivalence_no_go"],
            {
                "printed_euler": 288,
                "target_euler": 96,
                "reason": "derived autoequivalences keep the Euler pairing unchanged",
            },
        )

    def test_hodge_witness_is_restored(self) -> None:
        printed = hodge_contraction(printed_character(3))
        corrected = hodge_contraction(corrected_character(3))
        self.assertEqual(printed["full_rank"], 12)
        self.assertEqual(printed["full_kernel_dimension"], 16)
        self.assertEqual(printed["mixed_kernel_vector_A_Bantisym_C"], [-7, -2, 1])
        self.assertEqual(corrected["mixed_kernel_vector_A_Bantisym_C"], [-3, 0, 1])

    def test_certificate_has_no_inferred_group_action(self) -> None:
        result = certificate(3)
        self.assertTrue(result["normalized_character_preserved"])
        self.assertIn("group action", result["claim_boundary"])
        self.assertNotIn("effective_group_order", json.dumps(result))

    def test_hash_manifest_closes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            paths = write_results(Path(directory), 3)
            manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
            artifact = manifest["artifacts"][0]
            self.assertEqual(
                hashlib.sha256(paths["certificate"].read_bytes()).hexdigest(),
                artifact["sha256"],
            )


if __name__ == "__main__":
    unittest.main()
