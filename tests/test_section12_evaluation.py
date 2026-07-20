from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from thetatwist.core import corrected_character
from thetatwist.section12_evaluation import (
    cech_character,
    cech_contract,
    cyclic_distance,
    evaluation_handoff,
    mixed_kernel_witnesses,
    ordinary_kernel_witnesses,
    section12_incidence,
    write_evaluation_handoff,
)


class IncidenceTests(unittest.TestCase):
    def test_d3_incidence_ledger(self) -> None:
        ledger = section12_incidence(3)
        self.assertEqual(ledger["component_count"], 6)
        self.assertEqual(ledger["adjacent_pairs"], 6)
        self.assertEqual(ledger["point_pairs"], 9)
        self.assertEqual(ledger["triple_points"], 6)
        self.assertEqual(ledger["isolated_pairs"], 3)
        self.assertEqual(ledger["total_isolated_nodes"], 72)
        self.assertEqual(ledger["normalized_nodes"], 60)
        self.assertEqual(ledger["unnormalized_isolated_nodes"], 12)

    def test_general_count_formulas_close(self) -> None:
        for d in (3, 5, 9, 15):
            with self.subTest(d=d):
                ledger = section12_incidence(d)
                n = (d + 9) // 2
                self.assertEqual(ledger["adjacent_pairs"], n)
                self.assertEqual(ledger["point_pairs"], n * (n - 3) // 2)
                self.assertEqual(ledger["triple_points"], n)
                self.assertEqual(ledger["isolated_pairs"], n * (n - 5) // 2)
                self.assertEqual(
                    ledger["total_isolated_nodes"], 3 * (d + 9) * (d - 1)
                )
                self.assertEqual(
                    ledger["unnormalized_isolated_nodes"], (d + 9) * (d - 2)
                )

    def test_cycle_distance(self) -> None:
        self.assertEqual(cyclic_distance(0, 5, 6), 1)
        self.assertEqual(cyclic_distance(0, 4, 6), 2)
        self.assertEqual(cyclic_distance(0, 3, 6), 3)


class DerivedContractTests(unittest.TestCase):
    def test_cech_ledger_recovers_corrected_character(self) -> None:
        for d in (3, 5, 9):
            with self.subTest(d=d):
                self.assertEqual(cech_character(d), corrected_character(d))
                self.assertTrue(cech_contract(d)["character_matches_alpha_plus_3beta"])

    def test_d3_normalized_resolution_keeps_twelve_nodes(self) -> None:
        normalized = cech_contract(3)["normalized_incidence_resolution"]
        self.assertEqual(normalized["N0"]["summands"], 6)
        self.assertEqual(normalized["N1"]["adjacent_curve_summands"], 6)
        self.assertEqual(
            normalized["N1"]["consecutive_four_theta_summands"], 6
        )
        self.assertEqual(normalized["N1"]["residual_isolated_node_length"], 12)
        self.assertEqual(normalized["N2"]["consecutive_four_theta_summands"], 6)

    def test_kernel_split_is_exactly_ten_plus_six(self) -> None:
        ordinary = ordinary_kernel_witnesses()
        mixed = mixed_kernel_witnesses(3)
        self.assertEqual(len(ordinary), 10)
        self.assertEqual(len(mixed), 6)
        self.assertEqual(len(set(ordinary + mixed)), 16)
        self.assertIn("B_12+B_21", ordinary)
        self.assertIn("C_12-3A_12", mixed)

    def test_handoff_keeps_all_six_mixed_checks(self) -> None:
        handoff = evaluation_handoff(3)
        checks = handoff["completion_test"]["required_equalities"]
        self.assertEqual(len(checks), 6)
        self.assertIn("Do not reduce", handoff["symmetry_boundary"])
        self.assertNotIn("W_2", json.dumps(handoff))

    def test_parameter_incidence_closure_is_named_and_exact(self) -> None:
        closure = evaluation_handoff(3)["parameter_incidence_closure"]
        self.assertEqual(closure["name"], "parameter-incidence closure check")
        self.assertEqual(closure["symbolic_schema_status"], "passed")
        self.assertEqual(closure["instantiated_fiber_status"], "not yet computed")
        self.assertEqual(closure["component_parameter"], 6)
        self.assertEqual(closure["realized_surface_count"], 6)
        self.assertEqual(closure["normalization_parameter"], 60)
        self.assertEqual(closure["realized_normalized_node_count"], 60)
        self.assertEqual(closure["total_isolated_nodes"], 72)
        self.assertEqual(closure["realized_residual_node_count"], 12)
        self.assertEqual(closure["strata"]["adjacent_triple_theta_curves_T_i"], 6)
        self.assertEqual(closure["strata"]["consecutive_four_theta_schemes_Q_i"], 6)
        self.assertEqual(closure["strata"]["isolated_four_theta_schemes"], 3)

    def test_both_mixed_outcomes_are_declared(self) -> None:
        outcomes = evaluation_handoff(3)["preregistered_outcomes"]
        self.assertIn("all_six_zero", outcomes)
        self.assertIn("any_one_nonzero", outcomes)
        self.assertIn("does not assert surjectivity", outcomes["all_six_zero"])
        self.assertIn("recorded node choice", outcomes["any_one_nonzero"])
        self.assertIn("No mathematical verdict", outcomes["integrity_failure"])

    def test_execution_strategy_keeps_level_and_symmetry_boundaries(self) -> None:
        handoff = evaluation_handoff(3)
        strategy = handoff["execution_strategy"]
        symmetry = handoff["symmetry_analysis"]
        self.assertIn("not automatically a level-3", strategy["order_six_translation"])
        self.assertIn("level-6", strategy["order_six_translation"])
        self.assertIn("does not by itself", strategy["finite_field_selection"])
        self.assertTrue(any("6g=0" in check for check in strategy["fiber_checks"]))
        self.assertIn("do not merge", symmetry["translation_use"])
        self.assertIn("nontrivial linear action", symmetry["valid_reduction"])

    def test_manifest_closes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            paths = write_evaluation_handoff(Path(directory), 3)
            manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
            artifact = manifest["artifacts"][0]
            self.assertEqual(
                hashlib.sha256(paths["handoff"].read_bytes()).hexdigest(),
                artifact["sha256"],
            )


if __name__ == "__main__":
    unittest.main()
