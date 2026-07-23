import hashlib
import json
import unittest
from pathlib import Path

import numpy as np

from experiments.section12_instantiation import (
    cm_bridge_hs_lambda_overlap_transport as producer,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_lambda_overlap_transport_all_edges as all_edges,
)


class HSLambdaOverlapTransportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.output = producer.OUTPUT_DIRECTORY
        cls.certificate_path = (
            cls.output / "Z_0_E02_lambda_transport_certificate.json"
        )
        if not cls.certificate_path.exists():
            raise RuntimeError(
                "run cm_bridge_hs_lambda_overlap_transport.py before this test"
            )
        cls.certificate = json.loads(
            cls.certificate_path.read_text(encoding="utf-8")
        )

    def test_representative_edge_census_and_single_cochain_identity(self):
        counts = self.certificate["counts"]
        self.assertEqual(counts["edge_support_count"], 59)
        self.assertEqual(counts["computed_support_count"], 59)
        self.assertEqual(len(self.certificate["support_ledger"]), 59)
        for row in self.certificate["support_ledger"]:
            self.assertTrue(
                row["single_a_solves_K2_block_in_Hermite_quotient"]
            )
            self.assertTrue(
                row[
                    "same_a_solves_all_six_H0_blocks_in_Hermite_quotient"
                ]
            )
            self.assertTrue(row["lambda_Gamma_blocks_included"])
            self.assertTrue(
                row["hermite_residual_ledger"][
                    "all_Hermite_residuals_zero"
                ]
            )
            self.assertEqual(
                row["hermite_residual_ledger"]["verified_A_columns"], 34
            )
            self.assertEqual(
                row["hermite_residual_ledger"][
                    "verified_F31_scalar_coordinates"
                ],
                68,
            )

    def test_every_system_is_exact_and_canonically_solved(self):
        for row in self.certificate["support_ledger"]:
            system = row["linear_system"]
            self.assertEqual(system["matrix_shape"], [68, 16])
            self.assertEqual(system["A_unknowns"], 8)
            self.assertEqual(system["A_target_columns"], 34)
            self.assertEqual(system["F31_unknowns"], 16)
            self.assertEqual(system["F31_scalar_equations"], 68)
            self.assertEqual(system["rank"] + system["nullity"], 16)
            self.assertEqual(len(system["pivot_columns"]), system["rank"])
            self.assertEqual(len(system["free_columns"]), system["nullity"])
            self.assertTrue(system["canonical_free_variables_zero"])

    def test_orientation_and_weight_are_the_certified_values(self):
        orientation = self.certificate["orientation_normalization"]
        self.assertEqual(
            orientation["epsilon_C"], "(-1)^popcount(chart_bits)"
        )
        self.assertEqual(
            orientation["normalized_tau_relation"],
            "tau_s J8=t phi(tau_t)",
        )
        self.assertEqual(orientation["det_B"], "t^5")
        self.assertEqual(orientation["O4H_factor"], "t/det(B)=t^-4=h^4")

    def test_binary_term_schema_preserves_wedge_degree(self):
        path = self.output / self.certificate["binary_artifact"]["path"]
        with np.load(path) as artifact:
            self.assertTrue(
                np.array_equal(
                    artifact["support_positions"], np.arange(59, dtype=np.uint16)
                )
            )
            a_terms = np.asarray(
                artifact["a_target_source_coordinate_terms"], dtype=np.int64
            )
            self.assertEqual(a_terms.shape[1], 13)
            self.assertEqual(set(a_terms[:, 1]), {1})
            delta_terms = np.asarray(
                artifact["delta_K2_target_terms"], dtype=np.int64
            )
            self.assertEqual(delta_terms.shape[1], 13)
            self.assertEqual(set(delta_terms[:, 1]), {2})

    def test_manifest_and_artifact_hashes_verify(self):
        verified = producer.verify_artifacts(surface=0, edge_id=2)
        self.assertEqual(verified["counts"]["computed_support_count"], 59)
        manifest = json.loads(
            (self.output / "manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            manifest["producer_sha256"],
            hashlib.sha256(Path(producer.__file__).read_bytes()).hexdigest(),
        )
        for name, digest in manifest["artifacts"].items():
            self.assertEqual(
                hashlib.sha256((self.output / name).read_bytes()).hexdigest(),
                digest,
            )

    def test_all_edge_driver_plan_and_completed_registry_agree(self):
        plan = all_edges.edge_plan()
        self.assertEqual(len(plan), 192)
        self.assertEqual(sum(row["support_count"] for row in plan), 11665)
        if all_edges.OUTPUT_DIRECTORY.exists():
            registry = json.loads(
                (all_edges.OUTPUT_DIRECTORY / "edge_registry.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(registry["edge_count"], 192)
            self.assertEqual(registry["support_count"], 11665)

    def test_rref_solver_uses_zero_for_free_variables(self):
        system = np.asarray([[1, 1, 0], [0, 0, 1]], dtype=np.int64)
        rhs = np.asarray([7, 9], dtype=np.int64)
        solution, ledger = producer._rref_solve(system, rhs)
        self.assertTrue(np.array_equal(solution, np.asarray([7, 0, 9])))
        self.assertEqual(ledger["pivot_columns"], [0, 2])
        self.assertEqual(ledger["free_columns"], [1])


if __name__ == "__main__":
    unittest.main()
