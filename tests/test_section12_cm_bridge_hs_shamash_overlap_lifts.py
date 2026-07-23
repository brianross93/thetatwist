import csv
import hashlib
import json
import unittest
from pathlib import Path

import numpy as np

from experiments.section12_instantiation import (
    cm_bridge_hs_shamash_overlap_lifts as producer,
)


class HSShamashOverlapLiftsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.output = producer.OUTPUT_DIRECTORY
        certificate_path = cls.output / "shamash_overlap_lifts_certificate.json"
        if not certificate_path.exists():
            raise RuntimeError(
                "run cm_bridge_hs_shamash_overlap_lifts.py before this test"
            )
        cls.certificate = json.loads(
            certificate_path.read_text(encoding="utf-8")
        )

    @classmethod
    def csv_rows(cls, name):
        with (cls.output / name).open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    def test_complete_census_and_exact_component_ledgers(self):
        counts = self.certificate["counts"]
        self.assertEqual(counts["surfaces"], 6)
        self.assertEqual(counts["edges"], 192)
        self.assertEqual(counts["endpoint_maps"], 384)
        self.assertEqual(counts["source_endpoint_maps"], 192)
        self.assertEqual(counts["target_endpoint_maps"], 192)
        self.assertEqual(counts["overlap_support_points"], 11665)
        self.assertEqual(
            counts["source_equation_support_zero_checks"], 69990
        )
        self.assertEqual(counts["target_equation_support_checks"], 69990)
        self.assertEqual(counts["Hermite_J0_selector_checks"], 384)
        self.assertEqual(
            counts["rational_value_and_jet_transition_checks"], 11665
        )
        self.assertEqual(
            counts["full_J2_component_identities"], 7 * 11665
        )
        self.assertEqual(
            counts["full_J3_component_identities"], 13 * 11665
        )
        self.assertTrue(all(self.certificate["exact_checks"].values()))

    def test_all_endpoint_rows_and_edge_artifacts_are_hash_bound(self):
        registry = self.csv_rows("endpoint_registry.csv")
        self.assertEqual(len(registry), 384)
        self.assertEqual(
            {row["status"] for row in registry},
            {
                "EXACT_J0_J1_J2_J3_SOURCE_ZERO_CONTROL",
                "EXACT_J0_J1_J2_J3_TARGET_SHAMASH_COMPARISON",
            },
        )
        by_edge = {}
        for row in registry:
            key = (row["surface"], int(row["edge_id"]))
            by_edge.setdefault(key, []).append(row)
        self.assertEqual(len(by_edge), 192)
        for rows in by_edge.values():
            self.assertEqual({row["endpoint"] for row in rows}, {"source", "target"})
            binary = self.output / rows[0]["binary_artifact"]
            summary = self.output / rows[0]["summary_artifact"]
            self.assertEqual(
                hashlib.sha256(binary.read_bytes()).hexdigest(),
                rows[0]["binary_sha256"],
            )
            self.assertEqual(
                hashlib.sha256(summary.read_bytes()).hexdigest(),
                rows[0]["summary_sha256"],
            )

    def test_representative_materialized_maps_reconstruct_exactly(self):
        data = producer.load_problem_data()
        equations, _rows = producer._surface_equation_payload(data)
        overlap = data["overlaps"][(0, 0)]
        source = data["charts"][(0, overlap["source"])]
        target = data["charts"][(0, overlap["target"])]
        source_position = int(overlap["source_positions"][0])
        target_position = int(overlap["target_positions"][0])
        expected = producer.local_support_maps(
            source["values"][source_position],
            source["tangents"][source_position],
            target["values"][target_position],
            target["tangents"][target_position],
            equations[(0, overlap["source"])][0],
            equations[(0, overlap["target"])][0],
            overlap["changed"],
        )
        path = self.output / "Z_0_E00_target_shamash.npz"
        with np.load(path) as artifact:
            graph = np.asarray(
                artifact["target_graph_J0_J1_J2_J3_terms"], dtype=np.int64
            )
            gamma = np.asarray(artifact["target_Gamma_terms"], dtype=np.int64)
            omega = np.asarray(artifact["target_Omega_terms"], dtype=np.int64)
            expected_b = np.asarray(
                [
                    producer.equation_change_exponent(equation, overlap["changed"])
                    for equation in range(6)
                ],
                dtype=np.int64,
            )
            self.assertTrue(
                np.array_equal(
                    artifact["B_diagonal_laurent_exponents"], expected_b
                )
            )

        for degree in range(4):
            records = graph[(graph[:, 0] == 0) & (graph[:, 1] == degree)]
            reconstructed = producer.deserialize_block_map(
                records, degree, degree, 2
            )
            producer.assert_block_maps_equal(
                reconstructed,
                expected["graph_maps"][degree],
                f"stored graph J{degree}",
            )
        for equation in range(6):
            gamma_records = gamma[
                (gamma[:, 0] == 0) & (gamma[:, 1] == equation)
            ]
            omega_records = omega[
                (omega[:, 0] == 0) & (omega[:, 1] == equation)
            ]
            producer.assert_block_maps_equal(
                producer.deserialize_block_map(gamma_records, 2, 0, 2),
                expected["gamma"][equation],
                f"stored Gamma_{equation}",
            )
            producer.assert_block_maps_equal(
                producer.deserialize_block_map(omega_records, 3, 1, 2),
                expected["omega"][equation],
                f"stored Omega_{equation}",
            )

    def test_program_and_next_contract_preserve_the_claim_boundary(self):
        program = json.loads(
            (self.output / "comparison_program.json").read_text()
        )
        self.assertEqual(
            program["maps"]["J2"], "[[J2K,Gamma],[0,B tensor J0]]"
        )
        self.assertEqual(
            program["maps"]["J3"], "[[J3K,Omega],[0,B tensor J1K]]"
        )
        self.assertIn("H_edge^-1", program["jet_to_standard_basis_reconstruction"]["degree_q_graph"])
        pending = json.loads(
            (self.output / "higher_descent_next_contract.json").read_text()
        )
        self.assertIn(
            "A1*A3+A2^2+A3*A1+partial*A4+A4*partial=0",
            pending["totalization_relations"],
        )
        self.assertIn("square coherence", pending["not_yet_claimed"])
        self.assertIn("the Hodge conjecture", pending["not_yet_claimed"])

    def test_manifest_and_current_j1_generator_are_bound(self):
        manifest = json.loads((self.output / "manifest.json").read_text())
        self.assertEqual(
            manifest["generator"]["sha256"],
            hashlib.sha256(Path(producer.__file__).read_bytes()).hexdigest(),
        )
        self.assertEqual(len(manifest["files"]), 390)
        for name, record in manifest["files"].items():
            path = self.output / name
            self.assertEqual(path.stat().st_size, record["bytes"])
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(),
                record["sha256"],
            )
        j1_manifest_path = (
            producer.j1_lifts.OUTPUT_DIRECTORY / "manifest.json"
        )
        j1_manifest = json.loads(j1_manifest_path.read_text())
        self.assertEqual(
            j1_manifest["generator"]["sha256"],
            hashlib.sha256(
                Path(producer.j1_lifts.__file__).read_bytes()
            ).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()
