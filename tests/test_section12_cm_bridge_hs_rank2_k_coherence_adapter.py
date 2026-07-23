import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_coherence_adapter as adapter,
)


def generator_states(*ready_labels):
    ready_labels = set(ready_labels)
    result = {}
    for label, exponents in adapter.driver.GENERATOR_EXPONENTS.items():
        ready = label in ready_labels
        result[label] = {
            "generator": label,
            "generator_exponents": tuple(exponents),
            "output_directory": f"/tmp/{label}",
            "preflight_available": False,
            "preflight_is_prerequisite": False,
            "artifact_status": (
                "LOADED_SUPPORT_IDS_BOUND" if ready else "ARTIFACT_ABSENT"
            ),
            "artifact_loaded": ready,
            "internal_files_verified": ready,
            "source_binding": (
                "CURRENT_SOURCE_BOUND" if ready else None
            ),
            "support_ids_bound": ready,
            "support_binding_count": (
                adapter.driver.GENERATOR_SUPPORTS[label] if ready else 0
            ),
            "support_bindings_sha256": (
                f"{label}-support-hash" if ready else None
            ),
            "ready_for_common_basis_rebase": ready,
            "common_basis_rebased": False,
            "numeric_edge_matrix_emitted": False,
            "validation_error": None,
        }
    return result


class HSRank2KCoherenceAdapterTest(unittest.TestCase):
    def load_e0_or_skip(self):
        try:
            return adapter.load_generator_artifact("e0")
        except adapter.GeneratorArtifactUnavailable:
            self.skipTest("the optional large e0 artifact is not present")

    def test_completed_manifest_and_lazy_csr_records_verify(self):
        artifact = self.load_e0_or_skip()
        self.assertEqual(artifact.generator, "e0")
        self.assertEqual(artifact.exponents, (1, 0, 0, 0))
        self.assertEqual(artifact.support_count, 35)
        self.assertEqual(artifact.overlap_length, 70)
        self.assertIn(
            artifact.source_binding,
            {"CURRENT_SOURCE_BOUND", "LEGACY_CURRENT_SOURCE_MISMATCH"},
        )
        self.assertEqual(
            set(artifact.verified_files),
            {
                "counit_pq_cache.npz",
                "overlap_adapted_hermite.npz",
                "tree_contracted_support_maps.npz",
                "tree_contraction_certificate.json",
            },
        )

        value, tangent = adapter.load_support_matrix_pair(
            artifact, 0, "a"
        )
        self.assertEqual(value.shape, (1, 687))
        self.assertEqual(tangent.shape, (1, 687))
        self.assertGreater(value.nnz, 0)
        self.assertGreater(tangent.nnz, 0)
        self.assertTrue(np.all((0 < value.data) & (value.data < 31)))
        self.assertTrue(np.all((0 < tangent.data) & (tangent.data < 31)))

    def test_completed_new_generator_pair_ledgers_use_canonical_hashes(self):
        checked = []
        for generator in ("e1", "e2", "e3"):
            output = Path(adapter.driver.OUTPUT_DIRECTORIES[generator])
            manifest_path = output / "manifest.json"
            maps_path = output / "tree_contracted_support_maps.npz"
            certificate_path = output / "tree_contraction_certificate.json"
            if not (
                manifest_path.is_file()
                and maps_path.is_file()
                and certificate_path.is_file()
            ):
                continue
            manifest = adapter.load_json(manifest_path)
            certificate = adapter.load_json(certificate_path)
            normalization = adapter._canonical_normalization_certificate(
                output,
                manifest_path,
                manifest,
                certificate,
                generator,
            )
            self.assertEqual(
                normalization["pair_ledger_hash_encoding"],
                adapter.PAIR_LEDGER_HASH_ENCODING,
            )
            with np.load(maps_path, allow_pickle=False) as archive:
                for name, hash_key in (
                    adapter.NORMALIZATION_HASH_KEYS.items()
                ):
                    self.assertEqual(
                        adapter.pair_ledger_sha256(archive[name]),
                        normalization[hash_key],
                    )
            checked.append(generator)
        if not checked:
            self.skipTest("no completed e1/e2/e3 artifacts are present")

    def test_pair_ledger_hash_is_independent_of_stored_integer_width(self):
        value = np.asarray([[1, 2], [30, 0]], dtype=np.uint8)
        widened = value.astype("<i8")
        self.assertEqual(
            adapter.pair_ledger_sha256(value),
            adapter.pair_ledger_sha256(widened),
        )
        self.assertNotEqual(
            adapter.pair_ledger_sha256(value),
            adapter.hashlib.sha256(value.tobytes()).hexdigest(),
        )

    def test_csr_reconstruction_reduces_and_combines_records_mod_31(self):
        records = np.asarray(
            [
                [0, 1, 30],
                [0, 1, 2],
                [1, 0, 62],
                [1, 1, -1],
            ],
            dtype=np.int64,
        )
        matrix = adapter.csr_from_records(records, (2, 2))
        self.assertEqual(
            matrix.toarray().tolist(),
            [[0, 1], [0, 30]],
        )
        with self.assertRaises(adapter.HSRank2KCoherenceAdapterError):
            adapter.csr_from_records([[2, 0, 1]], (2, 2))

    def test_support_slots_are_bound_to_unique_tuple_pair_keys(self):
        bindings = adapter.regenerate_support_bindings(
            "e0",
            adapter.driver.GENERATOR_SUPPORTS["e0"],
            adapter.driver.OUTPUT_DIRECTORIES["e0"],
        )
        lookup = adapter.support_binding_lookup(bindings)
        self.assertEqual(len(bindings), 35)
        self.assertEqual(len(lookup), 35)
        self.assertEqual(bindings[0]["artifact_support_slot"], "s00")
        self.assertEqual(
            bindings[0]["source_support"],
            (1, 1, 5, 1),
        )
        self.assertEqual(
            bindings[0]["target_support"],
            (2, 1, 5, 16),
        )
        self.assertIs(
            lookup[((1, 1, 5, 1), (2, 1, 5, 16))],
            bindings[0],
        )
        self.assertEqual(
            len(adapter.support_bindings_sha256(bindings)),
            64,
        )

    def test_exact_edge_cycle_and_square_registries_match_endpoints(self):
        edges = adapter.enumerate_typed_edges(generator_states("e0"))
        cycles = adapter.enumerate_order_three_cycles(edges)
        squares = adapter.enumerate_commutation_squares(edges)

        self.assertEqual(len(edges), 324)
        self.assertEqual(len(cycles), 108)
        self.assertEqual(len(squares), 486)
        self.assertEqual(
            sum(row["ready_for_common_basis_rebase"] for row in edges),
            81,
        )
        self.assertEqual(
            sum(row["generator_templates_available"] for row in cycles),
            27,
        )
        self.assertEqual(
            sum(row["generator_templates_available"] for row in squares),
            0,
        )
        self.assertEqual(edges[0]["target_K_exponents"], (1, 0, 0, 0))
        self.assertEqual(cycles[0]["vertex_K_indices"], (0, 27, 54, 0))
        self.assertEqual(cycles[0]["edge_job_ids"], (0, 108, 216))
        self.assertEqual(
            squares[0]["left_then_right_edge_job_ids"],
            (0, 109),
        )
        self.assertEqual(
            squares[0]["right_then_left_edge_job_ids"],
            (1, 36),
        )
        self.assertTrue(
            adapter.compare_committed_cayley_registry(edges)["matches"]
        )
        self.assertFalse(any(row["coherence_closed"] for row in cycles))
        self.assertFalse(any(row["coherence_closed"] for row in squares))

    def test_new_generator_artifacts_change_readiness_without_code_changes(self):
        expected = {
            (): (0, 0, 0),
            ("e0",): (81, 27, 0),
            ("e0", "e1"): (162, 54, 81),
            ("e0", "e1", "e2", "e3"): (324, 108, 486),
        }
        for labels, counts in expected.items():
            edges = adapter.enumerate_typed_edges(generator_states(*labels))
            cycles = adapter.enumerate_order_three_cycles(edges)
            squares = adapter.enumerate_commutation_squares(edges)
            actual = (
                sum(
                    row["ready_for_common_basis_rebase"]
                    for row in edges
                ),
                sum(
                    row["generator_templates_available"]
                    for row in cycles
                ),
                sum(
                    row["generator_templates_available"]
                    for row in squares
                ),
            )
            self.assertEqual(actual, counts)
            self.assertFalse(
                any(row["numeric_edge_matrix_emitted"] for row in edges)
            )
            self.assertFalse(any(row["coherence_closed"] for row in cycles))
            self.assertFalse(any(row["coherence_closed"] for row in squares))

    def test_partial_build_and_written_certificate_keep_numeric_scope_open(self):
        states = generator_states("e0")
        fake_bindings = [
            {
                "generator": "e0",
                "generator_exponents": (1, 0, 0, 0),
                "artifact_support_slot": "s00",
                "source_support": (1, 1, 5, 1),
                "target_support": (2, 1, 5, 16),
                "support_pair_key": (
                    (1, 1, 5, 1),
                    (2, 1, 5, 16),
                ),
                "support_bindings_sha256": "e0-support-hash",
            }
        ]
        with mock.patch.object(
            adapter,
            "_inspect_generators",
            return_value=(states, {}, fake_bindings),
        ):
            result = adapter.build()

        counts = result["certificate"]["counts"]
        self.assertEqual(counts["generator_artifacts_loaded"], 1)
        self.assertEqual(
            counts["typed_edges_ready_for_common_basis_rebase"], 81
        )
        self.assertEqual(counts["typed_edges_common_basis_rebased"], 0)
        self.assertEqual(counts["order_three_cycles_closed"], 0)
        self.assertEqual(
            counts["generator_commutation_squares_closed"], 0
        )
        boundary = result["certificate"]["numeric_claim_boundary"]
        self.assertFalse(boundary["typed_edge_common_bases_constructed"])
        self.assertFalse(boundary["cycle_matrix_products_computed"])
        self.assertFalse(boundary["square_matrix_products_computed"])
        self.assertFalse(boundary["global_Phi_or_Yoneda_attempted"])
        policy = result["certificate"]["availability_policy"]
        self.assertFalse(
            policy["missing_generator_artifacts_block_registry_emission"]
        )
        self.assertFalse(policy["optional_preflight_is_a_prerequisite"])

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            with mock.patch.object(
                adapter,
                "build",
                return_value=result,
            ):
                certificate = adapter.write_artifacts(output)
            self.assertEqual(
                certificate["schema"],
                adapter.ADAPTER_CERTIFICATE_SCHEMA,
            )
            self.assertTrue((output / "manifest.json").is_file())
            with (
                output / "z2_k_typed_edges.csv"
            ).open(newline="", encoding="utf-8") as handle:
                edge_rows = list(csv.DictReader(handle))
            with (
                output / "z2_k_order_three_cycles.csv"
            ).open(newline="", encoding="utf-8") as handle:
                cycle_rows = list(csv.DictReader(handle))
            with (
                output / "z2_k_generator_commutation_squares.csv"
            ).open(newline="", encoding="utf-8") as handle:
                square_rows = list(csv.DictReader(handle))
            persisted = json.loads(
                (
                    output / "k_coherence_adapter_certificate.json"
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(len(edge_rows), 324)
            self.assertEqual(len(cycle_rows), 108)
            self.assertEqual(len(square_rows), 486)
            self.assertEqual(
                persisted["counts"]["strict_coherence_equalities_computed"],
                0,
            )

    def test_absent_or_incomplete_generator_is_not_loaded(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            with self.assertRaises(adapter.GeneratorArtifactUnavailable):
                adapter.load_generator_artifact("e1", output)

            (output / "manifest.json").write_text(
                json.dumps(
                    {
                        "schema": adapter.TREE_MANIFEST_SCHEMA,
                        "files": {},
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaises(adapter.GeneratorArtifactInvalid):
                adapter.load_generator_artifact("e1", output)


if __name__ == "__main__":
    unittest.main()
