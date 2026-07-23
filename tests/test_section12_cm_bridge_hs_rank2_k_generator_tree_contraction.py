import hashlib
import json
import unittest
from pathlib import Path

import numpy as np
from scipy import sparse

from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_tree_contraction as contraction,
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def csr_from_records(records, shape):
    records = np.asarray(records, dtype=np.int64)
    if not len(records):
        return sparse.csr_matrix(shape, dtype=np.int64)
    if (
        np.any(records[:, 0] < 0)
        or np.any(records[:, 0] >= shape[0])
        or np.any(records[:, 1] < 0)
        or np.any(records[:, 1] >= shape[1])
    ):
        raise AssertionError("a persisted CSR coordinate left its declared shape")
    matrix = sparse.csr_matrix(
        (
            records[:, 2] % contraction.P,
            (records[:, 0], records[:, 1]),
        ),
        shape=shape,
        dtype=np.int64,
    )
    matrix.data %= contraction.P
    matrix.eliminate_zeros()
    return matrix


class HSRank2KGeneratorTreeContractionArtifactTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.output = contraction.OUTPUT_DIRECTORY
        cls.certificate_path = cls.output / "tree_contraction_certificate.json"
        cls.maps_path = cls.output / "tree_contracted_support_maps.npz"
        cls.manifest_path = cls.output / "manifest.json"
        required = (
            cls.certificate_path,
            cls.maps_path,
            cls.manifest_path,
        )
        if not all(path.is_file() for path in required):
            raise unittest.SkipTest(
                "the large one-generator tree-contraction artifact is not present"
            )
        cls.certificate = json.loads(cls.certificate_path.read_text())
        cls.manifest = json.loads(cls.manifest_path.read_text())

    def test_certificate_schema_counts_shapes_and_scope(self):
        certificate = self.certificate
        self.assertEqual(
            certificate["schema"],
            "section12-HS-rank2-K-generator-tree-contraction-v3",
        )
        self.assertEqual(certificate["field"], "F_31")
        self.assertEqual(certificate["K_generator_exponents"], [1, 0, 0, 0])
        self.assertIn("one-generator", certificate["status"])
        self.assertIn("finite-overlap", certificate["status"])
        self.assertIn(
            "strict-Laurent-counit-lift-not-inferred",
            certificate["status"],
        )
        counts = certificate["counts"]
        self.assertEqual(counts["overlap_supports"], 35)
        self.assertEqual(counts["overlap_length"], 70)
        self.assertEqual(counts["target_length"], 98)
        self.assertEqual(counts["overlap_tree_ranks"], [1, 491, 1966, 6866])
        self.assertEqual(counts["target_tree_ranks"], [1, 687, 2750, 9610])
        self.assertEqual(counts["J1_shape"], [491, 687])
        self.assertEqual(counts["J2_shape"], [1966, 2750])
        self.assertEqual(counts["J3_shape"], [6866, 9610])
        self.assertEqual(counts["L0_shape"], [492, 688])
        self.assertEqual(counts["factorized_D3_chain_block_identities"], 455)
        checks = certificate["exact_checks"]
        self.assertEqual(checks["target_D2_value_rank_all_supports"], [685])
        self.assertEqual(
            checks["target_D2_free_directions_all_supports"], [2]
        )
        self.assertEqual(
            checks[
                "finite_overlap_lambda_edge_J2_minus_w_lambda_target_equals_a_D2_supports"
            ],
            35,
        )
        self.assertEqual(
            checks[
                "factorized_rectangular_D3_J3_equals_J2_D3_block_identities"
            ],
            455,
        )
        boundary = certificate["claim_boundary"]
        self.assertIn("factorized J3", boundary["proved"])
        self.assertNotIn("J3", boundary["not_claimed"])
        for unclaimed in (
            "strict Laurent/surface-ring",
            "other K generators",
            "K coherence",
            "cross-RHom/Yoneda",
            "Question 11.4",
            "Hodge",
        ):
            self.assertIn(unclaimed, boundary["not_claimed"])
        self.assertEqual(
            certificate["counit"]["strict_Laurent_tree_lift_status"],
            "OPEN",
        )

    def test_legacy_manifest_and_every_artifact_remain_internally_bound(self):
        manifest = self.manifest
        self.assertEqual(
            manifest["schema"],
            "section12-HS-rank2-K-generator-tree-contraction-manifest-v3",
        )
        generator = contraction.ROOT / manifest["generator"]["path"]
        self.assertEqual(generator.resolve(), Path(contraction.__file__).resolve())
        self.assertEqual(
            manifest["generator"]["sha256"],
            "f2848e704cf66ea5b99f1236f40ce21fedfeb89a5b90ecbd97ffaeeea4973cf7",
        )
        self.assertNotEqual(
            manifest["generator"]["sha256"],
            sha256(generator),
            "the e0 artifact is intentionally classified as legacy until the "
            "generic producer is final and e0 is refreshed",
        )

        expected_files = {
            "tree_contracted_support_maps.npz",
            "overlap_adapted_hermite.npz",
            "counit_pq_cache.npz",
            "tree_contraction_certificate.json",
        }
        self.assertEqual(set(manifest["files"]), expected_files)
        for name, binding in manifest["files"].items():
            path = self.output / name
            self.assertTrue(path.is_file())
            self.assertEqual(binding["bytes"], path.stat().st_size)
            self.assertEqual(binding["sha256"], sha256(path))

    def test_factorized_npz_has_exact_keys_for_all_35_supports(self):
        metadata = {
            "field",
            "generator_exponents",
            "J1_shape",
            "J2_shape",
            "J3_shape",
            "L0_shape",
        }
        matrix_factors = ("J1", "J2", "J3_G3", "L0", "L1", "a")
        expected = set(metadata)
        for support in range(35):
            prefix = f"s{support:02d}"
            for name in matrix_factors:
                for component in ("value", "tangent"):
                    expected.add(f"{prefix}_{name}_{component}")
            for equation in range(contraction.EQUATION_COUNT):
                for component in ("value", "tangent"):
                    expected.add(
                        f"{prefix}_J3_Omega_i1_e{equation}_{component}"
                    )
            expected.add(f"{prefix}_J3_B_value")
            expected.add(f"{prefix}_J3_B_tangent")

        with np.load(self.maps_path) as data:
            self.assertEqual(set(data.files), expected)
            self.assertEqual(data["field"].tolist(), [contraction.P])
            self.assertEqual(data["generator_exponents"].tolist(), [1, 0, 0, 0])
            self.assertEqual(data["J1_shape"].tolist(), [491, 687])
            self.assertEqual(data["J2_shape"].tolist(), [1966, 2750])
            self.assertEqual(data["J3_shape"].tolist(), [6866, 9610])
            self.assertEqual(data["L0_shape"].tolist(), [492, 688])
            for support in range(35):
                prefix = f"s{support:02d}"
                self.assertEqual(
                    data[f"{prefix}_J3_B_value"].shape,
                    (contraction.EQUATION_COUNT,),
                )
                self.assertEqual(
                    data[f"{prefix}_J3_B_tangent"].shape,
                    (contraction.EQUATION_COUNT,),
                )
                a_value = data[f"{prefix}_a_value"]
                a_tangent = data[f"{prefix}_a_tangent"]
                self.assertGreater(len(a_value), 0)
                self.assertGreater(len(a_tangent), 0)

    def test_support_zero_factorization_reconstructs_the_J3_hash(self):
        edge_f1, target_f1 = 491, 687
        edge_k3, target_k3 = 56 * 70, 56 * 98
        with np.load(self.maps_path) as data:
            self.assertEqual(
                tuple(map(int, data["J3_shape"])), (6866, 9610)
            )
            j1 = (
                csr_from_records(
                    data["s00_J1_value"], (edge_f1, target_f1)
                ),
                csr_from_records(
                    data["s00_J1_tangent"], (edge_f1, target_f1)
                ),
            )
            graph_j3 = (
                csr_from_records(
                    data["s00_J3_G3_value"], (edge_k3, target_k3)
                ),
                csr_from_records(
                    data["s00_J3_G3_tangent"], (edge_k3, target_k3)
                ),
            )
            omega_columns = []
            for equation in range(contraction.EQUATION_COUNT):
                omega_columns.append(
                    (
                        csr_from_records(
                            data[
                                f"s00_J3_Omega_i1_e{equation}_value"
                            ],
                            (edge_k3, target_f1),
                        ),
                        csr_from_records(
                            data[
                                f"s00_J3_Omega_i1_e{equation}_tangent"
                            ],
                            (edge_k3, target_f1),
                        ),
                    )
                )
            b_value = np.asarray(data["s00_J3_B_value"], dtype=np.int64)
            b_tangent = np.asarray(
                data["s00_J3_B_tangent"], dtype=np.int64
            )
            self.assertEqual(len(omega_columns), contraction.EQUATION_COUNT)
            self.assertEqual(
                b_value.shape, (contraction.EQUATION_COUNT,)
            )
            self.assertEqual(
                b_tangent.shape, (contraction.EQUATION_COUNT,)
            )
            b_j1 = [
                contraction.pair_left_scalar(
                    j1, (int(b_value[equation]), int(b_tangent[equation]))
                )
                for equation in range(contraction.EQUATION_COUNT)
            ]
            a_value = csr_from_records(
                data["s00_a_value"], (1, target_f1)
            )
            a_tangent = csr_from_records(
                data["s00_a_tangent"], (1, target_f1)
            )
        self.assertGreater(a_value.nnz, 0)
        self.assertGreater(a_tangent.nnz, 0)
        j3 = contraction.assemble_j3(graph_j3, omega_columns, b_j1)
        self.assertEqual(j3[0].shape, (6866, 9610))
        self.assertEqual(
            contraction.pair_hash(j3),
            self.certificate["residual_hashes"]["factorized_J3"][0],
        )


if __name__ == "__main__":
    unittest.main()
