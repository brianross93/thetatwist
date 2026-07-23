import unittest

import numpy as np
from scipy import sparse

from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_residue_transport as transport,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_tree_contraction as tree,
)


class HSRank2KGeneratorResidueTransportTest(unittest.TestCase):
    def test_dual_frobenius_solver_is_exact(self):
        target = (7, 5)
        multiplier = (11, 13)
        transported = (
            np.asarray(target, dtype=np.int64).reshape(1, 2)
            @ transport.dual_matrix(multiplier)
        )[0] % transport.P
        self.assertEqual(
            transport.frobenius_multiplier(target, transported),
            multiplier,
        )

    def test_dual_factorization_is_exact(self):
        left = (3, 4)
        right = (5, 7)
        product = transport.dual_multiply(left, right)
        self.assertEqual(transport.dual_divide(product, left), right)
        self.assertTrue(
            np.array_equal(
                transport.dual_matrix(product),
                transport.dual_matrix(left) @ transport.dual_matrix(right)
                % transport.P,
            )
        )

    def test_orientation_coboundary_formula_is_one_law(self):
        orientation_scalar = (8, 0)
        source_orientation = (9, 4)
        target_orientation = (7, 12)
        top_graph = (5, 6)
        coboundary = transport.dual_divide(
            source_orientation, target_orientation
        )
        residue = transport.dual_multiply(
            transport.dual_multiply(orientation_scalar, coboundary),
            top_graph,
        )
        source_residue = (
            orientation_scalar[0] * source_orientation[1] % transport.P,
            orientation_scalar[0] * source_orientation[0] % transport.P,
        )
        target_residue = (
            target_orientation[1],
            target_orientation[0],
        )
        transported = (
            np.asarray(source_residue, dtype=np.int64).reshape(1, 2)
            @ transport.dual_matrix(top_graph)
        ) % transport.P
        replay = (
            np.asarray(target_residue, dtype=np.int64).reshape(1, 2)
            @ transport.dual_matrix(residue)
        ) % transport.P
        self.assertTrue(np.array_equal(transported, replay))

    def test_e0_constant_25_regression_hash_is_fixed(self):
        pairs = [(25, 0)] * 35
        self.assertEqual(
            transport.array_sha256(pairs),
            "be478c7b87e7eee0d5dda50be164c17ea0638f42a1efaa88b907eede5ee87aef",
        )

    def test_target_normalization_scales_p_without_mixing_q_columns(self):
        target_dimension = 4
        column_count = 28 * target_dimension
        raw_pairs = []
        dense_inputs = []
        for base in range(2):
            value = (
                np.arange(column_count, dtype=np.int64)
                + 3
                + 7 * base
            ) % transport.P
            tangent = (
                2 * np.arange(column_count, dtype=np.int64)
                + 5
                + 11 * base
            ) % transport.P
            raw_pairs.append(
                (
                    sparse.csr_matrix(value.reshape(1, -1)),
                    sparse.csr_matrix(tangent.reshape(1, -1)),
                )
            )
            dense_inputs.append((value, tangent))
        normalizations = [(5, 7), (9, 4)]
        transformed = tree.transform_target_lambda(
            raw_pairs,
            np.eye(target_dimension, dtype=np.int64),
            normalizations,
        )
        for base, ((value, tangent), (scale, derivative)) in enumerate(
            zip(dense_inputs, normalizations)
        ):
            output_value = transformed[base][0].toarray()[0, :column_count]
            output_tangent = transformed[base][1].toarray()[0, :column_count]
            self.assertTrue(
                np.array_equal(
                    output_value,
                    scale * value % transport.P,
                )
            )
            self.assertTrue(
                np.array_equal(
                    output_tangent,
                    (derivative * value + scale * tangent) % transport.P,
                )
            )


if __name__ == "__main__":
    unittest.main()
