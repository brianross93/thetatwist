import csv
import hashlib
import json
import tempfile
import unittest
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from experiments.section12_instantiation import (
    cm_bridge_hs_graph_koszul_j1_overlap_lifts as producer,
)


class HSGraphKoszulJ1OverlapLiftsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        certificate_path = (
            producer.OUTPUT_DIRECTORY / "j1_overlap_lift_certificate.json"
        )
        if not certificate_path.exists():
            producer.write_artifacts()
        cls.certificate = json.loads(certificate_path.read_text(encoding="utf-8"))

    @staticmethod
    def rows(name):
        with (producer.OUTPUT_DIRECTORY / name).open(
            newline="", encoding="utf-8"
        ) as handle:
            return list(csv.DictReader(handle))

    def test_all_endpoint_and_generator_chain_identities_closed(self):
        counts = self.certificate["counts"]
        self.assertEqual(counts["edges"], 192)
        self.assertEqual(counts["endpoint_maps"], 384)
        self.assertEqual(counts["source_endpoint_maps"], 192)
        self.assertEqual(counts["target_endpoint_maps"], 192)
        self.assertEqual(counts["graph_generator_chain_identities"], 3072)
        checks = self.certificate["exact_checks"]
        self.assertTrue(
            checks["all_3072_Laurent_graph_generator_residuals_are_zero"]
        )
        verification = self.rows("j1_endpoint_verification.csv")
        self.assertEqual(len(verification), 384)
        self.assertEqual(
            {row["status"] for row in verification},
            {"EXACT_CHAIN_IDENTITY_VERIFIED"},
        )
        self.assertEqual(
            {row["identity"] for row in verification},
            {"d_edge*J1=J0*d_chart"},
        )
        self.assertTrue(
            all(int(row["residual_laurent_terms"]) == 0 for row in verification)
        )

    def test_factored_program_is_the_exact_rational_changed_block(self):
        rows = self.rows("j1_factor_program.csv")
        self.assertEqual(len(rows), 3264)
        by_map = defaultdict(list)
        for row in rows:
            by_map[row["map_id"]].append(row)
        self.assertEqual(len(by_map), 384)
        counts = Counter(row["endpoint"] for row in rows)
        self.assertEqual(counts, {"source": 1536, "target": 1728})

        for map_rows in by_map.values():
            endpoint = map_rows[0]["endpoint"]
            changed = int(map_rows[0]["changed_factor"])
            u = 2 * changed
            v = u + 1
            if endpoint == "source":
                self.assertEqual(len(map_rows), 8)
                self.assertEqual(
                    {
                        (
                            int(row["domain_generator"]),
                            int(row["codomain_generator"]),
                            tuple(json.loads(row["laurent_exponents_source_coordinates"])),
                            int(row["scalar_mod_31"]),
                            row["coefficient_block_kind"],
                        )
                        for row in map_rows
                    },
                    {
                        (index, index, (0,) * 8, 1, "J0")
                        for index in range(8)
                    },
                )
                continue

            self.assertEqual(len(map_rows), 9)
            changed_rows = [
                row
                for row in map_rows
                if int(row["domain_generator"]) in (u, v)
            ]
            denominator = [0] * 8
            denominator[v] = -1
            denominator = tuple(denominator)
            self.assertEqual(
                {
                    (
                        int(row["domain_generator"]),
                        int(row["codomain_generator"]),
                        tuple(json.loads(row["laurent_exponents_source_coordinates"])),
                        int(row["scalar_mod_31"]),
                        row["coefficient_block_kind"],
                        int(row["chart_multiplication_index"]),
                    )
                    for row in changed_rows
                },
                {
                    (u, u, denominator, 1, "J0", -1),
                    (u, v, denominator, 30, "J0_TIMES_M_CHART", u),
                    (v, v, denominator, 30, "J0_TIMES_M_CHART", v),
                },
            )

    def test_every_coefficient_block_is_materialized_and_hash_bound(self):
        rows = self.rows("j1_factor_program.csv")
        with np.load(
            producer.OUTPUT_DIRECTORY / "j1_coefficient_blocks.npz"
        ) as blocks:
            self.assertEqual(len(blocks.files), 768)
            for row in rows:
                key = row["coefficient_block_key"]
                self.assertIn(key, blocks.files)
                self.assertEqual(
                    producer.array_sha256(blocks[key]),
                    row["coefficient_block_sha256"],
                )

    def test_full_shamash_boundary_is_not_overclaimed(self):
        pending = json.loads(
            (
                producer.OUTPUT_DIRECTORY
                / "pending_j2_j3_shamash_contract.json"
            ).read_text(encoding="utf-8")
        )
        self.assertTrue(
            pending["graph_exterior_observation"][
                "wedge2_and_wedge3_are_formally_determined_by_J1"
            ]
        )
        self.assertTrue(pending["missing_maps_are_not_zero"])
        self.assertIn(
            "changed curve entry",
            pending["known_regular_sequence_change"]["B"],
        )
        self.assertIn(
            "Laurent t is not",
            pending["known_regular_sequence_change"]["warning"],
        )
        self.assertIn("full J2", pending["not_emitted"])
        self.assertIn("full J3", pending["not_emitted"])
        self.assertIn("transported HS lambda", pending["not_emitted"])
        self.assertIn("Yoneda product", pending["not_emitted"])

    def test_temporary_regeneration_and_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            certificate, result = producer.write_artifacts(output)
            self.assertEqual(certificate["counts"]["input_files_bound"], 576)
            self.assertEqual(len(result["input_rows"]), 576)
            manifest = json.loads((output / "manifest.json").read_text())
            self.assertEqual(
                manifest["schema"],
                "section12-HS-graph-Koszul-J1-overlap-lifts-manifest-v1",
            )
            self.assertEqual(
                manifest["generator"]["sha256"],
                hashlib.sha256(Path(producer.__file__).read_bytes()).hexdigest(),
            )
            self.assertEqual(len(manifest["files"]), 6)
            for name, record in manifest["files"].items():
                path = output / name
                self.assertEqual(path.stat().st_size, record["bytes"])
                self.assertEqual(
                    hashlib.sha256(path.read_bytes()).hexdigest(),
                    record["sha256"],
                )


if __name__ == "__main__":
    unittest.main()
