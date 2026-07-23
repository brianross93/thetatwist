import os
import tempfile
import unittest
from pathlib import Path

from experiments.section12_instantiation import section12_pc_handoff_manifest as handoff


class Section12PCHandoffManifestTests(unittest.TestCase):
    def make_frozen_inputs(self, root):
        for index, relative_path in enumerate(handoff.FROZEN_INPUTS):
            path = root / relative_path
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(f"frozen-input-{index}\n".encode("ascii"))

    def test_emit_parse_and_verify_exact_selected_inventory(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_frozen_inputs(root)
            result_directory = (
                root / handoff.GENERATOR_RESULT_DIRECTORIES["e2"]
            )
            (result_directory / "nested").mkdir(parents=True)
            (result_directory / "z-last.json").write_text(
                '{"status":"PASS"}\n', encoding="utf-8"
            )
            (result_directory / "nested/a-first.bin").write_bytes(b"\x00\x01")

            records = handoff.build_records(root, ("e2",))
            paths = [record.relative_path for record in records]
            self.assertEqual(paths, sorted(paths))
            self.assertEqual(len(records), len(handoff.FROZEN_INPUTS) + 2)

            text = handoff.render_manifest(records)
            parsed = handoff.parse_manifest(text)
            self.assertEqual(parsed, records)
            count, total_bytes = handoff.verify_records(root, parsed, ("e2",))
            self.assertEqual(count, len(records))
            self.assertEqual(
                total_bytes, sum(record.byte_count for record in records)
            )

    def test_verification_rejects_payload_tampering(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.make_frozen_inputs(root)
            records = handoff.build_records(root)
            (root / handoff.FROZEN_INPUTS[0]).write_bytes(b"changed\n")
            with self.assertRaisesRegex(
                handoff.HandoffManifestError, "changed"
            ):
                handoff.verify_records(root, records)

    def test_parser_rejects_unsorted_and_unsafe_paths(self):
        first = handoff.ManifestRecord("0" * 64, 1, "z/file")
        second = handoff.ManifestRecord("1" * 64, 1, "a/file")
        with self.assertRaisesRegex(
            handoff.HandoffManifestError, "not sorted"
        ):
            handoff.parse_manifest(first.line() + "\n" + second.line() + "\n")
        with self.assertRaisesRegex(
            handoff.HandoffManifestError, "unsafe component"
        ):
            handoff.parse_manifest(
                handoff.ManifestRecord("0" * 64, 1, "../escape").line() + "\n"
            )

    def test_generator_symlink_is_rejected_without_following_it(self):
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as outside:
            root = Path(temporary)
            self.make_frozen_inputs(root)
            result_directory = (
                root / handoff.GENERATOR_RESULT_DIRECTORIES["e3"]
            )
            result_directory.mkdir(parents=True)
            target = Path(outside) / "outside.bin"
            target.write_bytes(b"must-not-be-read")
            link = result_directory / "outside-link.bin"
            try:
                os.symlink(target, link)
            except (OSError, NotImplementedError) as exc:
                self.skipTest(f"symlinks unavailable: {exc}")
            with self.assertRaisesRegex(
                handoff.HandoffManifestError, "symlinks are not allowed"
            ):
                handoff.build_records(root, ("e3",))


if __name__ == "__main__":
    unittest.main()
