"""Emit and verify the exact Section 12 PC-handoff file manifest.

The manifest always contains the seventeen frozen inputs needed by the
standard-generator tree-contraction driver.  Callers may additionally select
one or more generator result directories.  Records have the form

    sha256  bytes  repository-relative/path

and are sorted by repository-relative path.  Payload symlinks are rejected;
the utility never follows them.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import stat
import sys
from dataclasses import dataclass
from pathlib import Path


FROZEN_INPUTS = (
    "results/section12_instantiation/cm_bridge_global_hs_chart_resolution/"
    "global_hs_chart_resolution_certificate.json",
    "results/section12_instantiation/cm_bridge_global_hs_chart_resolution/"
    "degree2_global_generator_coefficients.coo.csv",
    "results/section12_instantiation/cm_bridge_global_hs_chart_resolution/"
    "product_theta_atlas.csv",
    "results/section12_instantiation/cm_bridge_global_hs_chart_resolution/"
    "product_theta_atlas_edges.csv",
    "results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/"
    "six_hs_local_seeds.csv",
    "results/section12_instantiation/cm_bridge_six_hs_local_seed_emitter/"
    "six_hs_support_K_orbits.csv",
    "results/section12_instantiation/descent_section/"
    "descent_section_certificate.json",
    "results/section12_instantiation/independent_incidence_fiber/"
    "candidate_section_coefficients.csv",
    "results/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit/"
    "K_degree2_action_elements.csv",
    "results/section12_instantiation/cm_bridge_hs_chart_local_syzygy_audit/"
    "K_stable_product_chart_refinement.csv",
    "results/section12_instantiation/cm_bridge_k_equivariant_edge_audit/"
    "k_refined_chart_intersections.csv",
    "results/section12_instantiation/cm_bridge_k_equivariant_edge_audit/"
    "k_refined_edge_pullbacks.csv",
    "results/section12_instantiation/"
    "cm_bridge_hs_extension_from_surface_presentation/"
    "hs_extension_certificate.json",
    "results/section12_instantiation/"
    "cm_bridge_hs_extension_from_surface_presentation/"
    "Z_2_C0010_hs_extension.npz",
    "results/section12_instantiation/cm_bridge_hs_border_to_local_frames/"
    "hs_completed_lci_supports.csv",
    "results/section12_instantiation/cm_bridge_iw_hermite_border_basis/"
    "Z_2_C0000.json",
    "results/section12_instantiation/cm_bridge_iw_hermite_border_basis/"
    "Z_2_C0010.json",
)

GENERATOR_RESULT_DIRECTORIES = {
    "e0": (
        "results/section12_instantiation/"
        "cm_bridge_hs_rank2_k_generator_tree_contraction"
    ),
    "e1": (
        "results/section12_instantiation/"
        "cm_bridge_hs_rank2_k_generator_tree_contraction_e1"
    ),
    "e2": (
        "results/section12_instantiation/"
        "cm_bridge_hs_rank2_k_generator_tree_contraction_e2"
    ),
    "e3": (
        "results/section12_instantiation/"
        "cm_bridge_hs_rank2_k_generator_tree_contraction_e3"
    ),
}

_RECORD_PATTERN = re.compile(r"^([0-9a-f]{64})  ([0-9]+)  (.+)$")
_READ_CHUNK_BYTES = 1024 * 1024


class HandoffManifestError(RuntimeError):
    """A manifest or selected payload is missing, unsafe, or inconsistent."""


@dataclass(frozen=True)
class ManifestRecord:
    sha256: str
    byte_count: int
    relative_path: str

    def line(self) -> str:
        return f"{self.sha256}  {self.byte_count}  {self.relative_path}"


def repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _relative_parts(relative_path: str) -> tuple[str, ...]:
    if not relative_path or "\\" in relative_path or relative_path.startswith("/"):
        raise HandoffManifestError(
            f"manifest path is not repository-relative POSIX text: {relative_path!r}"
        )
    parts = tuple(relative_path.split("/"))
    if any(part in {"", ".", ".."} for part in parts):
        raise HandoffManifestError(
            f"manifest path contains an unsafe component: {relative_path!r}"
        )
    return parts


def _checked_payload_path(
    root: Path, relative_path: str, expected_kind: str
) -> tuple[Path, os.stat_result]:
    parts = _relative_parts(relative_path)
    current = Path(root)
    try:
        root_status = os.lstat(current)
    except OSError as exc:
        raise HandoffManifestError(f"cannot inspect repository root {current}: {exc}") from exc
    if stat.S_ISLNK(root_status.st_mode) or not stat.S_ISDIR(root_status.st_mode):
        raise HandoffManifestError(f"repository root is not a real directory: {current}")

    status = root_status
    for position, part in enumerate(parts):
        current = current / part
        try:
            status = os.lstat(current)
        except OSError as exc:
            raise HandoffManifestError(
                f"required payload path is unavailable: {relative_path}: {exc}"
            ) from exc
        if stat.S_ISLNK(status.st_mode):
            raise HandoffManifestError(
                f"payload symlinks are not allowed: {relative_path}"
            )
        if position < len(parts) - 1 and not stat.S_ISDIR(status.st_mode):
            raise HandoffManifestError(
                f"payload parent is not a directory: {relative_path}"
            )

    if expected_kind == "file" and not stat.S_ISREG(status.st_mode):
        raise HandoffManifestError(f"payload is not a regular file: {relative_path}")
    if expected_kind == "directory" and not stat.S_ISDIR(status.st_mode):
        raise HandoffManifestError(f"payload is not a directory: {relative_path}")
    return current, status


def _directory_files(root: Path, relative_directory: str) -> list[str]:
    directory, _status = _checked_payload_path(
        root, relative_directory, expected_kind="directory"
    )
    files: list[str] = []

    def visit(path: Path, relative_path: str) -> None:
        try:
            entries = sorted(os.scandir(path), key=lambda entry: entry.name)
        except OSError as exc:
            raise HandoffManifestError(
                f"cannot enumerate generator result directory {relative_path}: {exc}"
            ) from exc
        for entry in entries:
            child_relative = f"{relative_path}/{entry.name}"
            try:
                child_status = entry.stat(follow_symlinks=False)
            except OSError as exc:
                raise HandoffManifestError(
                    f"cannot inspect generator result entry {child_relative}: {exc}"
                ) from exc
            if stat.S_ISLNK(child_status.st_mode):
                raise HandoffManifestError(
                    f"payload symlinks are not allowed: {child_relative}"
                )
            if stat.S_ISDIR(child_status.st_mode):
                visit(Path(entry.path), child_relative)
            elif stat.S_ISREG(child_status.st_mode):
                files.append(child_relative)
            else:
                raise HandoffManifestError(
                    f"payload entry is not a regular file or directory: {child_relative}"
                )

    visit(directory, relative_directory)
    return files


def inventory_paths(root: Path, generators: tuple[str, ...] = ()) -> tuple[str, ...]:
    unknown = sorted(set(generators) - set(GENERATOR_RESULT_DIRECTORIES))
    if unknown:
        raise HandoffManifestError(
            "unknown generator label(s): " + ", ".join(unknown)
        )

    paths = set(FROZEN_INPUTS)
    for relative_path in FROZEN_INPUTS:
        _checked_payload_path(root, relative_path, expected_kind="file")
    for generator in sorted(set(generators)):
        paths.update(
            _directory_files(root, GENERATOR_RESULT_DIRECTORIES[generator])
        )
    return tuple(sorted(paths))


def _hash_regular_file(root: Path, relative_path: str) -> ManifestRecord:
    path, before = _checked_payload_path(root, relative_path, expected_kind="file")
    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise HandoffManifestError(f"cannot open payload file {relative_path}: {exc}") from exc

    digest = hashlib.sha256()
    try:
        after = os.fstat(descriptor)
        if not stat.S_ISREG(after.st_mode):
            raise HandoffManifestError(
                f"payload changed to a non-regular file: {relative_path}"
            )
        if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
            raise HandoffManifestError(
                f"payload changed while it was being opened: {relative_path}"
            )
        while True:
            block = os.read(descriptor, _READ_CHUNK_BYTES)
            if not block:
                break
            digest.update(block)
    finally:
        os.close(descriptor)
    return ManifestRecord(digest.hexdigest(), int(after.st_size), relative_path)


def build_records(
    root: Path, generators: tuple[str, ...] = ()
) -> tuple[ManifestRecord, ...]:
    return tuple(
        _hash_regular_file(root, relative_path)
        for relative_path in inventory_paths(root, generators)
    )


def render_manifest(records: tuple[ManifestRecord, ...]) -> str:
    ordered = sorted(records, key=lambda record: record.relative_path)
    return "".join(record.line() + "\n" for record in ordered)


def parse_manifest(text: str) -> tuple[ManifestRecord, ...]:
    records: list[ManifestRecord] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        match = _RECORD_PATTERN.fullmatch(line)
        if match is None:
            raise HandoffManifestError(
                f"invalid manifest record on line {line_number}"
            )
        digest, byte_count, relative_path = match.groups()
        _relative_parts(relative_path)
        records.append(ManifestRecord(digest, int(byte_count), relative_path))

    paths = [record.relative_path for record in records]
    if paths != sorted(paths):
        raise HandoffManifestError("manifest records are not sorted by path")
    if len(paths) != len(set(paths)):
        raise HandoffManifestError("manifest contains duplicate paths")
    return tuple(records)


def _read_text_without_following_symlinks(path: Path) -> str:
    try:
        before = os.lstat(path)
    except OSError as exc:
        raise HandoffManifestError(f"cannot inspect manifest {path}: {exc}") from exc
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise HandoffManifestError(f"manifest is not a regular file: {path}")

    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise HandoffManifestError(f"cannot open manifest {path}: {exc}") from exc
    chunks: list[bytes] = []
    try:
        after = os.fstat(descriptor)
        if not stat.S_ISREG(after.st_mode) or (
            before.st_dev,
            before.st_ino,
        ) != (after.st_dev, after.st_ino):
            raise HandoffManifestError(f"manifest changed while opening: {path}")
        while True:
            block = os.read(descriptor, _READ_CHUNK_BYTES)
            if not block:
                break
            chunks.append(block)
    finally:
        os.close(descriptor)
    try:
        return b"".join(chunks).decode("utf-8")
    except UnicodeDecodeError as exc:
        raise HandoffManifestError(f"manifest is not UTF-8 text: {path}") from exc


def verify_records(
    root: Path,
    records: tuple[ManifestRecord, ...],
    generators: tuple[str, ...] = (),
) -> tuple[int, int]:
    expected_paths = inventory_paths(root, generators)
    recorded_paths = tuple(record.relative_path for record in records)
    if recorded_paths != expected_paths:
        missing = sorted(set(expected_paths) - set(recorded_paths))
        extra = sorted(set(recorded_paths) - set(expected_paths))
        details = []
        if missing:
            details.append("missing=" + ",".join(missing))
        if extra:
            details.append("extra=" + ",".join(extra))
        raise HandoffManifestError(
            "manifest inventory does not match the selected payload"
            + (": " + "; ".join(details) if details else "")
        )

    total_bytes = 0
    for recorded in records:
        actual = _hash_regular_file(root, recorded.relative_path)
        if actual.byte_count != recorded.byte_count:
            raise HandoffManifestError(
                f"byte count changed for {recorded.relative_path}: "
                f"expected {recorded.byte_count}, found {actual.byte_count}"
            )
        if actual.sha256 != recorded.sha256:
            raise HandoffManifestError(
                f"SHA-256 changed for {recorded.relative_path}"
            )
        total_bytes += actual.byte_count
    return len(records), total_bytes


def _add_generator_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--generator",
        action="append",
        default=[],
        choices=tuple(sorted(GENERATOR_RESULT_DIRECTORIES)),
        help="include one generator result directory; repeat for multiple generators",
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    emit = subparsers.add_parser("emit", help="write manifest records to stdout")
    _add_generator_arguments(emit)

    verify = subparsers.add_parser("verify", help="verify an existing manifest")
    verify.add_argument("manifest", type=Path)
    _add_generator_arguments(verify)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    root = repository_root()
    generators = tuple(args.generator)
    try:
        if args.command == "emit":
            sys.stdout.write(render_manifest(build_records(root, generators)))
            return 0
        text = _read_text_without_following_symlinks(args.manifest)
        count, total_bytes = verify_records(
            root, parse_manifest(text), generators
        )
    except HandoffManifestError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    labels = ",".join(sorted(set(generators))) if generators else "inputs-only"
    print(
        f"PASS records={count} bytes={total_bytes} generators={labels}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
