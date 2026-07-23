"""Index the finite-overlap K-generator maps without claiming K coherence.

The four generator tree-contraction producers emit support-local rectangular
maps.  Their NPZ keys are ordered as ``s00``, ``s01``, ..., but the files do
not currently carry the corresponding source and target support coordinates.
This adapter supplies that missing identity layer by rebuilding the generator
overlap under the same bound generator context and attaching explicit support
tuples to every stored slot.

It also emits the exact Cayley-cell registry for ``K=(Z/3)^4``:

* 324 directed generator edges;
* 108 distinct directed order-three cycles; and
* 486 generator-commutation squares.

These registries are not numerical coherence certificates.  The stored tree
maps are rectangular restrictions into generator-dependent finite overlaps,
not invertible transitions in a common basis.  Consequently this module marks
an edge as ready for common-basis rebasing when its generator artifact and
support identities are available, but it never marks an edge matrix, cycle, or
square closed.  A later adapter must construct common bases and perform the
actual matrix products (or emit the required chain homotopies).

The adapter is partial-availability aware.  Re-running it discovers any newly
completed generator manifests; an optional preflight file is recorded as
evidence but is never treated as a prerequisite.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_rank2_k_chain_comparison as k_chain,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_rank2_k_generator_overlap_j1 as overlap_j1,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_rank2_k_generator_tree_contraction as tree_contraction,
)
from experiments.section12_instantiation import (  # noqa: E402
    cm_bridge_hs_rank2_k_generator_tree_contraction_driver as driver,
)


P = 31
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/"
    "cm_bridge_hs_rank2_k_coherence_adapter"
)
COMMITTED_CAYLEY_JOBS = k_chain.OUTPUT_DIRECTORY / (
    "z2_k_cayley_chain_comparison_jobs.csv"
)

TREE_MANIFEST_SCHEMA = (
    "section12-HS-rank2-K-generator-tree-contraction-manifest-v3"
)
TREE_CERTIFICATE_SCHEMA = (
    "section12-HS-rank2-K-generator-tree-contraction-v3"
)
ADAPTER_CERTIFICATE_SCHEMA = "section12-HS-rank2-K-coherence-adapter-v1"
ADAPTER_MANIFEST_SCHEMA = (
    "section12-HS-rank2-K-coherence-adapter-manifest-v1"
)
PAIR_LEDGER_HASH_ENCODING = "signed-little-endian-int64-C-order"
METADATA_REPAIR_SCHEMA = (
    "section12-HS-rank2-K-generator-metadata-repair-v1"
)

CORE_MATRIX_NAMES = ("J1", "J2", "J3_G3", "L0", "L1", "a")
MATRIX_COMPONENTS = ("value", "tangent")
NORMALIZATION_HASH_KEYS = {
    "target_lambda_normalization_pairs": (
        "target_lambda_normalization_pairs_sha256"
    ),
    "residue_transport_pairs": "residue_transport_pairs_sha256",
    "raw_counit_weight_pairs": "raw_counit_weight_pairs_sha256",
    "orientation_coboundary_pairs": (
        "orientation_coboundary_pairs_sha256"
    ),
}
GENERATOR_LABELS = tuple(driver.GENERATOR_EXPONENTS)
GENERATOR_BY_EXPONENTS = {
    tuple(exponents): label
    for label, exponents in driver.GENERATOR_EXPONENTS.items()
}
VERTEX_EXPONENTS = tuple(itertools.product(range(3), repeat=4))
VERTEX_INDEX_BY_EXPONENTS = {
    tuple(exponents): index
    for index, exponents in enumerate(VERTEX_EXPONENTS)
}


class HSRank2KCoherenceAdapterError(RuntimeError):
    """An artifact, support identity, or deterministic K registry is invalid."""


class GeneratorArtifactUnavailable(HSRank2KCoherenceAdapterError):
    """A generator has no completed manifest yet."""


class GeneratorArtifactInvalid(HSRank2KCoherenceAdapterError):
    """A completed generator manifest does not verify internally."""


@dataclass(frozen=True)
class GeneratorArtifact:
    """Verified metadata for one lazily loaded generator map artifact."""

    generator: str
    exponents: tuple[int, int, int, int]
    output_directory: Path
    maps_path: Path
    certificate_path: Path
    hermite_path: Path
    support_count: int
    overlap_length: int
    matrix_shapes: Mapping[str, tuple[int, int]]
    source_binding: str
    verified_files: tuple[str, ...]
    normalization_arrays: tuple[str, ...]
    maps_sha256: str


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Hash a possibly multi-gigabyte artifact without reading it at once."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def pair_ledger_sha256(value) -> str:
    """Hash pair ledgers with the producer's canonical integer encoding."""

    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(value, dtype="<i8")).tobytes()
    ).hexdigest()


def json_compact(value) -> str:
    return json.dumps(value, separators=(",", ":"))


def write_json_atomic(path: Path, value) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_csv_atomic(path: Path, rows, fieldnames: Sequence[str]) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def load_json(path: Path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise GeneratorArtifactInvalid(
            f"could not read JSON artifact {path}: {error}"
        ) from error


def add_exponents(left, right):
    return tuple(
        (int(left[index]) + int(right[index])) % 3 for index in range(4)
    )


def generator_step(generator: int) -> tuple[int, int, int, int]:
    return tuple(int(index == generator) for index in range(4))


def _expected_matrix_shapes(generator: str):
    dimensions = driver.predicted_dimensions(generator)
    return {
        "J1": (
            dimensions["source_F1"],
            dimensions["target_F1"],
        ),
        "J2": (
            dimensions["source_F2"],
            dimensions["target_F2"],
        ),
        "J3_G3": (
            56 * dimensions["source_F0"],
            56 * dimensions["target_F0"],
        ),
        "L0": (
            dimensions["source_L0_rows"],
            dimensions["target_L0_columns"],
        ),
        "L1": (
            dimensions["source_F2"],
            dimensions["target_F2"],
        ),
        "a": (1, dimensions["target_F1"]),
    }


def _verify_manifest_files(output: Path, manifest) -> tuple[str, ...]:
    files = manifest.get("files")
    if not isinstance(files, dict):
        raise GeneratorArtifactInvalid("the generator manifest has no file ledger")
    required = {
        "tree_contracted_support_maps.npz",
        "overlap_adapted_hermite.npz",
        "tree_contraction_certificate.json",
    }
    if not required.issubset(files):
        missing = sorted(required - set(files))
        raise GeneratorArtifactInvalid(
            f"the generator manifest omits required files {missing}"
        )

    verified = []
    for name, record in sorted(files.items()):
        path = output / name
        if not path.is_file():
            raise GeneratorArtifactInvalid(
                f"manifest-listed generator file is absent: {path}"
            )
        if not isinstance(record, dict):
            raise GeneratorArtifactInvalid(
                f"the manifest record for {name} is not a mapping"
            )
        try:
            expected_bytes = int(record.get("bytes", -1))
        except (TypeError, ValueError) as error:
            raise GeneratorArtifactInvalid(
                f"the manifest byte count for {name} is invalid"
            ) from error
        if path.stat().st_size != expected_bytes:
            raise GeneratorArtifactInvalid(
                f"the byte count changed for {path}"
            )
        expected_sha256 = record.get("sha256")
        if not expected_sha256 or sha256(path) != expected_sha256:
            raise GeneratorArtifactInvalid(f"the SHA-256 changed for {path}")
        verified.append(name)
    return tuple(verified)


def _require_npz_scalar_vector(archive, name, expected) -> None:
    if name not in archive.files:
        raise GeneratorArtifactInvalid(f"the map NPZ omits {name}")
    actual = tuple(map(int, np.asarray(archive[name]).reshape(-1)))
    if actual != tuple(map(int, expected)):
        raise GeneratorArtifactInvalid(
            f"the map NPZ value for {name} is {actual}, expected {tuple(expected)}"
        )


def _canonical_normalization_certificate(
    output: Path,
    manifest_path: Path,
    manifest: Mapping,
    certificate: Mapping,
    generator: str,
) -> Mapping:
    """Return canonical pair-ledger hashes, using the e1 repair sidecar.

    The first full e1 artifact recorded hashes of the stored uint8 arrays.
    Its non-destructive repair sidecar binds the same manifest and arrays to
    the canonical signed little-endian int64 encoding used by later producer
    certificates.  No numeric artifact is rewritten here.
    """

    normalization = certificate.get("generator_normalization", {})
    if not isinstance(normalization, dict):
        raise GeneratorArtifactInvalid(
            "the generator normalization certificate is not a mapping"
        )
    encoding = normalization.get("pair_ledger_hash_encoding")
    if encoding == PAIR_LEDGER_HASH_ENCODING:
        return normalization
    if encoding is not None:
        raise GeneratorArtifactInvalid(
            f"unsupported pair-ledger hash encoding {encoding!r}"
        )

    repair_path = output / "metadata_repair_certificate.json"
    if not repair_path.is_file():
        return normalization
    repair = load_json(repair_path)
    if (
        repair.get("schema") != METADATA_REPAIR_SCHEMA
        or repair.get("status") != "PASS"
        or repair.get("generator") != generator
        or tuple(repair.get("K_generator_exponents", ()))
        != tuple(driver.GENERATOR_EXPONENTS[generator])
    ):
        raise GeneratorArtifactInvalid(
            "the generator metadata-repair sidecar identity is invalid"
        )
    original = repair.get("original_artifact", {})
    if (
        original.get("manifest", {}).get("sha256")
        != sha256(manifest_path)
        or original.get("verified_files") != manifest.get("files")
    ):
        raise GeneratorArtifactInvalid(
            "the metadata-repair sidecar does not bind the current manifest"
        )

    ledgers = repair.get("pair_ledgers")
    if not isinstance(ledgers, dict):
        raise GeneratorArtifactInvalid(
            "the metadata-repair sidecar has no pair-ledger mapping"
        )
    canonical = {
        "pair_ledger_hash_encoding": PAIR_LEDGER_HASH_ENCODING,
    }
    for name, certificate_key in NORMALIZATION_HASH_KEYS.items():
        record = ledgers.get(name)
        if (
            not isinstance(record, dict)
            or record.get("canonical_encoding")
            != PAIR_LEDGER_HASH_ENCODING
            or not record.get("canonical_sha256")
        ):
            raise GeneratorArtifactInvalid(
                f"the metadata-repair sidecar has no canonical hash for {name}"
            )
        canonical[certificate_key] = record["canonical_sha256"]
    return canonical


def _verify_map_npz(
    maps_path: Path,
    generator: str,
    support_count: int,
    normalization_certificate: Mapping | None = None,
) -> tuple[tuple[str, ...], Mapping[str, tuple[int, int]]]:
    shapes = _expected_matrix_shapes(generator)
    dimensions = driver.predicted_dimensions(generator)
    try:
        with np.load(maps_path, allow_pickle=False) as archive:
            _require_npz_scalar_vector(archive, "field", (P,))
            _require_npz_scalar_vector(
                archive,
                "generator_exponents",
                driver.GENERATOR_EXPONENTS[generator],
            )
            _require_npz_scalar_vector(archive, "J1_shape", shapes["J1"])
            _require_npz_scalar_vector(archive, "J2_shape", shapes["J2"])
            _require_npz_scalar_vector(
                archive,
                "J3_shape",
                (
                    dimensions["source_F3"],
                    dimensions["target_F3"],
                ),
            )
            _require_npz_scalar_vector(archive, "L0_shape", shapes["L0"])

            keys = set(archive.files)
            required_keys = set()
            for support in range(support_count):
                prefix = f"s{support:02d}"
                for name in CORE_MATRIX_NAMES:
                    for component in MATRIX_COMPONENTS:
                        required_keys.add(f"{prefix}_{name}_{component}")
                for equation in range(6):
                    for component in MATRIX_COMPONENTS:
                        required_keys.add(
                            f"{prefix}_J3_Omega_i1_e{equation}_{component}"
                        )
                required_keys.add(f"{prefix}_J3_B_value")
                required_keys.add(f"{prefix}_J3_B_tangent")
            missing = sorted(required_keys - keys)
            if missing:
                preview = missing[:5]
                raise GeneratorArtifactInvalid(
                    "the map NPZ omits support-local records "
                    f"{preview}{' ...' if len(missing) > len(preview) else ''}"
                )

            normalization_arrays = []
            hash_encoding = (normalization_certificate or {}).get(
                "pair_ledger_hash_encoding"
            )
            if hash_encoding not in (None, PAIR_LEDGER_HASH_ENCODING):
                raise GeneratorArtifactInvalid(
                    "the pair-ledger hash encoding is unsupported"
                )
            for name, certificate_key in NORMALIZATION_HASH_KEYS.items():
                if name not in keys:
                    continue
                value = np.asarray(archive[name])
                if value.shape != (support_count, 2):
                    raise GeneratorArtifactInvalid(
                        f"{name} has shape {value.shape}, expected "
                        f"({support_count}, 2)"
                    )
                expected_hash = (normalization_certificate or {}).get(
                    certificate_key
                )
                if hash_encoding == PAIR_LEDGER_HASH_ENCODING:
                    if not expected_hash:
                        raise GeneratorArtifactInvalid(
                            f"the canonical certificate omits {certificate_key}"
                        )
                    actual_hash = pair_ledger_sha256(value)
                else:
                    actual_hash = hashlib.sha256(
                        np.ascontiguousarray(value).tobytes()
                    ).hexdigest()
                if expected_hash and actual_hash != expected_hash:
                    raise GeneratorArtifactInvalid(
                        f"{name} does not replay its certificate hash"
                    )
                normalization_arrays.append(name)
    except GeneratorArtifactInvalid:
        raise
    except (OSError, ValueError, KeyError) as error:
        raise GeneratorArtifactInvalid(
            f"could not inspect map NPZ {maps_path}: {error}"
        ) from error
    return tuple(normalization_arrays), shapes


def load_generator_artifact(
    generator: str,
    output_directory: Path | str | None = None,
) -> GeneratorArtifact:
    """Verify one completed generator artifact and return lazy map metadata.

    File hashes, certificate metadata, global NPZ shape records, and the
    support-local key census are checked.  The sparse matrices remain lazy
    because one generator NPZ expands to multiple gigabytes.  Use
    :func:`load_support_matrix_pair` or :func:`iter_support_matrix_pairs` to
    reconstruct CSR matrices on demand.
    """

    if generator not in driver.GENERATOR_EXPONENTS:
        raise GeneratorArtifactInvalid(f"unknown generator {generator!r}")
    output = (
        Path(output_directory).expanduser().resolve()
        if output_directory is not None
        else Path(driver.OUTPUT_DIRECTORIES[generator]).resolve()
    )
    manifest_path = output / "manifest.json"
    if not manifest_path.is_file():
        raise GeneratorArtifactUnavailable(
            f"generator {generator} has no completed manifest in {output}"
        )

    manifest = load_json(manifest_path)
    if manifest.get("schema") != TREE_MANIFEST_SCHEMA:
        raise GeneratorArtifactInvalid(
            f"unsupported generator manifest schema {manifest.get('schema')!r}"
        )
    if not isinstance(manifest.get("generator"), dict):
        raise GeneratorArtifactInvalid(
            "the generator manifest has no producer source record"
        )
    verified_files = _verify_manifest_files(output, manifest)

    certificate_path = output / "tree_contraction_certificate.json"
    certificate = load_json(certificate_path)
    if certificate.get("schema") != TREE_CERTIFICATE_SCHEMA:
        raise GeneratorArtifactInvalid(
            "the generator tree-contraction certificate schema changed"
        )
    if certificate.get("field") != "F_31":
        raise GeneratorArtifactInvalid("the generator certificate field changed")
    expected_exponents = tuple(driver.GENERATOR_EXPONENTS[generator])
    if tuple(certificate.get("K_generator_exponents", ())) != expected_exponents:
        raise GeneratorArtifactInvalid(
            "the certificate generator exponents do not match its directory"
        )

    counts = certificate.get("counts", {})
    dimensions = driver.predicted_dimensions(generator)
    expected_counts = {
        "overlap_supports": dimensions["overlap_factors"],
        "overlap_length": dimensions["overlap_length"],
        "target_length": dimensions["target_F0"],
        "J1_shape": list(_expected_matrix_shapes(generator)["J1"]),
        "J2_shape": list(_expected_matrix_shapes(generator)["J2"]),
        "J3_shape": [
            dimensions["source_F3"],
            dimensions["target_F3"],
        ],
        "L0_shape": list(_expected_matrix_shapes(generator)["L0"]),
    }
    for name, expected in expected_counts.items():
        if counts.get(name) != expected:
            raise GeneratorArtifactInvalid(
                f"certificate count {name} is {counts.get(name)!r}, "
                f"expected {expected!r}"
            )
    checks = certificate.get("exact_checks", {})
    if checks.get("all_tree_chain_residuals_zero") is not True:
        raise GeneratorArtifactInvalid(
            "the generator certificate does not certify zero tree residuals"
        )
    if checks.get("rectangular_map_is_not_asserted_invertible") is not True:
        raise GeneratorArtifactInvalid(
            "the generator certificate lost its rectangular-map boundary"
        )

    maps_path = output / "tree_contracted_support_maps.npz"
    normalization_certificate = _canonical_normalization_certificate(
        output,
        manifest_path,
        manifest,
        certificate,
        generator,
    )
    normalization_arrays, matrix_shapes = _verify_map_npz(
        maps_path,
        generator,
        dimensions["overlap_factors"],
        normalization_certificate,
    )
    recorded_source_hash = manifest.get("generator", {}).get("sha256")
    current_source_hash = sha256(Path(tree_contraction.__file__).resolve())
    source_binding = (
        "CURRENT_SOURCE_BOUND"
        if recorded_source_hash == current_source_hash
        else "LEGACY_CURRENT_SOURCE_MISMATCH"
    )
    return GeneratorArtifact(
        generator=generator,
        exponents=expected_exponents,
        output_directory=output,
        maps_path=maps_path,
        certificate_path=certificate_path,
        hermite_path=output / "overlap_adapted_hermite.npz",
        support_count=dimensions["overlap_factors"],
        overlap_length=dimensions["overlap_length"],
        matrix_shapes=matrix_shapes,
        source_binding=source_binding,
        verified_files=verified_files,
        normalization_arrays=normalization_arrays,
        maps_sha256=manifest["files"][maps_path.name]["sha256"],
    )


def csr_from_records(records, shape) -> sparse.csr_matrix:
    """Reconstruct one exact F_31 CSR matrix from ``(row,column,value)`` rows."""

    shape = tuple(map(int, shape))
    if len(shape) != 2 or min(shape) < 0:
        raise HSRank2KCoherenceAdapterError(
            f"invalid sparse matrix shape {shape}"
        )
    records = np.asarray(records)
    if records.ndim != 2 or records.shape[1] != 3:
        raise HSRank2KCoherenceAdapterError(
            "a CSR record array must have shape (n,3)"
        )
    if not len(records):
        return sparse.csr_matrix(shape, dtype=np.int64)

    rows = np.asarray(records[:, 0], dtype=np.int64)
    columns = np.asarray(records[:, 1], dtype=np.int64)
    values = np.asarray(records[:, 2], dtype=np.int64) % P
    if (
        np.any(rows < 0)
        or np.any(columns < 0)
        or np.any(rows >= shape[0])
        or np.any(columns >= shape[1])
    ):
        raise HSRank2KCoherenceAdapterError(
            "a sparse record index lies outside its declared shape"
        )
    matrix = sparse.csr_matrix(
        (values, (rows, columns)),
        shape=shape,
        dtype=np.int64,
    )
    matrix.sum_duplicates()
    matrix.data %= P
    matrix.eliminate_zeros()
    matrix.sort_indices()
    return matrix


def load_support_matrix_pair(
    artifact: GeneratorArtifact,
    support: int,
    matrix_name: str,
):
    """Lazily reconstruct one support-local value/tangent CSR pair."""

    support = int(support)
    if support < 0 or support >= artifact.support_count:
        raise HSRank2KCoherenceAdapterError(
            f"support {support} is outside 0..{artifact.support_count - 1}"
        )
    if matrix_name not in CORE_MATRIX_NAMES:
        raise HSRank2KCoherenceAdapterError(
            f"unsupported support matrix {matrix_name!r}"
        )
    shape = artifact.matrix_shapes[matrix_name]
    prefix = f"s{support:02d}_{matrix_name}"
    try:
        with np.load(artifact.maps_path, allow_pickle=False) as archive:
            value = csr_from_records(archive[f"{prefix}_value"], shape)
            tangent = csr_from_records(archive[f"{prefix}_tangent"], shape)
    except HSRank2KCoherenceAdapterError:
        raise
    except (OSError, ValueError, KeyError) as error:
        raise HSRank2KCoherenceAdapterError(
            f"could not reconstruct {prefix}: {error}"
        ) from error
    return value, tangent


def iter_support_matrix_pairs(
    artifact: GeneratorArtifact,
    matrix_names: Iterable[str] = CORE_MATRIX_NAMES,
) -> Iterator[tuple[int, str, tuple[sparse.csr_matrix, sparse.csr_matrix]]]:
    """Yield reconstructed pairs without retaining a whole multi-GB NPZ."""

    matrix_names = tuple(matrix_names)
    unknown = sorted(set(matrix_names) - set(CORE_MATRIX_NAMES))
    if unknown:
        raise HSRank2KCoherenceAdapterError(
            f"unsupported support matrices {unknown}"
        )
    try:
        with np.load(artifact.maps_path, allow_pickle=False) as archive:
            for support in range(artifact.support_count):
                for matrix_name in matrix_names:
                    prefix = f"s{support:02d}_{matrix_name}"
                    shape = artifact.matrix_shapes[matrix_name]
                    yield support, matrix_name, (
                        csr_from_records(archive[f"{prefix}_value"], shape),
                        csr_from_records(archive[f"{prefix}_tangent"], shape),
                    )
    except HSRank2KCoherenceAdapterError:
        raise
    except (OSError, ValueError, KeyError) as error:
        raise HSRank2KCoherenceAdapterError(
            f"could not iterate support matrices: {error}"
        ) from error


def _decode_support_tuple(value, field: str):
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except ValueError as error:
            raise HSRank2KCoherenceAdapterError(
                f"could not decode {field}: {error}"
            ) from error
    try:
        result = tuple(map(int, value))
    except (TypeError, ValueError) as error:
        raise HSRank2KCoherenceAdapterError(
            f"could not decode {field}"
        ) from error
    if len(result) != 4:
        raise HSRank2KCoherenceAdapterError(
            f"{field} is not a four-factor support tuple"
        )
    return result


def regenerate_support_bindings(
    generator: str,
    support_count: int | None = None,
    output_directory: Path | str | None = None,
):
    """Rebuild and tuple-bind the identities behind ``s00...sNN``.

    The artifact slot is necessarily attached in the producer's deterministic
    point-row order because the legacy NPZ has no coordinate arrays.  Every
    downstream lookup returned by this function is keyed by the explicit
    ``(source_support, target_support)`` tuple pair, never by a support index.
    """

    if generator not in driver.GENERATOR_EXPONENTS:
        raise HSRank2KCoherenceAdapterError(
            f"unknown generator {generator!r}"
        )
    output = (
        Path(output_directory).expanduser().resolve()
        if output_directory is not None
        else Path(driver.OUTPUT_DIRECTORIES[generator]).resolve()
    )
    with driver.bind_generator(generator, output):
        graph_data = overlap_j1.build_data()
    rows = graph_data.get("point_rows", ())
    expected = (
        driver.predicted_dimensions(generator)["overlap_factors"]
        if support_count is None
        else int(support_count)
    )
    if len(rows) != expected:
        raise HSRank2KCoherenceAdapterError(
            f"generator {generator} regenerated {len(rows)} supports, "
            f"expected {expected}"
        )

    bindings = []
    lookup = {}
    for slot, row in enumerate(rows):
        source = _decode_support_tuple(
            row.get("source_support"), "source_support"
        )
        target = _decode_support_tuple(
            row.get("target_support"), "target_support"
        )
        pair_key = (source, target)
        if pair_key in lookup:
            raise HSRank2KCoherenceAdapterError(
                f"generator {generator} repeated support pair {pair_key}"
            )
        binding = {
            "generator": generator,
            "generator_exponents": tuple(driver.GENERATOR_EXPONENTS[generator]),
            "artifact_support_slot": f"s{slot:02d}",
            "source_support": source,
            "target_support": target,
            "support_pair_key": pair_key,
        }
        lookup[pair_key] = binding
        bindings.append(binding)
    if len(lookup) != expected:
        raise HSRank2KCoherenceAdapterError(
            f"generator {generator} tuple lookup changed cardinality"
        )
    return bindings


def support_binding_lookup(bindings):
    """Return the explicit tuple-keyed support join used downstream."""

    lookup = {}
    for binding in bindings:
        source = tuple(binding["source_support"])
        target = tuple(binding["target_support"])
        key = (source, target)
        if key in lookup:
            raise HSRank2KCoherenceAdapterError(
                f"duplicate support tuple join key {key}"
            )
        lookup[key] = binding
    return lookup


def support_bindings_sha256(bindings) -> str:
    payload = [
        [
            list(map(int, binding["source_support"])),
            list(map(int, binding["target_support"])),
        ]
        for binding in bindings
    ]
    return hashlib.sha256(json_compact(payload).encode("utf-8")).hexdigest()


def enumerate_typed_edges(generator_states: Mapping | None = None):
    """Enumerate the exact 81-by-4 directed Cayley edges."""

    states = generator_states or {}
    states_by_exponents = {
        tuple(driver.GENERATOR_EXPONENTS[label]): states.get(label, {})
        for label in GENERATOR_LABELS
    }
    rows = []
    for source_index, source in enumerate(VERTEX_EXPONENTS):
        for generator in range(4):
            step = generator_step(generator)
            label = GENERATOR_BY_EXPONENTS[step]
            target = add_exponents(source, step)
            target_index = VERTEX_INDEX_BY_EXPONENTS[target]
            state = states_by_exponents.get(step, {})
            ready = bool(state.get("ready_for_common_basis_rebase", False))
            rows.append(
                {
                    "edge_job_id": len(rows),
                    "source_K_index": source_index,
                    "target_K_index": target_index,
                    "source_K_exponents": source,
                    "target_K_exponents": target,
                    "generator_index": generator,
                    "generator": label,
                    "generator_exponents": step,
                    "generator_artifact_loaded": bool(
                        state.get("artifact_loaded", False)
                    ),
                    "support_ids_bound": bool(
                        state.get("support_ids_bound", False)
                    ),
                    "support_binding_count": int(
                        state.get("support_binding_count", 0)
                    ),
                    "support_bindings_sha256": state.get(
                        "support_bindings_sha256"
                    ),
                    "ready_for_common_basis_rebase": ready,
                    "common_basis_rebased": False,
                    "numeric_edge_matrix_emitted": False,
                    "edge_chain_map_certified": False,
                    "status": (
                        "READY_FOR_COMMON_BASIS_REBASE"
                        if ready
                        else "WAITING_FOR_GENERATOR_ARTIFACT_OR_SUPPORT_IDS"
                    ),
                }
            )
    if len(rows) != 324 or len(
        {
            (
                row["source_K_exponents"],
                row["generator_exponents"],
            )
            for row in rows
        }
    ) != 324:
        raise HSRank2KCoherenceAdapterError(
            "the typed Cayley-edge registry changed"
        )
    return rows


def enumerate_order_three_cycles(edges):
    """Enumerate one directed three-cycle per generator coset."""

    edge_lookup = {
        (tuple(row["source_K_exponents"]), int(row["generator_index"])): row
        for row in edges
    }
    cycles = []
    for generator in range(4):
        step = generator_step(generator)
        for base in VERTEX_EXPONENTS:
            if base[generator] != 0:
                continue
            vertices = [base]
            edge_rows = []
            current = base
            for _ in range(3):
                edge = edge_lookup[(current, generator)]
                edge_rows.append(edge)
                current = add_exponents(current, step)
                vertices.append(current)
            if current != base:
                raise HSRank2KCoherenceAdapterError(
                    "a generator order-three path did not close combinatorially"
                )
            inputs_ready = all(
                row["ready_for_common_basis_rebase"] for row in edge_rows
            )
            cycles.append(
                {
                    "cycle_id": len(cycles),
                    "generator_index": generator,
                    "generator": GENERATOR_BY_EXPONENTS[step],
                    "generator_exponents": step,
                    "base_K_index": VERTEX_INDEX_BY_EXPONENTS[base],
                    "base_K_exponents": base,
                    "vertex_K_indices": tuple(
                        VERTEX_INDEX_BY_EXPONENTS[vertex]
                        for vertex in vertices
                    ),
                    "edge_job_ids": tuple(
                        int(row["edge_job_id"]) for row in edge_rows
                    ),
                    "generator_templates_available": inputs_ready,
                    "common_basis_composition_emitted": False,
                    "strict_equality_computed": False,
                    "chain_homotopy_emitted": False,
                    "coherence_closed": False,
                    "status": (
                        "OPEN_COMMON_BASIS_COMPOSITION_NOT_EMITTED"
                        if inputs_ready
                        else "WAITING_FOR_GENERATOR_ARTIFACT_OR_SUPPORT_IDS"
                    ),
                }
            )
    if len(cycles) != 108:
        raise HSRank2KCoherenceAdapterError(
            "the order-three-cycle census changed"
        )
    return cycles


def enumerate_commutation_squares(edges):
    """Enumerate all 81 squares for each unordered generator pair."""

    edge_lookup = {
        (tuple(row["source_K_exponents"]), int(row["generator_index"])): row
        for row in edges
    }
    squares = []
    for base in VERTEX_EXPONENTS:
        for left_generator, right_generator in itertools.combinations(
            range(4), 2
        ):
            left_step = generator_step(left_generator)
            right_step = generator_step(right_generator)
            after_left = add_exponents(base, left_step)
            after_right = add_exponents(base, right_step)
            endpoint_left_right = add_exponents(after_left, right_step)
            endpoint_right_left = add_exponents(after_right, left_step)
            if endpoint_left_right != endpoint_right_left:
                raise HSRank2KCoherenceAdapterError(
                    "a generator square does not share an endpoint"
                )
            path_left_right = (
                edge_lookup[(base, left_generator)],
                edge_lookup[(after_left, right_generator)],
            )
            path_right_left = (
                edge_lookup[(base, right_generator)],
                edge_lookup[(after_right, left_generator)],
            )
            inputs_ready = all(
                row["ready_for_common_basis_rebase"]
                for row in path_left_right + path_right_left
            )
            squares.append(
                {
                    "square_id": len(squares),
                    "base_K_index": VERTEX_INDEX_BY_EXPONENTS[base],
                    "base_K_exponents": base,
                    "left_generator_index": left_generator,
                    "right_generator_index": right_generator,
                    "left_generator": GENERATOR_BY_EXPONENTS[left_step],
                    "right_generator": GENERATOR_BY_EXPONENTS[right_step],
                    "endpoint_K_index": VERTEX_INDEX_BY_EXPONENTS[
                        endpoint_left_right
                    ],
                    "endpoint_K_exponents": endpoint_left_right,
                    "left_then_right_edge_job_ids": tuple(
                        int(row["edge_job_id"]) for row in path_left_right
                    ),
                    "right_then_left_edge_job_ids": tuple(
                        int(row["edge_job_id"]) for row in path_right_left
                    ),
                    "generator_templates_available": inputs_ready,
                    "common_basis_compositions_emitted": False,
                    "strict_equality_computed": False,
                    "chain_homotopy_emitted": False,
                    "coherence_closed": False,
                    "status": (
                        "OPEN_COMMON_BASIS_COMPOSITIONS_NOT_EMITTED"
                        if inputs_ready
                        else "WAITING_FOR_GENERATOR_ARTIFACTS_OR_SUPPORT_IDS"
                    ),
                }
            )
    if len(squares) != 486:
        raise HSRank2KCoherenceAdapterError(
            "the generator-commutation-square census changed"
        )
    return squares


def compare_committed_cayley_registry(
    edges,
    path: Path | str = COMMITTED_CAYLEY_JOBS,
):
    """Compare the emitted registry to the committed job table by tuples."""

    path = Path(path)
    if not path.is_file():
        return {
            "checked": False,
            "matches": None,
            "path": str(path),
            "reason": "committed Cayley job table is absent",
        }
    with path.open(newline="", encoding="utf-8") as handle:
        committed = list(csv.DictReader(handle))
    if len(committed) != len(edges):
        return {
            "checked": True,
            "matches": False,
            "path": str(path),
            "reason": "edge row count differs",
        }
    for emitted, stored in zip(edges, committed):
        try:
            stored_source = tuple(
                map(int, json.loads(stored["source_K_exponents"]))
            )
            stored_target = tuple(
                map(int, json.loads(stored["target_K_exponents"]))
            )
            stored_generator = int(stored["generator"])
            stored_job = int(stored["job_id"])
        except (KeyError, TypeError, ValueError):
            return {
                "checked": True,
                "matches": False,
                "path": str(path),
                "reason": "a committed row could not be decoded",
            }
        if (
            stored_job != emitted["edge_job_id"]
            or stored_source != emitted["source_K_exponents"]
            or stored_target != emitted["target_K_exponents"]
            or generator_step(stored_generator)
            != emitted["generator_exponents"]
        ):
            return {
                "checked": True,
                "matches": False,
                "path": str(path),
                "reason": f"tuple mismatch at edge {emitted['edge_job_id']}",
            }
    return {
        "checked": True,
        "matches": True,
        "path": str(path),
        "reason": None,
    }


def _inspect_generators(
    generator_directories: Mapping[str, Path | str] | None = None,
):
    directories = dict(driver.OUTPUT_DIRECTORIES)
    if generator_directories:
        directories.update(generator_directories)

    states = {}
    artifacts = {}
    all_bindings = []
    for generator in GENERATOR_LABELS:
        output = Path(directories[generator]).expanduser().resolve()
        preflight_available = driver.preflight_path(
            generator, output
        ).is_file()
        state = {
            "generator": generator,
            "generator_exponents": tuple(
                driver.GENERATOR_EXPONENTS[generator]
            ),
            "output_directory": str(output),
            "preflight_available": preflight_available,
            "preflight_is_prerequisite": False,
            "artifact_status": (
                "PREFLIGHT_ONLY"
                if preflight_available
                else "ARTIFACT_ABSENT"
            ),
            "artifact_loaded": False,
            "internal_files_verified": False,
            "source_binding": None,
            "support_ids_bound": False,
            "support_binding_count": 0,
            "support_bindings_sha256": None,
            "ready_for_common_basis_rebase": False,
            "common_basis_rebased": False,
            "numeric_edge_matrix_emitted": False,
            "validation_error": None,
        }
        try:
            artifact = load_generator_artifact(generator, output)
        except GeneratorArtifactUnavailable:
            states[generator] = state
            continue
        except GeneratorArtifactInvalid as error:
            state["artifact_status"] = "INVALID_AVAILABLE_ARTIFACT"
            state["validation_error"] = str(error)
            states[generator] = state
            continue

        artifacts[generator] = artifact
        state.update(
            {
                "artifact_status": "LOADED_VERIFIED",
                "artifact_loaded": True,
                "internal_files_verified": True,
                "source_binding": artifact.source_binding,
                "verified_file_count": len(artifact.verified_files),
                "maps_sha256": artifact.maps_sha256,
                "normalization_arrays": artifact.normalization_arrays,
                "support_count_in_artifact": artifact.support_count,
                "matrix_shapes": dict(artifact.matrix_shapes),
            }
        )
        try:
            bindings = regenerate_support_bindings(
                generator,
                artifact.support_count,
                output,
            )
            tuple_lookup = support_binding_lookup(bindings)
            if len(tuple_lookup) != artifact.support_count:
                raise HSRank2KCoherenceAdapterError(
                    "the tuple support lookup has the wrong size"
                )
        except (
            HSRank2KCoherenceAdapterError,
            overlap_j1.HSKGeneratorOverlapJ1Error,
            driver.KGeneratorTreeContractionDriverError,
            OSError,
            ValueError,
            KeyError,
        ) as error:
            state["artifact_status"] = "LOADED_SUPPORT_ID_BINDING_FAILED"
            state["validation_error"] = str(error)
            states[generator] = state
            continue

        binding_hash = support_bindings_sha256(bindings)
        state.update(
            {
                "artifact_status": "LOADED_SUPPORT_IDS_BOUND",
                "support_ids_bound": True,
                "support_binding_count": len(bindings),
                "support_bindings_sha256": binding_hash,
                "ready_for_common_basis_rebase": True,
            }
        )
        for binding in bindings:
            binding["support_bindings_sha256"] = binding_hash
        all_bindings.extend(bindings)
        states[generator] = state
    return states, artifacts, all_bindings


def build(
    generator_directories: Mapping[str, Path | str] | None = None,
    committed_cayley_jobs: Path | str = COMMITTED_CAYLEY_JOBS,
):
    """Build a partial-aware artifact census and exact open-cell registry."""

    states, artifacts, support_bindings = _inspect_generators(
        generator_directories
    )
    edges = enumerate_typed_edges(states)
    cycles = enumerate_order_three_cycles(edges)
    squares = enumerate_commutation_squares(edges)
    committed_comparison = compare_committed_cayley_registry(
        edges, committed_cayley_jobs
    )

    loaded = sum(state["artifact_loaded"] for state in states.values())
    support_bound = sum(
        state["support_ids_bound"] for state in states.values()
    )
    edges_ready = sum(
        row["ready_for_common_basis_rebase"] for row in edges
    )
    cycles_ready = sum(
        row["generator_templates_available"] for row in cycles
    )
    squares_ready = sum(
        row["generator_templates_available"] for row in squares
    )
    certificate = {
        "schema": ADAPTER_CERTIFICATE_SCHEMA,
        "field": "F_31",
        "status": (
            "ALL_GENERATOR_TEMPLATES_LOADED;"
            "COMMON_BASIS_NUMERIC_COMPOSITION_OPEN"
            if support_bound == 4
            else "PARTIAL_GENERATOR_TEMPLATES_LOADED;"
            "COMMON_BASIS_NUMERIC_COMPOSITION_OPEN"
        ),
        "scope": (
            "verified generator artifacts, tuple-bound support identities, "
            "and deterministic K edge/cycle/square registries only"
        ),
        "counts": {
            "K_vertices": len(VERTEX_EXPONENTS),
            "standard_generators": 4,
            "generator_artifacts_loaded": loaded,
            "generators_with_support_ids_bound": support_bound,
            "support_identity_bindings": len(support_bindings),
            "typed_edges": len(edges),
            "typed_edges_ready_for_common_basis_rebase": edges_ready,
            "typed_edges_common_basis_rebased": 0,
            "numeric_edge_chain_maps_emitted": 0,
            "distinct_order_three_cycles": len(cycles),
            "cycles_with_generator_templates_available": cycles_ready,
            "order_three_cycles_closed": 0,
            "generator_commutation_squares": len(squares),
            "squares_with_generator_templates_available": squares_ready,
            "generator_commutation_squares_closed": 0,
            "strict_coherence_equalities_computed": 0,
            "coherence_chain_homotopies_emitted": 0,
        },
        "generator_artifacts": states,
        "exact_registry_checks": {
            "vertex_census_is_81": len(VERTEX_EXPONENTS) == 81,
            "typed_edge_census_is_324": len(edges) == 324,
            "order_three_cycle_census_is_108": len(cycles) == 108,
            "commutation_square_census_is_486": len(squares) == 486,
            "committed_Cayley_job_comparison": committed_comparison,
            "support_joins_use_explicit_tuple_pairs": True,
            "generator_template_joins_use_exponent_tuples": True,
        },
        "availability_policy": {
            "missing_generator_artifacts_block_registry_emission": False,
            "optional_preflight_is_a_prerequisite": False,
            "legacy_source_hash_mismatch_blocks_internal_artifact_loading": False,
            "rerun_discovers_new_completed_manifests": True,
        },
        "numeric_claim_boundary": {
            "support_local_rectangular_CSR_records_are_loadable_on_demand": (
                loaded > 0
            ),
            "support_local_maps_asserted_invertible": False,
            "typed_edge_common_bases_constructed": False,
            "typed_edge_matrices_rebased": False,
            "cycle_matrix_products_computed": False,
            "square_matrix_products_computed": False,
            "strict_K_coherence_asserted": False,
            "K_coherence_up_to_chain_homotopy_asserted": False,
            "global_Phi_or_Yoneda_attempted": False,
        },
        "next_numeric_scope": (
            "construct common source/target bases for the typed edges, rebase "
            "the rectangular generator templates, and only then compute the "
            "108 cycle and 486 square products"
        ),
    }
    return {
        "states": states,
        "artifacts": artifacts,
        "support_bindings": support_bindings,
        "edges": edges,
        "cycles": cycles,
        "squares": squares,
        "certificate": certificate,
    }


GENERATOR_STATUS_FIELDS = (
    "generator",
    "generator_exponents",
    "output_directory",
    "preflight_available",
    "preflight_is_prerequisite",
    "artifact_status",
    "artifact_loaded",
    "internal_files_verified",
    "source_binding",
    "support_ids_bound",
    "support_binding_count",
    "support_bindings_sha256",
    "ready_for_common_basis_rebase",
    "common_basis_rebased",
    "numeric_edge_matrix_emitted",
    "validation_error",
)
SUPPORT_BINDING_FIELDS = (
    "generator",
    "generator_exponents",
    "artifact_support_slot",
    "source_support",
    "target_support",
    "support_pair_key",
    "support_bindings_sha256",
)
EDGE_FIELDS = (
    "edge_job_id",
    "source_K_index",
    "target_K_index",
    "source_K_exponents",
    "target_K_exponents",
    "generator_index",
    "generator",
    "generator_exponents",
    "generator_artifact_loaded",
    "support_ids_bound",
    "support_binding_count",
    "support_bindings_sha256",
    "ready_for_common_basis_rebase",
    "common_basis_rebased",
    "numeric_edge_matrix_emitted",
    "edge_chain_map_certified",
    "status",
)
CYCLE_FIELDS = (
    "cycle_id",
    "generator_index",
    "generator",
    "generator_exponents",
    "base_K_index",
    "base_K_exponents",
    "vertex_K_indices",
    "edge_job_ids",
    "generator_templates_available",
    "common_basis_composition_emitted",
    "strict_equality_computed",
    "chain_homotopy_emitted",
    "coherence_closed",
    "status",
)
SQUARE_FIELDS = (
    "square_id",
    "base_K_index",
    "base_K_exponents",
    "left_generator_index",
    "right_generator_index",
    "left_generator",
    "right_generator",
    "endpoint_K_index",
    "endpoint_K_exponents",
    "left_then_right_edge_job_ids",
    "right_then_left_edge_job_ids",
    "generator_templates_available",
    "common_basis_compositions_emitted",
    "strict_equality_computed",
    "chain_homotopy_emitted",
    "coherence_closed",
    "status",
)


def _csv_value(value):
    if isinstance(value, (tuple, list, dict)):
        return json_compact(value)
    return value


def _csv_rows(rows, fields):
    return [
        {field: _csv_value(row.get(field)) for field in fields}
        for row in rows
    ]


def write_artifacts(
    output_directory: Path | str = OUTPUT_DIRECTORY,
    generator_directories: Mapping[str, Path | str] | None = None,
    committed_cayley_jobs: Path | str = COMMITTED_CAYLEY_JOBS,
):
    """Write the partial-aware registries and honest open-status certificate."""

    output = Path(output_directory).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    result = build(generator_directories, committed_cayley_jobs)

    generator_rows = [
        result["states"][generator] for generator in GENERATOR_LABELS
    ]
    files_and_rows = (
        (
            "generator_artifact_status.csv",
            _csv_rows(generator_rows, GENERATOR_STATUS_FIELDS),
            GENERATOR_STATUS_FIELDS,
        ),
        (
            "generator_support_bindings.csv",
            _csv_rows(result["support_bindings"], SUPPORT_BINDING_FIELDS),
            SUPPORT_BINDING_FIELDS,
        ),
        (
            "z2_k_typed_edges.csv",
            _csv_rows(result["edges"], EDGE_FIELDS),
            EDGE_FIELDS,
        ),
        (
            "z2_k_order_three_cycles.csv",
            _csv_rows(result["cycles"], CYCLE_FIELDS),
            CYCLE_FIELDS,
        ),
        (
            "z2_k_generator_commutation_squares.csv",
            _csv_rows(result["squares"], SQUARE_FIELDS),
            SQUARE_FIELDS,
        ),
    )
    written = []
    for name, rows, fields in files_and_rows:
        path = output / name
        write_csv_atomic(path, rows, fields)
        written.append(path)

    certificate_path = output / "k_coherence_adapter_certificate.json"
    certificate = dict(result["certificate"])
    certificate["files"] = {
        "generator_artifact_status": "generator_artifact_status.csv",
        "generator_support_bindings": "generator_support_bindings.csv",
        "typed_edges": "z2_k_typed_edges.csv",
        "order_three_cycles": "z2_k_order_three_cycles.csv",
        "generator_commutation_squares": (
            "z2_k_generator_commutation_squares.csv"
        ),
    }
    write_json_atomic(certificate_path, certificate)
    written.append(certificate_path)

    manifest = {
        "schema": ADAPTER_MANIFEST_SCHEMA,
        "producer": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "files": {
            path.name: {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in written
        },
    }
    write_json_atomic(output / "manifest.json", manifest)
    return certificate


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Index available K-generator tree artifacts and emit open "
            "edge/cycle/square coherence registries."
        )
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=OUTPUT_DIRECTORY,
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    certificate = write_artifacts(args.output_directory)
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
