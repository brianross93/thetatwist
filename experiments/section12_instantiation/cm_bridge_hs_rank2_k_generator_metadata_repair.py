"""Write a non-destructive pair-ledger hash repair sidecar.

The first full ``e1`` artifact hashed identical dual-number pair ledgers in two
different numeric encodings.  Its numeric arrays and manifest file bindings are
unchanged.  This module reads those arrays, verifies the original manifest,
and records both the original uint8 hash and one canonical signed little-endian
int64 hash.  It does not rewrite any original artifact.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_residue_transport as residue_transport,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_shamash as generator_shamash,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_tree_contraction as tree_contraction,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_tree_contraction_driver as driver,
)


PAIR_KEYS = (
    "target_lambda_normalization_pairs",
    "residue_transport_pairs",
    "raw_counit_weight_pairs",
    "orientation_coboundary_pairs",
)
CERTIFICATE_HASH_KEYS = {
    "target_lambda_normalization_pairs": (
        "target_lambda_normalization_pairs_sha256"
    ),
    "residue_transport_pairs": "residue_transport_pairs_sha256",
    "raw_counit_weight_pairs": "raw_counit_weight_pairs_sha256",
    "orientation_coboundary_pairs": (
        "orientation_coboundary_pairs_sha256"
    ),
}
PREFLIGHT_HASH_KEYS = {
    "target_lambda_normalization_pairs": (
        "target_normalization_pairs_sha256"
    ),
    "residue_transport_pairs": "residue_transport_pairs_sha256",
    "orientation_coboundary_pairs": (
        "orientation_coboundary_pairs_sha256"
    ),
}


class GeneratorMetadataRepairError(RuntimeError):
    """An original artifact or its pair-ledger binding changed."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json_atomic(path: Path, value) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def build_certificate(output_directory: Path, generator: str):
    output = Path(output_directory).expanduser().resolve()
    if generator not in driver.GENERATOR_EXPONENTS:
        raise GeneratorMetadataRepairError(
            f"unknown standard generator {generator!r}"
        )
    manifest_path = output / "manifest.json"
    certificate_path = output / "tree_contraction_certificate.json"
    maps_path = output / "tree_contracted_support_maps.npz"
    preflight_path = output / "generator_preflight_certificate.json"
    for path in (manifest_path, certificate_path, maps_path, preflight_path):
        if not path.is_file():
            raise GeneratorMetadataRepairError(
                f"required original artifact is absent: {path}"
            )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    original_certificate = json.loads(
        certificate_path.read_text(encoding="utf-8")
    )
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    expected_exponents = list(driver.GENERATOR_EXPONENTS[generator])
    if original_certificate["K_generator_exponents"] != expected_exponents:
        raise GeneratorMetadataRepairError(
            "the original certificate names a different generator"
        )
    if preflight["K_generator_exponents"] != expected_exponents:
        raise GeneratorMetadataRepairError(
            "the current preflight names a different generator"
        )

    verified_files = {}
    for name, binding in manifest["files"].items():
        path = output / name
        if not path.is_file():
            raise GeneratorMetadataRepairError(
                f"manifest-bound file is absent: {path}"
            )
        actual = {
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        if actual != binding:
            raise GeneratorMetadataRepairError(
                f"manifest binding changed for {name}"
            )
        verified_files[name] = actual

    support_count = int(original_certificate["counts"]["overlap_supports"])
    ledgers = {}
    with np.load(maps_path, allow_pickle=False) as maps:
        if maps["generator_exponents"].astype(int).tolist() != expected_exponents:
            raise GeneratorMetadataRepairError(
                "the numeric map artifact names a different generator"
            )
        for key in PAIR_KEYS:
            array = np.asarray(maps[key])
            if array.shape != (support_count, 2):
                raise GeneratorMetadataRepairError(
                    f"{key} has shape {array.shape}, expected "
                    f"{(support_count, 2)}"
                )
            original_uint8_hash = tree_contraction.array_sha256(array)
            canonical_hash = generator_shamash.pair_ledger_sha256(array)
            certificate_hash = original_certificate[
                "generator_normalization"
            ][CERTIFICATE_HASH_KEYS[key]]
            if original_uint8_hash != certificate_hash:
                raise GeneratorMetadataRepairError(
                    f"the original uint8 certificate hash changed for {key}"
                )
            preflight_key = PREFLIGHT_HASH_KEYS.get(key)
            preflight_hash = (
                preflight["normalization"][preflight_key]
                if preflight_key is not None
                else None
            )
            if preflight_hash is not None and canonical_hash != preflight_hash:
                raise GeneratorMetadataRepairError(
                    f"the canonical preflight hash changed for {key}"
                )
            ledgers[key] = {
                "shape": list(array.shape),
                "stored_dtype": str(array.dtype),
                "original_uint8_sha256": original_uint8_hash,
                "canonical_sha256": canonical_hash,
                "canonical_encoding": (
                    "signed-little-endian-int64-C-order"
                ),
                "matches_current_preflight": (
                    None
                    if preflight_hash is None
                    else canonical_hash == preflight_hash
                ),
            }

    target_hash = ledgers[
        "target_lambda_normalization_pairs"
    ]["canonical_sha256"]
    if (
        target_hash
        != original_certificate["counts"][
            "target_lambda_normalization_pairs_sha256"
        ]
    ):
        raise GeneratorMetadataRepairError(
            "the canonical target-normalization hash does not match counts"
        )

    current_sources = {
        "residue_transport": {
            "path": str(
                Path(residue_transport.__file__).resolve().relative_to(ROOT)
            ),
            "sha256": sha256(Path(residue_transport.__file__).resolve()),
        },
        "shamash": {
            "path": str(
                Path(generator_shamash.__file__).resolve().relative_to(ROOT)
            ),
            "sha256": sha256(Path(generator_shamash.__file__).resolve()),
        },
        "tree": {
            "path": str(
                Path(tree_contraction.__file__).resolve().relative_to(ROOT)
            ),
            "sha256": sha256(Path(tree_contraction.__file__).resolve()),
        },
        "driver": {
            "path": str(Path(driver.__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(driver.__file__).resolve()),
        },
    }
    return {
        "schema": "section12-HS-rank2-K-generator-metadata-repair-v1",
        "status": "PASS",
        "generator": generator,
        "K_generator_exponents": expected_exponents,
        "repair_scope": (
            "pair-ledger hash encoding only; all original numeric artifacts, "
            "certificate, and manifest remain byte-for-byte unchanged"
        ),
        "original_artifact": {
            "directory": str(output),
            "manifest": {
                "path": manifest_path.name,
                "sha256": sha256(manifest_path),
            },
            "generator_source_binding": manifest["generator"],
            "verified_files": verified_files,
        },
        "current_sources": current_sources,
        "current_preflight": {
            "path": preflight_path.name,
            "sha256": sha256(preflight_path),
            "source_sha256": preflight["source_sha256"],
        },
        "pair_ledgers": ledgers,
        "consistency": {
            "manifest_files_verified": len(verified_files),
            "all_original_uint8_hashes_match": True,
            "all_available_canonical_preflight_hashes_match": True,
            "target_hash_matches_original_counts": True,
            "numeric_artifacts_rewritten": False,
        },
        "strict_Laurent_lift": "OPEN",
    }


def write_certificate(output_directory: Path, generator: str):
    output = Path(output_directory).expanduser().resolve()
    certificate = build_certificate(output, generator)
    path = output / "metadata_repair_certificate.json"
    write_json_atomic(path, certificate)
    return path, certificate


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Write a non-destructive pair-ledger metadata repair."
    )
    parser.add_argument(
        "--generator",
        required=True,
        choices=tuple(driver.GENERATOR_EXPONENTS),
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    path, certificate = write_certificate(args.output, args.generator)
    print(
        json.dumps(
            {
                "path": str(path),
                "certificate": certificate,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
