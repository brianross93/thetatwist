"""Run the certified tree contraction for one standard K generator at a time.

The underlying producers select a K edge through their module-level
``GENERATOR_EXPONENTS`` value.  This thin driver binds the same exponent tuple
in all three producers for the duration of one run and restores the previous
state afterward.  It does not relax or replace any check in those producers.

Examples
--------

Inspect the planned output without running the expensive calculation::

    python ..._driver.py --generator e1 --dry-run

Run any standard generator::

    python ..._driver.py --generator e0

The standard generators have different overlap lengths.  The shared
generator-general residue transport handles those differences.  ``--preflight``
is an optional, lightweight transport check; it is advisory evidence and is
not a prerequisite for the full producer.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from contextlib import contextmanager
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_overlap_j1 as overlap_j1,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_shamash as generator_shamash,
)
from experiments.section12_instantiation import (
    cm_bridge_hs_rank2_k_generator_tree_contraction as tree_contraction,
)


GENERATOR_EXPONENTS = {
    "e0": (1, 0, 0, 0),
    "e1": (0, 1, 0, 0),
    "e2": (0, 0, 1, 0),
    "e3": (0, 0, 0, 1),
}
GENERATOR_SUPPORTS = {
    "e0": 35,
    "e1": 31,
    "e2": 35,
    "e3": 32,
}
GENERATOR_NORMALIZATION_DATA = {
    "e0": {
        "inverse_adapted_variables": (1, 7),
        "adapted_orientation_scalar": 6,
        "B_determinant_scalar": 6,
        "H_normalization": 17,
        "H_normalization_fourth_power": 7,
    },
    "e1": {
        "inverse_adapted_variables": (1, 5, 7),
        "adapted_orientation_scalar": 8,
        "B_determinant_scalar": 6,
        "H_normalization": 11,
        "H_normalization_fourth_power": 9,
    },
    "e2": {
        "inverse_adapted_variables": (3, 5),
        "adapted_orientation_scalar": 18,
        "B_determinant_scalar": 6,
        "H_normalization": 17,
        "H_normalization_fourth_power": 7,
    },
    "e3": {
        "inverse_adapted_variables": (3, 5, 7),
        "adapted_orientation_scalar": 12,
        "B_determinant_scalar": 1,
        "H_normalization": 4,
        "H_normalization_fourth_power": 8,
    },
}
RESULTS_PARENT = (
    tree_contraction.ROOT / "results" / "section12_instantiation"
)
OUTPUT_DIRECTORIES = {
    "e0": tree_contraction.OUTPUT_DIRECTORY,
    "e1": RESULTS_PARENT
    / "cm_bridge_hs_rank2_k_generator_tree_contraction_e1",
    "e2": RESULTS_PARENT
    / "cm_bridge_hs_rank2_k_generator_tree_contraction_e2",
    "e3": RESULTS_PARENT
    / "cm_bridge_hs_rank2_k_generator_tree_contraction_e3",
}


class KGeneratorTreeContractionDriverError(RuntimeError):
    """A generator label or cross-module binding is inconsistent."""


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def preflight_path(generator, output_directory=None):
    output = (
        Path(output_directory).expanduser().resolve()
        if output_directory is not None
        else OUTPUT_DIRECTORIES[generator]
    )
    return output / "generator_preflight_certificate.json"


def preflight_source_hashes():
    residue_path = (
        ROOT
        / "experiments/section12_instantiation"
        / "cm_bridge_hs_rank2_k_generator_residue_transport.py"
    )
    return {
        "overlap": sha256(overlap_j1.__file__),
        "shamash": sha256(generator_shamash.__file__),
        "tree": sha256(tree_contraction.__file__),
        "residue_transport": sha256(residue_path),
        "driver": sha256(__file__),
    }


def valid_preflight(generator, output_directory=None):
    path = preflight_path(generator, output_directory)
    if not path.is_file():
        return False
    try:
        certificate = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    dimensions = predicted_dimensions(generator)
    return (
        certificate.get("schema")
        == "section12-HS-rank2-K-generator-preflight-v1"
        and certificate.get("status") == "PASS"
        and certificate.get("generator") == generator
        and tuple(certificate.get("K_generator_exponents", ()))
        == GENERATOR_EXPONENTS[generator]
        and int(certificate.get("overlap_factors", -1))
        == dimensions["overlap_factors"]
        and int(certificate.get("overlap_length", -1))
        == dimensions["overlap_length"]
        and certificate.get("source_sha256") == preflight_source_hashes()
    )


def predicted_dimensions(generator):
    """Return the dimension-driven source and fixed target complex sizes."""

    factors = GENERATOR_SUPPORTS[generator]
    overlap_length = 2 * factors
    source_f1 = 7 * overlap_length + 1
    source_f2 = 28 * overlap_length + 6
    source_f3 = 56 * overlap_length + 6 * source_f1
    return {
        "overlap_factors": factors,
        "overlap_length": overlap_length,
        "source_F0": overlap_length,
        "source_F1": source_f1,
        "source_F2": source_f2,
        "source_F3": source_f3,
        "source_L0_rows": source_f1 + 1,
        "target_F0": 98,
        "target_F1": 687,
        "target_F2": 2750,
        "target_F3": 9610,
        "target_L0_columns": 688,
    }


def existing_e0_bytes():
    output = OUTPUT_DIRECTORIES["e0"]
    if not output.is_dir():
        return None
    return sum(path.stat().st_size for path in output.iterdir() if path.is_file())


def existing_artifact_status(generator, output_directory=None):
    """Classify an existing artifact without turning the result into a gate."""

    output = (
        Path(output_directory).expanduser().resolve()
        if output_directory is not None
        else OUTPUT_DIRECTORIES[generator]
    )
    manifest_path = output / "manifest.json"
    if not manifest_path.is_file():
        return "ABSENT"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        recorded = manifest["generator"]["sha256"]
    except (KeyError, OSError, TypeError, ValueError):
        return "UNREADABLE_MANIFEST"
    if recorded == sha256(tree_contraction.__file__):
        return "CURRENT_SOURCE_BOUND"
    return "LEGACY_CURRENT_SOURCE_MISMATCH"


def run_plan(generator, output_directory=None):
    if generator not in GENERATOR_EXPONENTS:
        raise KGeneratorTreeContractionDriverError(
            f"unknown standard generator {generator!r}"
        )
    output = (
        Path(output_directory).expanduser().resolve()
        if output_directory is not None
        else OUTPUT_DIRECTORIES[generator]
    )
    reference_bytes = existing_e0_bytes()
    preflight_current = valid_preflight(generator, output)
    return {
        "generator": generator,
        "exponents": list(GENERATOR_EXPONENTS[generator]),
        "output_directory": str(output),
        "runs_requested": 1,
        "predicted_dimensions": predicted_dimensions(generator),
        "normalization_data": GENERATOR_NORMALIZATION_DATA[generator],
        "execution_ready": True,
        "execution_blocker": None,
        "optional_preflight": {
            "required": False,
            "current_source_hash_bound": preflight_current,
            "path": str(preflight_path(generator, output)),
        },
        "existing_artifact_status": existing_artifact_status(generator, output),
        "reference_e0_artifact_bytes": reference_bytes,
        "estimated_output_bytes": reference_bytes,
        "estimate_basis": (
            "current completed e0 artifact footprint"
            if reference_bytes is not None
            else "unavailable until the e0 artifact exists"
        ),
        "certificate_scope": (
            "identical one-generator finite-overlap checks; no K coherence, "
            "cross-RHom/Yoneda, Question 11.4, or Hodge claim"
        ),
    }


@contextmanager
def bind_generator(generator, output_directory):
    """Bind one exponent tuple across the three existing exact producers."""

    exponents = GENERATOR_EXPONENTS[generator]
    output_directory = Path(output_directory)
    previous = {
        "overlap_exponents": overlap_j1.GENERATOR_EXPONENTS,
        "shamash_exponents": generator_shamash.GENERATOR_EXPONENTS,
        "tree_exponents": tree_contraction.GENERATOR_EXPONENTS,
        "tree_output": tree_contraction.OUTPUT_DIRECTORY,
    }
    overlap_j1.GENERATOR_EXPONENTS = exponents
    generator_shamash.GENERATOR_EXPONENTS = exponents
    tree_contraction.GENERATOR_EXPONENTS = exponents
    tree_contraction.OUTPUT_DIRECTORY = output_directory
    try:
        if not (
            overlap_j1.GENERATOR_EXPONENTS
            == generator_shamash.GENERATOR_EXPONENTS
            == tree_contraction.GENERATOR_EXPONENTS
            == exponents
        ):
            raise KGeneratorTreeContractionDriverError(
                "the generator exponent binding split across producers"
            )
        yield
    finally:
        overlap_j1.GENERATOR_EXPONENTS = previous["overlap_exponents"]
        generator_shamash.GENERATOR_EXPONENTS = previous["shamash_exponents"]
        tree_contraction.GENERATOR_EXPONENTS = previous["tree_exponents"]
        tree_contraction.OUTPUT_DIRECTORY = previous["tree_output"]


def run_one(generator, output_directory=None):
    plan = run_plan(generator, output_directory)
    output = Path(plan["output_directory"])
    with bind_generator(generator, output):
        certificate = tree_contraction.write_artifacts(output)
    if tuple(certificate["K_generator_exponents"]) != GENERATOR_EXPONENTS[generator]:
        raise KGeneratorTreeContractionDriverError(
            "the emitted certificate names the wrong generator"
        )
    dimensions = predicted_dimensions(generator)
    counts = certificate.get("counts", {})
    normalization = certificate.get("generator_normalization", {})
    if (
        int(counts.get("overlap_supports", -1))
        != dimensions["overlap_factors"]
        or int(counts.get("overlap_length", -1))
        != dimensions["overlap_length"]
        or int(normalization.get("support_pairs", -1))
        != dimensions["overlap_factors"]
        or not normalization.get(
            "target_lambda_normalization_pairs_sha256"
        )
    ):
        raise KGeneratorTreeContractionDriverError(
            "the emitted support-count or normalization metadata is wrong"
        )
    return certificate


def preflight_one(generator, output_directory=None):
    """Run transport/dimension checks only; never launch the tree producer."""

    if generator not in GENERATOR_EXPONENTS:
        raise KGeneratorTreeContractionDriverError(
            f"unknown standard generator {generator!r}"
        )
    # Lazy import avoids the residue module's import of this label registry
    # during normal driver startup.
    from experiments.section12_instantiation import (
        cm_bridge_hs_rank2_k_generator_residue_transport as residue_transport,
    )

    output = (
        Path(output_directory).expanduser().resolve()
        if output_directory is not None
        else OUTPUT_DIRECTORIES[generator]
    )
    result = residue_transport.build(generator)
    transport = result["certificate"]
    dimensions = predicted_dimensions(generator)
    if (
        int(transport["overlap_factors"])
        != dimensions["overlap_factors"]
        or int(transport["overlap_length"])
        != dimensions["overlap_length"]
        or not transport["whole_overlap_top_residue_identity"]
        or not transport["whole_overlap_normalization_factorization"]
        or not transport["orientation_coboundary_formula_exact"]
        or not transport["target_normalization_is_unit"]
    ):
        raise KGeneratorTreeContractionDriverError(
            "the residue transport failed the generator preflight"
        )
    certificate = {
        "schema": "section12-HS-rank2-K-generator-preflight-v1",
        "status": "PASS",
        "generator": generator,
        "K_generator_exponents": list(GENERATOR_EXPONENTS[generator]),
        "overlap_factors": dimensions["overlap_factors"],
        "overlap_length": dimensions["overlap_length"],
        "predicted_dimensions": dimensions,
        "normalization": {
            "target_normalization_pairs_sha256": transport[
                "target_normalization_pairs_sha256"
            ],
            "residue_transport_pairs_sha256": transport[
                "residue_transport_pairs_sha256"
            ],
            "orientation_coboundary_pairs_sha256": transport[
                "orientation_coboundary_pairs_sha256"
            ],
            "pair_ledger_hash_encoding": transport[
                "pair_ledger_hash_encoding"
            ],
            "strict_Laurent_lift": "OPEN",
        },
        "source_sha256": preflight_source_hashes(),
        "full_tree_producer_launched": False,
    }
    output.mkdir(parents=True, exist_ok=True)
    preflight_path(generator, output).write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if not valid_preflight(generator, output):
        raise KGeneratorTreeContractionDriverError(
            "the persisted preflight did not validate"
        )
    return certificate


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Run the exact finite-overlap tree contraction for one standard "
            "K generator."
        )
    )
    parser.add_argument(
        "--generator",
        required=True,
        choices=tuple(GENERATOR_EXPONENTS),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="optional distinct output directory",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the bound generator, output path, and disk estimate only",
    )
    parser.add_argument(
        "--preflight",
        action="store_true",
        help="run transport and dimension checks without launching the tree producer",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    plan = run_plan(args.generator, args.output)
    if args.dry_run:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return
    if args.preflight:
        certificate = preflight_one(args.generator, args.output)
        print(json.dumps(certificate, indent=2, sort_keys=True))
        return
    certificate = run_one(args.generator, args.output)
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
