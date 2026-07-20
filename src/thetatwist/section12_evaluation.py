"""Exact handoff for the corrected Section 12 evaluation problem.

This module does not invent the missing Atiyah map.  It records the exact
theta-divisor incidence complex, proves its counting ledger, separates the ten
ordinary polarized deformation directions from the six mixed Poisson--gerbe
directions, and emits the finite input contract for the remaining calculation.
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

from .core import Character, corrected_character, exponential, product


SCHEMA_VERSION = "thetatwist-section12-evaluation-handoff-2"


def _check_d(d: int) -> int:
    if d < 3 or d % 2 == 0:
        raise ValueError("Section 12 assumes odd d >= 3")
    return (d + 9) // 2


def _character(values: Iterable[int | Fraction]) -> Character:
    result = tuple(Fraction(value) for value in values)
    if len(result) != 5:
        raise ValueError("a theta character has five coefficients")
    return result  # type: ignore[return-value]


def _add(left: Character, right: Character, scale: int | Fraction = 1) -> Character:
    factor = Fraction(scale)
    return _character(x + factor * y for x, y in zip(left, right, strict=True))


def _complete_intersection_character(codimension: int) -> Character:
    """Return ``ch(O_(D_1 intersection ... intersection D_c))``.

    Every divisor has class ``Theta``.  A Koszul resolution gives
    ``(1-exp(-Theta))^codimension``.
    """

    if not 1 <= codimension <= 4:
        raise ValueError("codimension must be between 1 and 4")
    factor = _add(_character((1, 0, 0, 0, 0)), exponential(-1), -1)
    result = _character((1, 0, 0, 0, 0))
    for _ in range(codimension):
        result = product(result, factor)
    return result


def cyclic_distance(i: int, j: int, order: int) -> int:
    """Return the unoriented distance between two vertices of a cycle."""

    if order < 3:
        raise ValueError("order must be at least three")
    if not 0 <= i < order or not 0 <= j < order:
        raise ValueError("cycle index is outside the declared order")
    delta = (i - j) % order
    return min(delta, order - delta)


def section12_incidence(d: int = 3) -> dict[str, Any]:
    """Enumerate the incidence types of the corrected theta-cycle.

    The components are ``Z_i=D_i intersection D_(i+1)``.  Genericity makes an
    intersection empty when the union of the divisor labels has size greater
    than four.
    """

    n = _check_d(d)
    edges = tuple(frozenset((i, (i + 1) % n)) for i in range(n))

    adjacent_pairs: list[tuple[int, int]] = []
    point_pairs: list[tuple[int, int]] = []
    isolated_pairs: list[tuple[int, int]] = []
    for i, j in combinations(range(n), 2):
        divisor_count = len(edges[i] | edges[j])
        if divisor_count == 3:
            adjacent_pairs.append((i, j))
        elif divisor_count == 4:
            point_pairs.append((i, j))
            if cyclic_distance(i, j, n) > 2:
                isolated_pairs.append((i, j))

    triple_points: list[tuple[int, int, int]] = []
    for indices in combinations(range(n), 3):
        divisor_count = len(frozenset().union(*(edges[index] for index in indices)))
        if divisor_count <= 4:
            if divisor_count != 4:
                raise ArithmeticError("an unexpected triple-incidence type occurred")
            triple_points.append(indices)

    point_length = math.factorial(4)
    total_isolated_nodes = point_length * len(isolated_pairs)
    normalized_nodes = (d + 9) * (2 * d - 1)
    unnormalized_isolated_nodes = total_isolated_nodes - normalized_nodes

    expected = {
        "adjacent_pairs": n,
        "point_pairs": n * (n - 3) // 2,
        "triple_points": n,
        "isolated_pairs": n * (n - 5) // 2,
        "total_isolated_nodes": 3 * (d + 9) * (d - 1),
        "normalized_nodes": (d + 9) * (2 * d - 1),
        "unnormalized_isolated_nodes": (d + 9) * (d - 2),
    }
    observed = {
        "adjacent_pairs": len(adjacent_pairs),
        "point_pairs": len(point_pairs),
        "triple_points": len(triple_points),
        "isolated_pairs": len(isolated_pairs),
        "total_isolated_nodes": total_isolated_nodes,
        "normalized_nodes": normalized_nodes,
        "unnormalized_isolated_nodes": unnormalized_isolated_nodes,
    }
    if observed != expected:
        raise ArithmeticError("the Section 12 incidence ledger did not close")

    return {
        "d": d,
        "component_count": n,
        "components": [sorted(edge) for edge in edges],
        "adjacent_component_pairs": [list(pair) for pair in adjacent_pairs],
        "point_intersection_pairs": [list(pair) for pair in point_pairs],
        "isolated_double_point_pairs": [list(pair) for pair in isolated_pairs],
        "consecutive_triple_intersections": [list(row) for row in triple_points],
        "point_set_length": point_length,
        **observed,
        "selection_boundary": (
            "The source specifies only the total number of normalized isolated "
            "nodes. The evaluation calculation must supply the selected nodes "
            "and the normalization extension maps."
        ),
    }


def cech_character(d: int = 3) -> Character:
    """Recover the corrected character from the exact Cech incidence ledger."""

    incidence = section12_incidence(d)
    n = incidence["component_count"]
    surface = _complete_intersection_character(2)
    curve = _complete_intersection_character(3)
    point_set = _complete_intersection_character(4)
    point = _character((0, 0, 0, 0, Fraction(1, math.factorial(4))))

    # 0 -> O_Z -> C0 -> C1 -> C2 -> 0.
    structure_z = _add(_character((0, 0, 0, 0, 0)), surface, n)
    structure_z = _add(structure_z, curve, -incidence["adjacent_pairs"])
    structure_z = _add(structure_z, point_set, -incidence["point_pairs"])
    structure_z = _add(structure_z, point_set, incidence["triple_points"])

    # 0 -> O_Z -> nu_* O_Ztilde -> Q_R -> 0.
    normalized = _add(structure_z, point, incidence["normalized_nodes"])
    untwisted_object = _add(_character((1, 0, 0, 0, 0)), normalized, -1)
    return product(untwisted_object, exponential(3))


def cech_contract(d: int = 3) -> dict[str, Any]:
    """Return the finite derived presentation that the CAS must materialize."""

    incidence = section12_incidence(d)
    return {
        "union_augmentation": "0 -> O_Z -> C0 -> C1 -> C2 -> 0",
        "union_terms": [
            {
                "id": "C0",
                "summands": incidence["component_count"],
                "type": "O_(D_i intersection D_(i+1))",
                "koszul_codimension": 2,
            },
            {
                "id": "C1_curves",
                "summands": incidence["adjacent_pairs"],
                "type": "O_(three theta divisors)",
                "koszul_codimension": 3,
            },
            {
                "id": "C1_points",
                "summands": incidence["point_pairs"],
                "type": "O_(four theta divisors), length 24 per summand",
                "koszul_codimension": 4,
            },
            {
                "id": "C2_points",
                "summands": incidence["triple_points"],
                "type": "O_(four theta divisors), Cech triple correction",
                "koszul_codimension": 4,
            },
        ],
        "partial_normalization": {
            "sequence": "0 -> O_Z -> nu_*O_Ztilde -> Q_R -> 0",
            "length_Q_R": incidence["normalized_nodes"],
            "required_data": (
                "selected reduced isolated nodes and the extension/lift maps; "
                "the K-class alone does not determine the derived object"
            ),
        },
        "normalized_incidence_resolution": {
            "sequence": "0 -> nu_*O_Ztilde -> N0 -> N1 -> N2 -> 0",
            "N0": {
                "summands": incidence["component_count"],
                "type": "O_(D_i intersection D_(i+1))",
            },
            "N1": {
                "adjacent_curve_summands": incidence["adjacent_pairs"],
                "consecutive_four_theta_summands": incidence["triple_points"],
                "residual_isolated_node_length": incidence[
                    "unnormalized_isolated_nodes"
                ],
                "type": (
                    "adjacent curves plus the six four-theta Q_i schemes plus "
                    "the selected residual node scheme R_keep"
                ),
            },
            "N2": {
                "consecutive_four_theta_summands": incidence["triple_points"],
                "type": "Cech triple corrections on the Q_i schemes",
            },
            "selection_required": (
                "R_keep is the reduced residual node scheme. Its complement in "
                "the isolated nodes is normalized."
            ),
        },
        "corrected_object": "E=[O_X -> nu_*O_Ztilde] tensor O_X(3 Theta)",
        "character_matches_alpha_plus_3beta": cech_character(d)
        == corrected_character(d),
    }


def ordinary_kernel_witnesses() -> tuple[str, ...]:
    """Return the ten polarized-deformation witnesses in ``H^1(T_X)``."""

    diagonal = tuple(f"B_{i}{i}" for i in range(1, 5))
    symmetric = tuple(
        f"B_{i}{j}+B_{j}{i}" for i, j in combinations(range(1, 5), 2)
    )
    return diagonal + symmetric


def mixed_kernel_witnesses(d: int = 3) -> tuple[str, ...]:
    """Return the six mixed Poisson--gerbe witnesses still to evaluate."""

    _check_d(d)
    return tuple(f"C_{i}{j}-{d}A_{i}{j}" for i, j in combinations(range(1, 5), 2))


def evaluation_handoff(d: int = 3) -> dict[str, Any]:
    """Build the corrected-only handoff for the missing evaluation map."""

    incidence = section12_incidence(d)
    ordinary = ordinary_kernel_witnesses()
    mixed = mixed_kernel_witnesses(d)
    if len(ordinary) != 10 or len(mixed) != 6:
        raise ArithmeticError("the 10+6 Hodge-kernel split did not close")

    return {
        "schema_version": SCHEMA_VERSION,
        "source": "arXiv:2509.23403v2, Section 12 and Question 11.4",
        "scope": "corrected Section 12 alpha+3beta object only",
        "d": d,
        "object": "E=[O_X -> nu_*O_Ztilde] tensor O_X(3 Theta)",
        "incidence": incidence,
        "parameter_incidence_closure": {
            "name": "parameter-incidence closure check",
            "symbolic_schema_status": "passed",
            "instantiated_fiber_status": "not yet computed",
            "component_equation": "n=(d+9)/2",
            "component_parameter": incidence["component_count"],
            "realized_surface_count": len(incidence["components"]),
            "normalization_equation": "r=(d+9)(2d-1)",
            "normalization_parameter": incidence["normalized_nodes"],
            "realized_normalized_node_count": incidence["normalized_nodes"],
            "total_isolated_nodes": incidence["total_isolated_nodes"],
            "residual_equation": "total-r=(d+9)(d-2)",
            "realized_residual_node_count": incidence[
                "unnormalized_isolated_nodes"
            ],
            "strata": {
                "surfaces_Z_i": incidence["component_count"],
                "adjacent_triple_theta_curves_T_i": incidence["adjacent_pairs"],
                "consecutive_four_theta_schemes_Q_i": incidence["triple_points"],
                "length_each_Q_i": incidence["point_set_length"],
                "isolated_four_theta_schemes": incidence["isolated_pairs"],
                "length_each_isolated_scheme": incidence["point_set_length"],
            },
            "character_realization": "Cech-Koszul character is alpha+3beta",
            "interpretation": (
                "The corrected character, the printed n and r values, and the "
                "exact incidence schema agree. This is a separate realization "
                "check, not a claim that the printed inputs are logically "
                "independent or that an explicit finite-field fiber is already "
                "constructed."
            ),
        },
        "cech_contract": cech_contract(d),
        "hodge_kernel": {
            "dimension": 16,
            "ordinary_polarized_dimension": len(ordinary),
            "ordinary_witnesses": list(ordinary),
            "ordinary_status": (
                "closed geometrically: extend the polarized abelian scheme, "
                "translation sections, relative theta divisors, selected etale "
                "node sections, and the relative partial normalization"
            ),
            "mixed_dimension": len(mixed),
            "mixed_witnesses": list(mixed),
            "mixed_status": "open: compute each image in Ext^2(E,E)",
        },
        "completion_test": {
            "required_equalities": [f"ev_E({witness})=0" for witness in mixed],
            "consequence": (
                "Together with ker(ev_E) contained in ker(contraction by ch(E)), "
                "these six equalities give equality of the two kernels and the "
                "weaker criterion in Question 11.4."
            ),
        },
        "preregistered_outcomes": {
            "all_six_zero": (
                "For this corrected object and recorded node choice, the ten "
                "geometric vanishings and six computed vanishings give equality "
                "of the evaluation and Chern-contraction kernels. This proves "
                "the weaker criterion in Question 11.4 for this object; it does "
                "not assert surjectivity."
            ),
            "any_one_nonzero": (
                "For this corrected object and recorded node choice, a nonzero "
                "mixed evaluation proves failure of the weaker criterion. It "
                "does not cover other residual-node choices unless independence "
                "from that choice is also proved."
            ),
            "certificate_forms": {
                "zero": "an exact End(P) cochain h with D_End(h)=z_ij",
                "nonzero": (
                    "a dual cocycle ell with ell*D_End=0 and ell(z_ij) nonzero"
                ),
            },
            "integrity_failure": (
                "No mathematical verdict: the run did not instantiate the "
                "declared object."
            ),
            "field_scope": (
                "An all-zero finite-field result is a result over that field. "
                "A complex conclusion requires the fixed-model lift or exact "
                "reconstruction."
            ),
        },
        "symmetry_boundary": (
            "Do not reduce the six mixed checks to one representative unless a "
            "chain-level action on this fixed divisor and node configuration is proved."
        ),
        "symmetry_analysis": {
            "translation_use": (
                "Use cyclic translations to reuse component and node assembly. "
                "Translations act trivially on the HT^2 source and do not merge "
                "the six mixed witnesses."
            ),
            "inversion_use": (
                "Inversion fixes the two-factor mixed classes and does not merge "
                "the six mixed witnesses."
            ),
            "valid_reduction": (
                "A reduction requires a proved automorphism with nontrivial "
                "linear action on H^1 and a compatible chain-level action. "
                "Compute its actual orbits on the six mixed witnesses."
            ),
        },
        "execution_strategy": {
            "global_model": (
                "Do not eliminate the 3Theta embedding in P^80; h0(3Theta)=81."
            ),
            "coordinate_model": (
                "Use exact theta-coordinate multiplication, restriction, "
                "point-evaluation, and normalization matrices."
            ),
            "order_six_translation": (
                "An exact order-six translation is not automatically a level-3 "
                "clock action because the theta group of 3Theta covers X[3]. "
                "Use direct_sum_i H0(t_(ig)^*3Theta), with translation permuting "
                "the six summands, or use a compatible level-6 theta structure. "
                "Neither construction automatically linearizes E."
            ),
            "finite_field_selection": (
                "p congruent to 1 modulo 12 supplies useful roots of unity but "
                "does not by itself make the torsion or theta structure rational."
            ),
            "fiber_checks": [
                "good reduction",
                "required rational theta data",
                "6g=0, 2g nonzero, and 3g nonzero",
                "the expected 72-node incidence ledger",
                "Frobenius stability of the selected residual node scheme",
            ],
            "optional_cyclic_node_choice": (
                "If node-level cyclic symmetry is used, choose R_keep as two "
                "full free six-element orbits and record them."
            ),
            "characteristic_zero_provenance": (
                "Reduce one fixed integral model at independent good primes; "
                "reconstruct from two primes and verify at a third."
            ),
            "shared_totalization": (
                "Build and cache one total complex and one Atiyah cocycle; "
                "contract all six witnesses against that common data."
            ),
            "macaulay2_role": (
                "Manage supplied finite matrices, complexes, Hom complexes, "
                "Ext reduction, and syzygy certificates rather than global "
                "P^80 elimination."
            ),
        },
        "matrix_inputs": [
            "equations or transition data for the six translated theta divisors",
            "the Cech restriction maps for all listed intersections",
            "the selected normalized nodes and normalization extension maps",
            "a locally free totalization of [O_X -> nu_*O_Ztilde](3Theta)",
            "the principal-parts Atiyah cocycle of that totalization",
            "compatible bases for HT^2 and Ext^2",
        ],
        "claim_boundary": (
            "The incidence and 10+6 reduction are exact. This handoff does not "
            "supply the six mixed Ext^2 values and does not prove the Hodge conjecture."
        ),
    }


def write_evaluation_handoff(directory: str | Path, d: int = 3) -> dict[str, Path]:
    """Write the handoff and a SHA-256 manifest."""

    output = Path(directory)
    output.mkdir(parents=True, exist_ok=True)
    handoff_path = output / f"d{d}_section12_evaluation_handoff.json"
    handoff_path.write_text(
        json.dumps(evaluation_handoff(d), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    manifest_path = output / "section12_evaluation_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema_version": SCHEMA_VERSION,
                "artifacts": [
                    {
                        "path": handoff_path.name,
                        "sha256": hashlib.sha256(handoff_path.read_bytes()).hexdigest(),
                        "bytes": handoff_path.stat().st_size,
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return {"handoff": handoff_path, "manifest": manifest_path}
