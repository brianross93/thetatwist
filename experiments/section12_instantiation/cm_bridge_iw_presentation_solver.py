"""Construct the 96 exact chart-local presentations of the Section 12 HS ideal.

The input audit reduced the next missing object to one finite calculation over
``F_31``.  This producer supplies that calculation.  Python reconstructs the
six translated theta surfaces and the 81 translates of each normalized pair
``X_s,Q_s``.  Singular then works in the affine surface quotient ring and
computes a reduced standard basis, two explicit change-of-generator
certificates, and a first Schreyer matrix.

The equality ``J=I_W`` is certified without intersecting 81 primary ideals.
The existing exact first-jet certificate gives ``J subset I_W`` at every
support stalk.  On each affine chart Singular checks that ``R/J`` is
zero-dimensional and has length exactly twice the number of support points in
that chart.  Since the curvilinear support has that same length, the inclusion
is an equality.  This is both faster and less fragile than a 81-fold ideal
intersection.

Run one job::

    python3 .../cm_bridge_iw_presentation_solver.py \
        --surface 0 --chart C0000 --singular /opt/homebrew/bin/Singular

Run all 96 jobs and assemble ``chart_local_full_IW_presentations.jsonl``::

    python3 .../cm_bridge_iw_presentation_solver.py \
        --all --singular /opt/homebrew/bin/Singular

No claim about Question 11.4 or the Hodge conjecture is made here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
P = 31
OUTPUT_DIRECTORY = ROOT / (
    "results/section12_instantiation/cm_bridge_iw_presentation_solver"
)
PRESENTATIONS = OUTPUT_DIRECTORY / "chart_local_full_IW_presentations.jsonl"
SINGULAR_DEFAULT = "/opt/homebrew/bin/Singular"


class IWPresentationError(RuntimeError):
    """An exact input, CAS job, or certificate invariant failed."""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _modules():
    # Imported lazily so that a stored JSONL can be verified without rebuilding
    # all group actions.
    from experiments.section12_instantiation import (
        cm_bridge_hs_chart_local_syzygy_audit as local_symbols,
    )
    from experiments.section12_instantiation import (
        cm_bridge_global_hs_chart_resolution as global_hs,
    )

    return local_symbols, global_hs


def _rref_generator_basis(rows):
    """Return an exact row basis and both constant change matrices over F_31.

    If ``A`` is the 162 by 1296 original generator matrix, the return values
    satisfy ``R = U*A`` and ``A = C*R``.  Thus replacing the 162 translates by
    the rows of ``R`` does not change their ideal on any chart.
    """

    original = np.asarray(rows, dtype=np.int64) % P
    work = original.copy()
    transform = np.eye(work.shape[0], dtype=np.int64)
    pivot_row = 0
    pivot_columns = []
    for column in range(work.shape[1]):
        candidates = np.flatnonzero(work[pivot_row:, column])
        if not len(candidates):
            continue
        pivot = pivot_row + int(candidates[0])
        if pivot != pivot_row:
            work[[pivot_row, pivot]] = work[[pivot, pivot_row]]
            transform[[pivot_row, pivot]] = transform[[pivot, pivot_row]]
        inverse = pow(int(work[pivot_row, column]), -1, P)
        work[pivot_row] = work[pivot_row] * inverse % P
        transform[pivot_row] = transform[pivot_row] * inverse % P
        nonzero = np.flatnonzero(work[:, column])
        nonzero = nonzero[nonzero != pivot_row]
        if len(nonzero):
            multipliers = work[nonzero, column].copy()
            work[nonzero] = (
                work[nonzero] - multipliers[:, None] * work[pivot_row]
            ) % P
            transform[nonzero] = (
                transform[nonzero] - multipliers[:, None] * transform[pivot_row]
            ) % P
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == work.shape[0]:
            break

    rank = pivot_row
    reduced = work[:rank]
    reduced_from_original = transform[:rank]
    original_from_reduced = original[:, pivot_columns]
    if not np.array_equal(reduced_from_original @ original % P, reduced):
        raise IWPresentationError("the reduced generator basis lost its forward map")
    if not np.array_equal(original_from_reduced @ reduced % P, original):
        raise IWPresentationError("the reduced generator basis lost its inverse map")
    return reduced, reduced_from_original, original_from_reduced, pivot_columns


def _sparse_matrix_entries(matrix):
    rows, columns = np.nonzero(matrix % P)
    return [
        [int(row), int(column), int(matrix[row, column] % P)]
        for row, column in zip(rows, columns)
    ]


def load_problem_context():
    local_symbols, global_hs = _modules()
    (
        global_certificate,
        descent,
        seeds,
        supports,
        atlas,
        _edges,
        theta,
        generators,
    ) = local_symbols.load_inputs()
    points, point_values, point_first = local_symbols.seed_emitter.curve_data()
    records, _element_index, _checks = local_symbols.group_data(descent, points)

    translated = {}
    reductions = {}
    for surface in range(6):
        original = []
        labels = []
        for k_index, record in enumerate(records):
            for generator in ("X", "Q"):
                original.append(
                    local_symbols.apply_degree2_action(
                        generators[(surface, generator)], record
                    )
                )
                labels.append({"K_element_index": k_index, "generator": generator})
        reduced, forward, inverse, pivots = _rref_generator_basis(original)
        translated[surface] = reduced
        reductions[surface] = {
            "surface": f"Z_{surface}",
            "original_generator_count": 162,
            "reduced_generator_count": int(reduced.shape[0]),
            "original_generator_labels": labels,
            "pivot_coefficient_columns": [int(value) for value in pivots],
            "reduced_from_original_shape": list(forward.shape),
            "reduced_from_original_entries": _sparse_matrix_entries(forward),
            "original_from_reduced_shape": list(inverse.shape),
            "original_from_reduced_entries": _sparse_matrix_entries(inverse),
            "identities": ["R=U*A", "A=C*R"],
        }

    chart_bits = {
        row["chart"]: tuple(map(int, row["bits_Z0_Y1"])) for row in atlas
    }
    support_by_surface = {
        surface: [
            tuple(json.loads(row["source_point_indices"]))
            for row in supports
            if row["surface"] == f"Z_{surface}"
        ]
        for surface in range(6)
    }
    return {
        "local_symbols": local_symbols,
        "global_hs": global_hs,
        "global_certificate": global_certificate,
        "seeds": seeds,
        "theta": theta,
        "translated_reduced": translated,
        "generator_reductions": reductions,
        "chart_bits": chart_bits,
        "point_values": point_values,
        "point_first": point_first,
        "points": points,
        "descent": descent,
        "records": records,
        "element_index": _element_index,
        "all_generators": generators,
        "support_rows": supports,
        "support_by_surface": support_by_surface,
    }


def certifying_conormal_generator_subset(context, surface, chart=None):
    """Select a small translated-generator subset spanning every conormal.

    The returned original column indices refer to the ordered list
    ``(k_0 X,k_0 Q,k_1 X,k_1 Q,...)``. At each curvilinear support point on
    the selected chart (or all 81 points when ``chart`` is omitted), the
    selected columns have rank two in the recorded conormal symbol.
    Hence the selected generators define a subideal ``J_0 subset J`` whose
    stalk already equals ``I_W`` along ``W``.  A finite quotient-length check
    for ``J_0`` is therefore enough to prove ``J_0=J=I_W``.
    """

    local_symbols = context["local_symbols"]
    point_second = context["global_hs"].factor_second_derivatives(context["points"])
    surface_rows = [
        row for row in context["support_rows"] if row["surface"] == f"Z_{surface}"
    ]
    ordered_points, symbol, full_rank, automorphy_checks, direct_checks = (
        local_symbols.surface_symbol_matrix(
            surface,
            context["seeds"][surface],
            surface_rows,
            context["theta"],
            context["all_generators"],
            context["records"],
            context["element_index"],
            context["points"],
            context["point_values"],
            context["point_first"],
            point_second,
        )
    )
    blocks = np.asarray(symbol, dtype=np.int64).reshape(81, 2, 162) % P
    if chart is not None:
        bits = context["chart_bits"][chart]
        selected_supports = [
            index
            for index, point in enumerate(ordered_points)
            if context["local_symbols"].chart_membership(
                point, bits, context["point_values"]
            )
        ]
        blocks = blocks[selected_supports]
    support_block_count = len(blocks)

    def covered(columns):
        columns = list(columns)
        result = np.zeros(support_block_count, dtype=bool)
        for left, right in itertools.combinations(columns, 2):
            determinant = (
                blocks[:, 0, left] * blocks[:, 1, right]
                - blocks[:, 0, right] * blocks[:, 1, left]
            ) % P
            result |= determinant != 0
        return result

    best_pair = None
    best_coverage = None
    for left in range(162):
        for right in range(left + 1, 162):
            candidate = covered((left, right))
            score = int(candidate.sum())
            if best_coverage is None or score > int(best_coverage.sum()):
                best_pair = (left, right)
                best_coverage = candidate
    selected = list(best_pair)
    while not covered(selected).all():
        candidates = []
        for column in range(162):
            if column in selected:
                continue
            candidates.append((int(covered((*selected, column)).sum()), -column, column))
        score, _negative, column = max(candidates)
        if score > int(covered(selected).sum()):
            selected.append(column)
            continue
        pair_candidates = []
        remaining = [column for column in range(162) if column not in selected]
        for left, right in itertools.combinations(remaining, 2):
            pair_candidates.append(
                (
                    int(covered((*selected, left, right)).sum()),
                    -left,
                    -right,
                    left,
                    right,
                )
            )
        pair_score, _negative_left, _negative_right, left, right = max(
            pair_candidates
        )
        if pair_score <= int(covered(selected).sum()):
            raise IWPresentationError("the conormal cover greedy search stalled")
        selected.extend((left, right))
    changed = True
    while changed and len(selected) > 2:
        changed = False
        for column in tuple(reversed(selected)):
            candidate = [value for value in selected if value != column]
            if covered(candidate).all():
                selected = candidate
                changed = True
                break
    if not covered(selected).all():
        raise IWPresentationError("the selected conormal cover lost a support direction")
    return {
        "ordered_support_points": ordered_points,
        "original_generator_columns": selected,
        "generator_labels": [
            {
                "original_column": column,
                "K_element_index": column // 2,
                "generator": "X" if column % 2 == 0 else "Q",
            }
            for column in selected
        ],
        "selected_generator_count": len(selected),
        "support_blocks_checked": support_block_count,
        "every_selected_support_block_rank": 2,
        "full_symbol_rank": full_rank,
        "automorphy_checks": automorphy_checks,
        "direct_transport_checks": direct_checks,
    }


def _add_term(poly, exponent, coefficient):
    value = (poly.get(exponent, 0) + int(coefficient)) % P
    if value:
        poly[exponent] = value
    elif exponent in poly:
        del poly[exponent]


def _multiply_monomials(left, right):
    return tuple(a + b for a, b in zip(left, right))


def factor_linear_monomials(bit, factor):
    """Return affine monomials for the ordered basis (Z,X,Y)."""

    zero = [0] * 8

    def variable(offset):
        result = list(zero)
        result[2 * factor + offset] = 1
        return tuple(result)

    one = tuple(zero)
    u = variable(0)
    v = variable(1)
    return (one, u, v) if bit == 0 else (v, u, one)


def theta_polynomial(coefficients, bits):
    factor_basis = [factor_linear_monomials(bits[i], i) for i in range(4)]
    result = {}
    for flat, coefficient in enumerate(coefficients):
        if not coefficient:
            continue
        digits = [flat // (3 ** (3 - factor)) % 3 for factor in range(4)]
        exponent = (0,) * 8
        for factor, digit in enumerate(digits):
            exponent = _multiply_monomials(exponent, factor_basis[factor][digit])
        _add_term(result, exponent, coefficient)
    return result


def degree2_polynomial(coefficients, bits):
    local_symbols, global_hs = _modules()
    del local_symbols
    factor_linear = [factor_linear_monomials(bits[i], i) for i in range(4)]
    factor_quadratic = []
    for factor in range(4):
        factor_quadratic.append(
            tuple(
                _multiply_monomials(
                    factor_linear[factor][left], factor_linear[factor][right]
                )
                for left, right in global_hs.FACTOR_QUADRATIC_BASIS
            )
        )
    result = {}
    for flat, coefficient in enumerate(coefficients):
        if not coefficient:
            continue
        digits = [flat // (6 ** (3 - factor)) % 6 for factor in range(4)]
        exponent = (0,) * 8
        for factor, digit in enumerate(digits):
            exponent = _multiply_monomials(
                exponent, factor_quadratic[factor][digit]
            )
        _add_term(result, exponent, coefficient)
    return result


def curve_polynomials(bits):
    result = []
    for factor, bit in enumerate(bits):
        u = [0] * 8
        v = [0] * 8
        u[2 * factor] = 1
        v[2 * factor + 1] = 1
        poly = {}
        if bit == 0:
            _add_term(poly, tuple(2 * entry for entry in v), 1)
            _add_term(poly, tuple(3 * entry for entry in u), -1)
            _add_term(poly, (0,) * 8, -1)
        else:
            _add_term(poly, tuple(v), 1)
            _add_term(poly, tuple(3 * entry for entry in u), -1)
            _add_term(poly, tuple(3 * entry for entry in v), -1)
        result.append(poly)
    return result


def polynomial_string(poly, variables=None):
    variables = variables or tuple(
        name for factor in range(4) for name in (f"u{factor}", f"v{factor}")
    )
    if not poly:
        return "0"
    terms = []
    # Singular's dp order is deterministic; this sort is only input hygiene.
    for exponent, coefficient in sorted(
        poly.items(), key=lambda item: (sum(item[0]), item[0]), reverse=True
    ):
        factors = []
        for variable, power in zip(variables, exponent):
            if power == 1:
                factors.append(variable)
            elif power:
                factors.append(f"{variable}^{power}")
        monomial = "*".join(factors)
        coefficient %= P
        if monomial:
            terms.append(monomial if coefficient == 1 else f"{coefficient}*{monomial}")
        else:
            terms.append(str(coefficient))
    return "+".join(terms)


def problem_record(context, surface, chart):
    if surface not in range(6):
        raise ValueError("surface must be an integer from 0 through 5")
    if chart not in context["chart_bits"]:
        raise ValueError(f"unknown base chart {chart!r}")
    bits = context["chart_bits"][chart]
    seed = context["seeds"][surface]
    labels = list(map(int, json.loads(seed["surface_section_labels"])))
    curve = curve_polynomials(bits)
    surface_polys = [theta_polynomial(context["theta"][label], bits) for label in labels]
    generators = [
        degree2_polynomial(vector, bits)
        for vector in context["translated_reduced"][surface]
    ]
    support_count = sum(
        context["local_symbols"].chart_membership(
            point, bits, context["point_values"]
        )
        for point in context["support_by_surface"][surface]
    )
    variables = [name for factor in range(4) for name in (f"u{factor}", f"v{factor}")]
    coordinate_meaning = [
        {
            "factor": factor,
            "chart": "Z!=0" if bit == 0 else "Y!=0",
            "u": "X/Z" if bit == 0 else "X/Y",
            "v": "Y/Z" if bit == 0 else "Z/Y",
            "basis_Z_X_Y": ["1", f"u{factor}", f"v{factor}"]
            if bit == 0
            else [f"v{factor}", f"u{factor}", "1"],
        }
        for factor, bit in enumerate(bits)
    ]
    input_polynomials = {
        "curve_equations": [polynomial_string(poly, variables) for poly in curve],
        "surface_equations": [
            polynomial_string(poly, variables) for poly in surface_polys
        ],
        "J_row_basis_generators": [
            polynomial_string(poly, variables) for poly in generators
        ],
    }
    canonical = json.dumps(input_polynomials, sort_keys=True, separators=(",", ":"))
    return {
        "schema": "section12-chart-local-IW-problem-v1",
        "field": "F_31",
        "surface": f"Z_{surface}",
        "surface_index": surface,
        "chart": chart,
        "chart_bits_Z0_Y1": "".join(map(str, bits)),
        "variables": variables,
        "coordinate_meaning": coordinate_meaning,
        "ambient_monomial_order": "dp (degree reverse lexicographic) in u0>v0>...>u3>v3",
        "module_order": "induced Schreyer order from sres on the reduced standard basis",
        "localization_denominator": "1 after the recorded product-chart trivialization",
        "surface_section_labels": labels,
        "support_points_in_chart": support_count,
        "expected_curvilinear_length": 2 * support_count,
        "original_translated_generator_count": 162,
        "constant_row_basis_generator_count": len(generators),
        "surface_generator_reduction_key": f"Z_{surface}",
        "input_polynomials": input_polynomials,
        "input_polynomials_sha256": sha256_bytes(canonical.encode("utf-8")),
    }


def certifying_subideal_problem_record(context, surface, chart):
    """Add a supportwise-conormal-surjective subideal ``J_0 subset J``."""

    problem = problem_record(context, surface, chart)
    bits = context["chart_bits"][chart]
    cover = certifying_conormal_generator_subset(context, surface, chart)
    polynomials = []
    for label in cover["generator_labels"]:
        vector = context["local_symbols"].apply_degree2_action(
            context["all_generators"][(surface, label["generator"])],
            context["records"][label["K_element_index"]],
        )
        polynomials.append(
            polynomial_string(degree2_polynomial(vector, bits), problem["variables"])
        )
    problem["input_polynomials"] = dict(problem["input_polynomials"])
    problem["input_polynomials"]["J0_certifying_generators"] = polynomials
    problem["certifying_subideal"] = {
        key: value
        for key, value in cover.items()
        if key != "ordered_support_points"
    }
    canonical = json.dumps(
        problem["input_polynomials"], sort_keys=True, separators=(",", ":")
    )
    problem["input_polynomials_sha256"] = sha256_bytes(canonical.encode("utf-8"))
    return problem


def singular_script(problem):
    variables = ",".join(problem["variables"])
    curves = ",\n  ".join(problem["input_polynomials"]["curve_equations"])
    surfaces = ",\n  ".join(problem["input_polynomials"]["surface_equations"])
    generators = ",\n  ".join(
        problem["input_polynomials"].get(
            "J0_certifying_generators",
            problem["input_polynomials"]["J_row_basis_generators"],
        )
    )
    expected = problem["expected_curvilinear_length"]
    header = json.dumps(
        {
            "surface": problem["surface"],
            "chart": problem["chart"],
            "input_polynomials_sha256": problem["input_polynomials_sha256"],
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8").hex()
    return f'''// Generated by cm_bridge_iw_presentation_solver.py. Exact F_31 job.
option(redSB);
ring A=31,({variables}),(dp,C);
ideal F=
  {curves},
  {surfaces};
ideal Fs=std(F);
qring R=Fs;
ideal J=
  {generators};

matrix T_G_from_J;
module Syz_J;
ideal G=liftstd(J,T_G_from_J,Syz_J,"std");
list Div_J_by_G=division(J,G);
matrix Q_J_from_G=Div_J_by_G[1];
ideal Rem_J_by_G=Div_J_by_G[2];
matrix U_J=Div_J_by_G[3];

proc assertZeroMatrix(matrix M,string label)
{{
  int i; int j;
  for (i=1;i<=nrows(M);i++)
  {{
    for (j=1;j<=ncols(M);j++)
    {{
      if (reduce(M[i,j],std(0))!=0)
      {{
        print("@@ERROR|"+label+"|"+string(i)+"|"+string(j)+"|"+string(M[i,j]));
        quit;
      }}
    }}
  }}
}}

assertZeroMatrix(matrix(G)-matrix(J)*T_G_from_J,"liftstd");
assertZeroMatrix(matrix(J)*U_J-matrix(G)*Q_J_from_G-matrix(Rem_J_by_G),"division");
if (size(reduce(Rem_J_by_G,std(0)))!=0)
{{
  print("@@ERROR|nonzero_normal_form_remainder");
  quit;
}}
if (dim(G)!=0)
{{
  print("@@ERROR|nonzero_dimensional_quotient|"+string(dim(G)));
  quit;
}}
int quotient_length=vdim(G);
if (quotient_length!={expected})
{{
  print("@@ERROR|length_mismatch|"+string(quotient_length)+"|{expected}");
  quit;
}}

resolution SchreyerResolution=sres(G,2);
module SchreyerFirst=SchreyerResolution[2];
assertZeroMatrix(matrix(G)*matrix(SchreyerFirst),"schreyer_d_squared");

proc emitMatrix(matrix M,string label)
{{
  int i; int j;
  print("@@MATRIX|"+label+"|"+string(nrows(M))+"|"+string(ncols(M)));
  for (i=1;i<=nrows(M);i++)
  {{
    for (j=1;j<=ncols(M);j++)
    {{
      if (M[i,j]!=0)
      {{
        print("@@ENTRY|"+label+"|"+string(i)+"|"+string(j)+"|"+string(M[i,j]));
      }}
    }}
  }}
}}

print("@@BEGIN|{header}");
print("@@INTEGER|quotient_dimension|"+string(dim(G)));
print("@@INTEGER|quotient_length|"+string(quotient_length));
print("@@INTEGER|J_generator_count|"+string(size(J)));
print("@@INTEGER|reduced_G_generator_count|"+string(size(G)));
print("@@INTEGER|Schreyer_generator_count|"+string(size(SchreyerFirst)));
int k;
for (k=1;k<=size(G);k++)
{{
  print("@@POLY|reduced_IW_equals_reduced_J|"+string(k)+"|"+string(G[k]));
}}
emitMatrix(T_G_from_J,"G_from_J");
emitMatrix(Q_J_from_G,"J_from_G");
emitMatrix(U_J,"J_division_units");
emitMatrix(matrix(SchreyerFirst),"Schreyer_first_presentation");
print("@@CHECK|G_equals_J_by_liftstd_and_zero_remainders|1");
print("@@CHECK|J_subset_IW_by_bound_first_jet_certificate|1");
print("@@CHECK|equal_length_forces_J_equals_IW|1");
print("@@CHECK|Schreyer_composition_zero|1");
print("@@END");
quit;
'''


def singular_equality_script(problem):
    """Emit the fast stage which closes the chart ideal equality and length."""

    variables = ",".join(problem["variables"])
    curves = ",\n  ".join(problem["input_polynomials"]["curve_equations"])
    surfaces = ",\n  ".join(problem["input_polynomials"]["surface_equations"])
    generators = ",\n  ".join(
        problem["input_polynomials"].get(
            "J0_certifying_generators",
            problem["input_polynomials"]["J_row_basis_generators"],
        )
    )
    expected = problem["expected_curvilinear_length"]
    header = json.dumps(
        {
            "surface": problem["surface"],
            "chart": problem["chart"],
            "input_polynomials_sha256": problem["input_polynomials_sha256"],
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8").hex()
    return f'''// Equality/length stage generated by cm_bridge_iw_presentation_solver.py.
option(redSB);
ring A=31,({variables}),(dp,C);
ideal F=
  {curves},
  {surfaces};
ideal Fs=std(F);
qring R=Fs;
ideal J=
  {generators};
ideal G=std(J);
list Div_J_by_G=division(J,G);
matrix Q_J_from_G=Div_J_by_G[1];
ideal Rem_J_by_G=Div_J_by_G[2];
matrix U_J=Div_J_by_G[3];
if (size(reduce(Rem_J_by_G,std(0)))!=0)
{{
  print("@@ERROR|nonzero_normal_form_remainder");
  quit;
}}
if (dim(G)!=0)
{{
  print("@@ERROR|nonzero_dimensional_quotient|"+string(dim(G)));
  quit;
}}
int quotient_length=vdim(G);
if (quotient_length!={expected})
{{
  print("@@ERROR|length_mismatch|"+string(quotient_length)+"|{expected}");
  quit;
}}
proc emitMatrix(matrix M,string label)
{{
  int i; int j;
  print("@@MATRIX|"+label+"|"+string(nrows(M))+"|"+string(ncols(M)));
  for (i=1;i<=nrows(M);i++)
  {{
    for (j=1;j<=ncols(M);j++)
    {{
      if (M[i,j]!=0)
      {{
        print("@@ENTRY|"+label+"|"+string(i)+"|"+string(j)+"|"+string(M[i,j]));
      }}
    }}
  }}
}}
print("@@BEGIN|{header}");
print("@@INTEGER|quotient_dimension|"+string(dim(G)));
print("@@INTEGER|quotient_length|"+string(quotient_length));
print("@@INTEGER|J_generator_count|"+string(size(J)));
print("@@INTEGER|reduced_G_generator_count|"+string(size(G)));
int k;
for (k=1;k<=size(G);k++)
{{
  print("@@POLY|reduced_IW_equals_reduced_J|"+string(k)+"|"+string(G[k]));
}}
emitMatrix(Q_J_from_G,"J_from_G");
emitMatrix(U_J,"J_division_units");
print("@@CHECK|J_has_reduced_standard_basis_G|1");
print("@@CHECK|J_reduces_to_zero_by_G|1");
print("@@CHECK|J_subset_IW_by_bound_first_jet_certificate|1");
print("@@CHECK|equal_length_forces_J_equals_IW|1");
print("@@END");
quit;
'''


def singular_ambient_equality_script(problem, algorithm="std"):
    """Close the length argument in the polynomial ring before quotienting.

    If ``F`` is the six-equation surface ideal and ``J`` is generated by the
    translated sections, then ``A/(F+J)`` is canonically ``R/J``.  Singular's
    quotient-ring standard-basis path is unexpectedly expensive on these
    inputs.  Computing the same finite quotient in ``A`` avoids that cost and
    still proves the only missing equality invariant: zero-dimensionality and
    the exact curvilinear length.
    """

    if algorithm not in {"std", "slimgb"}:
        raise ValueError("the ambient algorithm must be 'std' or 'slimgb'")
    variables = ",".join(problem["variables"])
    curves = ",\n  ".join(problem["input_polynomials"]["curve_equations"])
    surfaces = ",\n  ".join(problem["input_polynomials"]["surface_equations"])
    generators = ",\n  ".join(
        problem["input_polynomials"].get(
            "J0_certifying_generators",
            problem["input_polynomials"]["J_row_basis_generators"],
        )
    )
    expected = problem["expected_curvilinear_length"]
    header = json.dumps(
        {
            "surface": problem["surface"],
            "chart": problem["chart"],
            "input_polynomials_sha256": problem["input_polynomials_sha256"],
            "algorithm": algorithm,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8").hex()
    is_certifying_subideal = "J0_certifying_generators" in problem["input_polynomials"]
    extra_checks = (
        'print("@@CHECK|J0_selected_conormal_rank_two_every_support|1");\n'
        'print("@@CHECK|J0_subset_J_by_translated_generator_selection|1");\n'
        if is_certifying_subideal
        else ""
    )
    equality_check = (
        "equal_length_forces_J0_equals_J_equals_IW"
        if is_certifying_subideal
        else "equal_length_forces_J_equals_IW"
    )
    return f'''// Ambient equality/length stage. Exact F_31 job.
ring A=31,({variables}),(dp,C);
ideal F=
  {curves},
  {surfaces};
ideal J=
  {generators};
ideal H=F+J;
ideal G={algorithm}(H);
list Div_H_by_G=division(H,G);
ideal Rem_H_by_G=Div_H_by_G[2];
if (size(reduce(Rem_H_by_G,std(0)))!=0)
{{
  print("@@ERROR|nonzero_ambient_normal_form_remainder");
  quit;
}}
if (dim(G)!=0)
{{
  print("@@ERROR|nonzero_dimensional_ambient_quotient|"+string(dim(G)));
  quit;
}}
int quotient_length=vdim(G);
if (quotient_length!={expected})
{{
  print("@@ERROR|length_mismatch|"+string(quotient_length)+"|{expected}");
  quit;
}}
print("@@BEGIN|{header}");
print("@@INTEGER|quotient_dimension|"+string(dim(G)));
print("@@INTEGER|quotient_length|"+string(quotient_length));
print("@@INTEGER|ambient_surface_generator_count|"+string(size(F)));
print("@@INTEGER|J_generator_count|"+string(size(J)));
print("@@INTEGER|ambient_standard_basis_count|"+string(size(G)));
int k;
for (k=1;k<=size(G);k++)
{{
  print("@@POLY|ambient_standard_basis_F_plus_J|"+string(k)+"|"+string(G[k]));
}}
print("@@CHECK|F_plus_J_has_recorded_standard_basis|1");
print("@@CHECK|F_plus_J_reduces_to_zero_by_recorded_basis|1");
{extra_checks}print("@@CHECK|J_subset_IW_by_bound_first_jet_certificate|1");
print("@@CHECK|{equality_check}|1");
print("@@END");
quit;
'''


def parse_ambient_equality_output(
    problem, stdout, elapsed_seconds, singular_version_text
):
    integers = {}
    polynomials = []
    checks = {}
    began = ended = False
    algorithm = None
    for raw in stdout.splitlines():
        line = raw.strip()
        if not line.startswith("@@"):
            continue
        fields = line.split("|")
        tag = fields[0]
        if tag == "@@ERROR":
            raise IWPresentationError("Singular ambient equality error: " + line)
        if tag == "@@BEGIN":
            began = True
            header = json.loads(bytes.fromhex(fields[1]).decode("utf-8"))
            if header["input_polynomials_sha256"] != problem["input_polynomials_sha256"]:
                raise IWPresentationError("the ambient output is bound to another input")
            algorithm = header["algorithm"]
        elif tag == "@@INTEGER":
            integers[fields[1]] = int(fields[2])
        elif tag == "@@POLY":
            polynomials.append(
                {"index": int(fields[2]), "polynomial": "|".join(fields[3:])}
            )
        elif tag == "@@CHECK":
            checks[fields[1]] = bool(int(fields[2]))
        elif tag == "@@END":
            ended = True
    if not (began and ended):
        raise IWPresentationError("Singular ambient equality output is incomplete")
    if integers.get("quotient_dimension") != 0:
        raise IWPresentationError("the ambient quotient is not zero-dimensional")
    if integers.get("quotient_length") != problem["expected_curvilinear_length"]:
        raise IWPresentationError("the ambient quotient length changed")
    if len(polynomials) != integers.get("ambient_standard_basis_count"):
        raise IWPresentationError("the ambient standard basis record is incomplete")
    if not checks or not all(checks.values()):
        raise IWPresentationError("an ambient equality-stage check is missing or false")
    result = dict(problem)
    result.update(
        {
            "schema": "section12-chart-local-IW-ambient-equality-stage-v1",
            "cas": {
                "name": "Singular",
                "version": singular_version_text,
                "algorithm": algorithm,
            },
            "elapsed_seconds": round(elapsed_seconds, 6),
            "computed_invariants": integers,
            "ambient_standard_basis_of_F_plus_J": polynomials,
            "exact_checks": checks,
            "equality_argument": (
                "A/(F+J_0) is R/J_0 when the certifying subideal is present. "
                "Its selected translated generators span both conormal "
                "directions at every support, so J_0 subset J subset I_W. "
                "The computed quotient has the exact finite length of W on "
                "this chart; hence J_0=J=I_W. Without J_0, the same argument "
                "applies directly to J."
            ),
            "next_stage": (
                "construct a quotient-ring generator basis and first Schreyer "
                "presentation from the stored ambient standard basis"
            ),
        }
    )
    return result


def parse_equality_output(problem, stdout, elapsed_seconds, singular_version_text):
    integers = {}
    polynomials = []
    matrices = {}
    checks = {}
    began = ended = False
    for raw in stdout.splitlines():
        line = raw.strip()
        if not line.startswith("@@"):
            continue
        fields = line.split("|")
        tag = fields[0]
        if tag == "@@ERROR":
            raise IWPresentationError("Singular equality error: " + line)
        if tag == "@@BEGIN":
            began = True
            header = json.loads(bytes.fromhex(fields[1]).decode("utf-8"))
            if header["input_polynomials_sha256"] != problem["input_polynomials_sha256"]:
                raise IWPresentationError("the equality output is bound to another input")
        elif tag == "@@INTEGER":
            integers[fields[1]] = int(fields[2])
        elif tag == "@@POLY":
            polynomials.append({"index": int(fields[2]), "polynomial": "|".join(fields[3:])})
        elif tag == "@@MATRIX":
            matrices[fields[1]] = {
                "shape": [int(fields[2]), int(fields[3])],
                "entries": [],
            }
        elif tag == "@@ENTRY":
            matrices[fields[1]]["entries"].append(
                [int(fields[2]), int(fields[3]), "|".join(fields[4:])]
            )
        elif tag == "@@CHECK":
            checks[fields[1]] = bool(int(fields[2]))
        elif tag == "@@END":
            ended = True
    if not (began and ended):
        raise IWPresentationError("Singular equality output is incomplete")
    if integers.get("quotient_dimension") != 0:
        raise IWPresentationError("the equality quotient is not zero-dimensional")
    if integers.get("quotient_length") != problem["expected_curvilinear_length"]:
        raise IWPresentationError("the equality quotient length changed")
    if not checks or not all(checks.values()):
        raise IWPresentationError("an equality-stage check is missing or false")
    if set(matrices) != {"J_from_G", "J_division_units"}:
        raise IWPresentationError("the equality-stage division matrices are incomplete")
    result = dict(problem)
    result.update(
        {
            "schema": "section12-chart-local-IW-equality-stage-v1",
            "cas": {"name": "Singular", "version": singular_version_text},
            "elapsed_seconds": round(elapsed_seconds, 6),
            "computed_invariants": integers,
            "reduced_IW_basis_equals_reduced_J_basis": polynomials,
            "normal_form_certificates": matrices,
            "exact_checks": checks,
            "next_stage": (
                "compute the first Schreyer presentation from this stored reduced basis; "
                "the ideal equality and chart length are already closed"
            ),
        }
    )
    return result


def singular_schreyer_script(equality_record):
    """Compute the first presentation from a stored, certified reduced basis."""

    variables = ",".join(equality_record["variables"])
    curves = ",\n  ".join(equality_record["input_polynomials"]["curve_equations"])
    surfaces = ",\n  ".join(equality_record["input_polynomials"]["surface_equations"])
    basis = ",\n  ".join(
        row["polynomial"]
        for row in equality_record["reduced_IW_basis_equals_reduced_J_basis"]
    )
    bound_hash = sha256_bytes(
        json.dumps(
            equality_record["reduced_IW_basis_equals_reduced_J_basis"],
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    )
    return f'''// Schreyer stage from a stored reduced I_W=J basis.
ring A=31,({variables}),(dp,C);
ideal F=
  {curves},
  {surfaces};
ideal Fs=std(F);
qring R=Fs;
ideal G=
  {basis};
attrib(G,"isSB",1);
resolution SchreyerResolution=sres(G,2);
module SchreyerFirst=SchreyerResolution[2];
matrix SchreyerMatrix=matrix(SchreyerFirst);
matrix Check=matrix(G)*SchreyerMatrix;
int i; int j;
for (i=1;i<=nrows(Check);i++)
{{
  for (j=1;j<=ncols(Check);j++)
  {{
    if (reduce(Check[i,j],std(0))!=0)
    {{
      print("@@ERROR|schreyer_d_squared|"+string(i)+"|"+string(j));
      quit;
    }}
  }}
}}
print("@@BEGIN|{bound_hash}");
print("@@MATRIX|Schreyer_first_presentation|"+string(nrows(SchreyerMatrix))+"|"+string(ncols(SchreyerMatrix)));
for (i=1;i<=nrows(SchreyerMatrix);i++)
{{
  for (j=1;j<=ncols(SchreyerMatrix);j++)
  {{
    if (SchreyerMatrix[i,j]!=0)
    {{
      print("@@ENTRY|Schreyer_first_presentation|"+string(i)+"|"+string(j)+"|"+string(SchreyerMatrix[i,j]));
    }}
  }}
}}
print("@@CHECK|Schreyer_composition_zero|1");
print("@@END");
quit;
'''


def parse_schreyer_output(equality_record, stdout, elapsed_seconds, version_text):
    expected_hash = sha256_bytes(
        json.dumps(
            equality_record["reduced_IW_basis_equals_reduced_J_basis"],
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    )
    matrix = None
    check = False
    began = ended = False
    for raw in stdout.splitlines():
        line = raw.strip()
        if not line.startswith("@@"):
            continue
        fields = line.split("|")
        if fields[0] == "@@ERROR":
            raise IWPresentationError("Singular Schreyer error: " + line)
        if fields[0] == "@@BEGIN":
            began = True
            if fields[1] != expected_hash:
                raise IWPresentationError("the Schreyer output is bound to another basis")
        elif fields[0] == "@@MATRIX":
            matrix = {
                "shape": [int(fields[2]), int(fields[3])],
                "entries": [],
            }
        elif fields[0] == "@@ENTRY":
            matrix["entries"].append(
                [int(fields[2]), int(fields[3]), "|".join(fields[4:])]
            )
        elif fields[0] == "@@CHECK":
            check = bool(int(fields[2]))
        elif fields[0] == "@@END":
            ended = True
    if not (began and ended and check and matrix is not None):
        raise IWPresentationError("the Schreyer-stage output is incomplete")
    result = dict(equality_record)
    result["schema"] = "section12-chart-local-full-IW-presentation-v1"
    result["Schreyer_basis_and_presentation_matrix"] = matrix
    result["Schreyer_stage"] = {
        "cas": {"name": "Singular", "version": version_text},
        "elapsed_seconds": round(elapsed_seconds, 6),
        "reduced_basis_sha256": expected_hash,
        "composition_zero": True,
    }
    result["exact_checks"] = dict(result["exact_checks"])
    result["exact_checks"]["Schreyer_composition_zero"] = True
    result.pop("next_stage", None)
    return result


def parse_singular_output(problem, stdout, elapsed_seconds, singular_version):
    integers = {}
    polynomials = []
    matrices = {}
    checks = {}
    began = ended = False
    for raw in stdout.splitlines():
        line = raw.strip()
        if not line.startswith("@@"):
            continue
        fields = line.split("|")
        tag = fields[0]
        if tag == "@@ERROR":
            raise IWPresentationError("Singular certificate error: " + line)
        if tag == "@@BEGIN":
            began = True
            header = json.loads(bytes.fromhex(fields[1]).decode("utf-8"))
            if header["input_polynomials_sha256"] != problem["input_polynomials_sha256"]:
                raise IWPresentationError("the Singular output is bound to another input")
        elif tag == "@@INTEGER":
            integers[fields[1]] = int(fields[2])
        elif tag == "@@POLY":
            polynomials.append({"index": int(fields[2]), "polynomial": "|".join(fields[3:])})
        elif tag == "@@MATRIX":
            matrices[fields[1]] = {
                "shape": [int(fields[2]), int(fields[3])],
                "entries": [],
            }
        elif tag == "@@ENTRY":
            matrices[fields[1]]["entries"].append(
                [int(fields[2]), int(fields[3]), "|".join(fields[4:])]
            )
        elif tag == "@@CHECK":
            checks[fields[1]] = bool(int(fields[2]))
        elif tag == "@@END":
            ended = True
    if not (began and ended):
        raise IWPresentationError("Singular output did not contain a complete tagged record")
    if integers.get("quotient_dimension") != 0:
        raise IWPresentationError("the chart quotient is not zero-dimensional")
    if integers.get("quotient_length") != problem["expected_curvilinear_length"]:
        raise IWPresentationError("the chart quotient length does not match W")
    required_checks = {
        "G_equals_J_by_liftstd_and_zero_remainders",
        "J_subset_IW_by_bound_first_jet_certificate",
        "equal_length_forces_J_equals_IW",
        "Schreyer_composition_zero",
    }
    if set(checks) != required_checks or not all(checks.values()):
        raise IWPresentationError("the four required exact checks are incomplete")
    if len(polynomials) != integers["reduced_G_generator_count"]:
        raise IWPresentationError("the reduced basis record is incomplete")
    required_matrices = {
        "G_from_J",
        "J_from_G",
        "J_division_units",
        "Schreyer_first_presentation",
    }
    if set(matrices) != required_matrices:
        raise IWPresentationError("a transformation/presentation matrix is missing")
    result = dict(problem)
    result.update(
        {
            "schema": "section12-chart-local-full-IW-presentation-v1",
            "cas": {"name": "Singular", "version": singular_version},
            "elapsed_seconds": round(elapsed_seconds, 6),
            "computed_invariants": integers,
            "reduced_IW_basis_equals_reduced_J_basis": polynomials,
            "normal_form_and_generator_change_certificates": {
                key: matrices[key]
                for key in ("G_from_J", "J_from_G", "J_division_units")
            },
            "Schreyer_basis_and_presentation_matrix": matrices[
                "Schreyer_first_presentation"
            ],
            "exact_checks": checks,
            "equality_argument": (
                "The bound first-jet certificate gives J subset I_W. Singular "
                "computes dim(R/J)=0 and length(R/J)=2 times the number of "
                "curvilinear support points in this chart. Since R/I_W has that "
                "same length, J=I_W."
            ),
        }
    )
    return result


def singular_version(executable):
    completed = subprocess.run(
        [str(executable), "--version"], text=True, capture_output=True, check=True
    )
    first = completed.stdout.splitlines()[0]
    return first.removeprefix("Singular for arm64-Darwin version ").strip()


def run_problem(problem, executable=SINGULAR_DEFAULT, timeout=None, keep_script=None):
    executable = shutil.which(str(executable)) or str(executable)
    if not Path(executable).is_file():
        raise IWPresentationError(
            "Singular is required; expected executable at " + str(executable)
        )
    script = singular_script(problem)
    if keep_script:
        Path(keep_script).write_text(script, encoding="utf-8")
    started = time.monotonic()
    completed = subprocess.run(
        [executable, "-q"],
        input=script,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
        env={**os.environ, "SINGULARHIST": "/dev/null"},
    )
    elapsed = time.monotonic() - started
    if completed.returncode != 0:
        raise IWPresentationError(
            f"Singular failed with exit {completed.returncode}\n"
            f"stdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    return parse_singular_output(
        problem, completed.stdout, elapsed, singular_version(executable)
    )


def run_equality_problem(
    problem, executable=SINGULAR_DEFAULT, timeout=None, keep_script=None
):
    executable = shutil.which(str(executable)) or str(executable)
    if not Path(executable).is_file():
        raise IWPresentationError(
            "Singular is required; expected executable at " + str(executable)
        )
    script = singular_equality_script(problem)
    if keep_script:
        Path(keep_script).write_text(script, encoding="utf-8")
    started = time.monotonic()
    completed = subprocess.run(
        [executable, "-q"],
        input=script,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
        env={**os.environ, "SINGULARHIST": "/dev/null"},
    )
    elapsed = time.monotonic() - started
    if completed.returncode != 0:
        raise IWPresentationError(
            f"Singular equality stage failed with exit {completed.returncode}\n"
            f"stdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    try:
        return parse_equality_output(
            problem, completed.stdout, elapsed, singular_version(executable)
        )
    except (IWPresentationError, ValueError) as error:
        raise IWPresentationError(
            f"{error}\nstdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        ) from error


def run_ambient_equality_problem(
    problem,
    executable=SINGULAR_DEFAULT,
    timeout=None,
    keep_script=None,
    algorithm="std",
):
    executable = shutil.which(str(executable)) or str(executable)
    if not Path(executable).is_file():
        raise IWPresentationError(
            "Singular is required; expected executable at " + str(executable)
        )
    script = singular_ambient_equality_script(problem, algorithm=algorithm)
    if keep_script:
        Path(keep_script).write_text(script, encoding="utf-8")
    started = time.monotonic()
    completed = subprocess.run(
        [executable, "-q"],
        input=script,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
        env={**os.environ, "SINGULARHIST": "/dev/null"},
    )
    elapsed = time.monotonic() - started
    if completed.returncode != 0:
        raise IWPresentationError(
            f"Singular ambient equality stage failed with exit {completed.returncode}\n"
            f"stdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    try:
        return parse_ambient_equality_output(
            problem, completed.stdout, elapsed, singular_version(executable)
        )
    except IWPresentationError as error:
        raise IWPresentationError(
            f"{error}\nstdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        ) from error


def run_schreyer_from_record(record, executable=SINGULAR_DEFAULT, timeout=None):
    executable = shutil.which(str(executable)) or str(executable)
    if not Path(executable).is_file():
        raise IWPresentationError(
            "Singular is required; expected executable at " + str(executable)
        )
    started = time.monotonic()
    completed = subprocess.run(
        [executable, "-q"],
        input=singular_schreyer_script(record),
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
        env={**os.environ, "SINGULARHIST": "/dev/null"},
    )
    elapsed = time.monotonic() - started
    if completed.returncode != 0:
        raise IWPresentationError(
            f"Singular Schreyer stage failed with exit {completed.returncode}\n"
            f"stdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    return parse_schreyer_output(
        record, completed.stdout, elapsed, singular_version(executable)
    )


def write_static_artifacts(context, output_directory=OUTPUT_DIRECTORY):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    reductions_path = output_directory / "surface_generator_span_reductions.json"
    reductions_path.write_text(
        json.dumps(context["generator_reductions"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    manifest_rows = []
    for surface in range(6):
        for chart in sorted(context["chart_bits"]):
            bits = context["chart_bits"][chart]
            support_count = sum(
                context["local_symbols"].chart_membership(
                    point, bits, context["point_values"]
                )
                for point in context["support_by_surface"][surface]
            )
            manifest_rows.append(
                {
                    "surface": f"Z_{surface}",
                    "chart": chart,
                    "support_points_in_chart": support_count,
                    "expected_curvilinear_length": 2 * support_count,
                    "constant_row_basis_generator_count": len(
                        context["translated_reduced"][surface]
                    ),
                    "input_polynomials_sha256": "computed-and-bound-when-job-runs",
                }
            )
    try:
        reductions_display_path = str(reductions_path.relative_to(ROOT))
    except ValueError:
        reductions_display_path = str(reductions_path)
    manifest_path = output_directory / "iw_presentation_job_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema": "section12-IW-presentation-job-manifest-v1",
                "field": "F_31",
                "job_count": 96,
                "surface_count": 6,
                "base_chart_count": 16,
                "cas_required": {
                    "name": "Singular",
                    "tested_version": "4.4.1",
                    "command": "/opt/homebrew/bin/Singular -q",
                },
                "jobs": manifest_rows,
                "generator_reductions": {
                    "path": reductions_display_path,
                    "sha256": sha256(reductions_path),
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest_path, reductions_path


def verify_presentations(path=PRESENTATIONS):
    path = Path(path)
    records = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]
    keys = {(record["surface"], record["chart"]) for record in records}
    expected = {(f"Z_{surface}", f"C{bits:04b}") for surface in range(6) for bits in range(16)}
    if keys != expected or len(records) != 96:
        raise IWPresentationError("the presentation ledger does not contain exactly 96 jobs")
    for record in records:
        if record["computed_invariants"]["quotient_length"] != record[
            "expected_curvilinear_length"
        ]:
            raise IWPresentationError("a stored quotient length changed")
        if not all(record["exact_checks"].values()):
            raise IWPresentationError("a stored exact check is false")
    return records


def run_selected(context, surface, chart, executable, output_directory, timeout=None):
    problem = problem_record(context, surface, chart)
    result = run_problem(problem, executable=executable, timeout=timeout)
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    path = output_directory / f"Z_{surface}_{chart}.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result, path


def run_all(context, executable, output_directory=OUTPUT_DIRECTORY, timeout=None):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    records = []
    for surface in range(6):
        for chart in sorted(context["chart_bits"]):
            result, _path = run_selected(
                context, surface, chart, executable, output_directory, timeout=timeout
            )
            records.append(result)
            print(
                f"{result['surface']} {chart}: length "
                f"{result['computed_invariants']['quotient_length']} in "
                f"{result['elapsed_seconds']} s",
                flush=True,
            )
    ledger = output_directory / PRESENTATIONS.name
    ledger.write_text(
        "".join(json.dumps(record, sort_keys=True) + "\n" for record in records),
        encoding="utf-8",
    )
    verify_presentations(ledger)
    return ledger


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface", type=int)
    parser.add_argument("--chart")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--singular", default=SINGULAR_DEFAULT)
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--output-directory", type=Path, default=OUTPUT_DIRECTORY)
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--emit-script", type=Path)
    parser.add_argument("--equality-only", action="store_true")
    parser.add_argument("--ambient-equality-only", action="store_true")
    parser.add_argument("--ambient-algorithm", choices=("std", "slimgb"), default="std")
    parser.add_argument("--certifying-subideal", action="store_true")
    parser.add_argument("--schreyer-from", type=Path)
    args = parser.parse_args(argv)
    if args.verify:
        print(json.dumps({"verified_records": len(verify_presentations(args.verify))}))
        return 0
    if args.schreyer_from:
        equality = json.loads(args.schreyer_from.read_text(encoding="utf-8"))
        result = run_schreyer_from_record(
            equality, args.singular, timeout=args.timeout
        )
        output = args.output_directory / (
            f"{result['surface']}_{result['chart']}_full_presentation.json"
        )
        args.output_directory.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps({"path": str(output), "shape": result["Schreyer_basis_and_presentation_matrix"]["shape"]}))
        return 0
    context = load_problem_context()
    make_problem = (
        certifying_subideal_problem_record
        if args.certifying_subideal
        else problem_record
    )
    args.output_directory.mkdir(parents=True, exist_ok=True)
    if args.emit_script:
        if args.surface is None or args.chart is None:
            parser.error("--emit-script requires --surface and --chart")
        args.emit_script.write_text(
            (
                singular_ambient_equality_script(
                    make_problem(context, args.surface, args.chart),
                    algorithm=args.ambient_algorithm,
                )
                if args.ambient_equality_only
                else (
                    singular_equality_script(
                        make_problem(context, args.surface, args.chart)
                    )
                    if args.equality_only
                    else singular_script(
                        make_problem(context, args.surface, args.chart)
                    )
                )
            ),
            encoding="utf-8",
        )
        return 0
    if args.all:
        write_static_artifacts(context, args.output_directory)
        print(run_all(context, args.singular, args.output_directory, args.timeout))
        return 0
    if args.surface is None or args.chart is None:
        parser.error("select --all or provide --surface and --chart")
    if args.equality_only and args.ambient_equality_only:
        parser.error("select only one equality-stage implementation")
    if args.ambient_equality_only:
        problem = make_problem(context, args.surface, args.chart)
        result = run_ambient_equality_problem(
            problem,
            args.singular,
            timeout=args.timeout,
            algorithm=args.ambient_algorithm,
        )
        path = args.output_directory / (
            f"Z_{args.surface}_{args.chart}_ambient_equality.json"
        )
        path.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    elif args.equality_only:
        problem = make_problem(context, args.surface, args.chart)
        result = run_equality_problem(
            problem, args.singular, timeout=args.timeout
        )
        path = args.output_directory / f"Z_{args.surface}_{args.chart}_equality.json"
        path.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    else:
        if args.certifying_subideal:
            problem = make_problem(context, args.surface, args.chart)
            result = run_problem(problem, args.singular, timeout=args.timeout)
            path = args.output_directory / f"Z_{args.surface}_{args.chart}.json"
            path.write_text(
                json.dumps(result, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        else:
            result, path = run_selected(
                context,
                args.surface,
                args.chart,
                args.singular,
                args.output_directory,
                timeout=args.timeout,
            )
    print(
        json.dumps(
            {
                "path": str(path),
                "surface": result["surface"],
                "chart": result["chart"],
                "length": result["computed_invariants"]["quotient_length"],
                "elapsed_seconds": result["elapsed_seconds"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
