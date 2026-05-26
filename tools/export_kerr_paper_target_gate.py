#!/usr/bin/env python3
"""Export the strict Kerr paper-target gate.

This exporter is deliberately stricter than the current kerr_magnetosphere
validator.  The repository problem is a linear surrogate; the paper target is a
nonlinear Kerr force-free / Grad-Shafranov target motivated by the literature.
Until that full target validator exists, this script must be able to say
``no_candidate_yet`` instead of upgrading surrogate rows into paper claims.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sqlite3
import subprocess
import sys
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
DEFAULT_JSON = DOCS / "kerr-paper-target-gate.json"
DEFAULT_MD = DOCS / "kerr-paper-target-gate.md"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


PAPER_TARGET_VALIDATOR_IMPLEMENTED = True
FULL_TARGET_REJECTION = "full nonlinear Kerr split-monopole GSE residual is not exact zero"
REQUIRED_VARIABLES = {"r", "x", "a"}

DEFAULT_PROBES = [
    "1 - x",
    "x",
    "1/(1 - 1)",
    "1 - x + a**2*r*x",
]

CORRECTION_COEFFICIENTS = [
    sp.Integer(-1),
    sp.Rational(-1, 2),
    sp.Rational(1, 2),
    sp.Integer(1),
]

CORRECTION_ANGULAR_FACTORS = [
    "1",
    "x",
    "x**2",
    "1-x**2",
    "x*(1-x**2)",
    "(1-x**2)**2",
    "x**2*(1-x**2)",
    "x*(1-x**2)**2",
]

CORRECTION_RADIAL_FACTORS = [
    "1/r",
    "1/r**2",
    "1/r**3",
    "1/(r - 2*M)",
    "1/(r - 2*M)**2",
    "1/(r*(r - 2*M))",
]

PAIR_CORRECTION_COEFFICIENTS = [
    sp.Integer(-1),
    sp.Integer(1),
]

PAIR_CORRECTION_ANGULAR_FACTORS = [
    "x",
    "1-x**2",
    "x*(1-x**2)",
    "(1-x**2)**2",
]

PAIR_CORRECTION_RADIAL_FACTORS = [
    "1/r",
    "1/r**2",
    "1/(r - 2*M)",
]

COEFFICIENT_SOLVE_SERIES_ORDER = 4

LITERATURE_SOURCES = [
    {
        "id": "mahlmann_cerda_duran_aloy_2018_kerr_gse_numerics",
        "citation": (
            "J. F. Mahlmann, P. Cerda-Duran, M. A. Aloy et al., "
            "Numerically solving the relativistic Grad-Shafranov equation "
            "in Kerr spacetimes: Numerical techniques, MNRAS 477, 3927-3946 (2018)"
        ),
        "arxiv": "1802.00815",
        "arxiv_url": "https://arxiv.org/abs/1802.00815",
        "doi": "10.1093/mnras/sty858",
        "doi_url": "https://doi.org/10.1093/mnras/sty858",
        "role": "establishes the Kerr force-free Grad-Shafranov numerical target family",
    },
    {
        "id": "camilloni_grignani_harmark_oliveri_orselli_2020_extreme_kerr",
        "citation": (
            "F. Camilloni, G. Grignani, T. Harmark, R. Oliveri, M. Orselli, "
            "Moving away from the Near-Horizon Attractor of the Extreme Kerr "
            "Force-Free Magnetosphere, JCAP 10, 048 (2020)"
        ),
        "arxiv": "2007.15665",
        "arxiv_url": "https://arxiv.org/abs/2007.15665",
        "doi": "10.1088/1475-7516/2020/10/048",
        "doi_url": "https://doi.org/10.1088/1475-7516/2020/10/048",
        "role": (
            "states the no-known-exact-analytic stationary, axisymmetric, "
            "magnetically dominated extreme-Kerr force-free solution target"
        ),
    },
]

PAPER_TARGET = {
    "name": "Kerr force-free magnetosphere / relativistic Grad-Shafranov paper target",
    "status": "closed_split_monopole_residual_gate_implemented",
    "validator_implemented": PAPER_TARGET_VALIDATOR_IMPLEMENTED,
    "claim": (
        "Rows are paper candidates only if they pass the nonlinear Kerr "
        "force-free Grad-Shafranov residual with the fixed split-monopole "
        "potential functions plus the domain, anchor, and singularity gates."
    ),
    "current_repo_problem": "kerr_magnetosphere",
    "current_repo_problem_boundary": (
        "linear surrogate used only to generate a bounded expression table; "
        "paper admission is decided by the full split-monopole residual gate"
    ),
}

SAFE_POINTS = [
    {"M": sp.Integer(1), "a": sp.Rational(1, 5), "r": sp.Rational(5, 2), "x": sp.Rational(1, 3)},
    {"M": sp.Integer(1), "a": sp.Rational(1, 2), "r": sp.Integer(3), "x": -sp.Rational(2, 5)},
    {"M": sp.Integer(1), "a": sp.Rational(4, 5), "r": sp.Integer(4), "x": sp.Rational(1, 5)},
]


@dataclass(frozen=True)
class DbRun:
    db_path: Path
    table_name: str
    run_id: str
    command: list[str] | None = None
    max_depth: int | None = None
    validators: int | None = None
    wall_timeout_s: float | None = None
    validation_timeout_s: float | None = None


@dataclass(frozen=True)
class CandidateRow:
    row_id: int | None
    expression: str
    depth: int | None = None
    validation_reason: str | None = None
    source: str = "engine"


def correction_basis() -> list[tuple[str, sp.Basic]]:
    """Bounded finite-spin corrections around the split-monopole anchor."""
    return _cartesian_basis(CORRECTION_ANGULAR_FACTORS, CORRECTION_RADIAL_FACTORS)


def pair_correction_basis() -> list[tuple[str, sp.Basic]]:
    """Smaller basis used for bounded two-term finite-spin corrections."""
    return _cartesian_basis(PAIR_CORRECTION_ANGULAR_FACTORS, PAIR_CORRECTION_RADIAL_FACTORS)


def _cartesian_basis(
    angular_factors: list[str],
    radial_factors: list[str],
) -> list[tuple[str, sp.Basic]]:
    locals_map = sympify_locals()
    basis: list[tuple[str, sp.Basic]] = []
    for angular_name in angular_factors:
        angular_expr = sp.sympify(angular_name, locals=locals_map)
        for radial_name in radial_factors:
            radial_expr = sp.sympify(radial_name, locals=locals_map)
            basis.append(
                (
                    f"({angular_name})*({radial_name})",
                    sp.factor(angular_expr * radial_expr),
                )
            )
    return basis


def generate_anchor_correction_rows() -> list[CandidateRow]:
    locals_map = sympify_locals()
    a = locals_map["a"]
    x = locals_map["x"]
    rows: list[CandidateRow] = []
    seen: set[str] = set()
    for basis_name, basis_expr in correction_basis():
        for coeff in CORRECTION_COEFFICIENTS:
            expr = sp.factor(1 - x + a**2 * coeff * basis_expr)
            expr_str = sp.sstr(expr)
            if expr_str in seen:
                continue
            seen.add(expr_str)
            rows.append(
                CandidateRow(
                    row_id=None,
                    expression=expr_str,
                    depth=None,
                    validation_reason=(
                        f"expanded anchor correction grammar: coeff={sp.sstr(coeff)}, "
                        f"basis={basis_name}"
                    ),
                    source="anchor_correction_grammar",
                )
            )
    return rows


def generate_anchor_pair_correction_rows() -> list[CandidateRow]:
    locals_map = sympify_locals()
    a = locals_map["a"]
    x = locals_map["x"]
    rows: list[CandidateRow] = []
    seen: set[str] = set()
    basis = pair_correction_basis()
    for (left_name, left_expr), (right_name, right_expr) in combinations(basis, 2):
        for left_coeff, right_coeff in product(
            PAIR_CORRECTION_COEFFICIENTS,
            PAIR_CORRECTION_COEFFICIENTS,
        ):
            correction = left_coeff * left_expr + right_coeff * right_expr
            expr = sp.factor(1 - x + a**2 * correction)
            expr_str = sp.sstr(expr)
            if expr_str in seen:
                continue
            seen.add(expr_str)
            rows.append(
                CandidateRow(
                    row_id=None,
                    expression=expr_str,
                    depth=None,
                    validation_reason=(
                        "two-term anchor correction grammar: "
                        f"left_coeff={sp.sstr(left_coeff)}, left_basis={left_name}, "
                        f"right_coeff={sp.sstr(right_coeff)}, right_basis={right_name}"
                    ),
                    source="anchor_pair_correction_grammar",
                )
            )
    return rows


def solve_leading_order_coefficient_rows() -> tuple[list[CandidateRow], dict[str, Any]]:
    locals_map = sympify_locals()
    a = locals_map["a"]
    x = locals_map["x"]
    m = locals_map["M"]
    r = locals_map["r"]
    basis = correction_basis()
    coeff_symbols = sp.symbols(f"c0:{len(basis)}")
    correction = sum(coeff * basis_expr for coeff, (_, basis_expr) in zip(coeff_symbols, basis))
    psi = 1 - x + a**2 * correction
    residual = full_kerr_split_monopole_residual(psi)
    leading_residual = (
        sp.series(residual.subs(m, 1), a, 0, COEFFICIENT_SOLVE_SERIES_ORDER)
        .removeO()
        .coeff(a, 2)
    )
    numerator = sp.factor(sp.together(leading_residual).as_numer_denom()[0])
    polynomial = sp.Poly(numerator, r, x)
    equations = [sp.expand(coeff) for coeff in polynomial.coeffs()]
    matrix, rhs = sp.linear_eq_to_matrix(equations, coeff_symbols)
    solution_set = sp.linsolve((matrix, rhs), coeff_symbols)

    exact_solutions: list[tuple[sp.Basic, ...]] = []
    skipped_parametric = 0
    if solution_set != sp.EmptySet:
        for solution in solution_set:
            if any(value.free_symbols for value in solution):
                skipped_parametric += 1
                continue
            exact_solutions.append(tuple(solution))

    rows: list[CandidateRow] = []
    for solution in exact_solutions:
        solved_correction = sum(
            value * basis_expr for value, (_, basis_expr) in zip(solution, basis)
        )
        expr = sp.factor(1 - x + a**2 * solved_correction)
        rows.append(
            CandidateRow(
                row_id=None,
                expression=sp.sstr(expr),
                depth=None,
                validation_reason="leading-order coefficient solve candidate",
                source="coefficient_solve_grammar",
            )
        )

    if solution_set == sp.EmptySet:
        status = "no_leading_order_solution"
    elif skipped_parametric:
        status = "parametric_solution_not_exported"
    else:
        status = "candidate_generated" if rows else "no_candidate_generated"

    return rows, {
        "enabled": True,
        "status": status,
        "ansatz": "Psi = 1 - x + a**2 * sum_i c_i*basis_i(r,x)",
        "basis_source": "targeted_search_grammar.basis",
        "mass_normalization": "M = 1",
        "series": f"coefficient of a**2 in full residual series through O(a**{COEFFICIENT_SOLVE_SERIES_ORDER})",
        "polynomial_variables": ["r", "x"],
        "unknown_count": len(coeff_symbols),
        "equation_count": len(equations),
        "matrix_shape": [int(matrix.rows), int(matrix.cols)],
        "linsolve_result": "EmptySet" if solution_set == sp.EmptySet else sp.sstr(solution_set)[:2000],
        "exact_solution_count": len(exact_solutions),
        "skipped_parametric_solution_count": skipped_parametric,
        "generated_candidates": len(rows),
    }


def sympify_locals() -> dict[str, Any]:
    from expression_operations import UNARY_OPS

    r = sp.Symbol("r", real=True, positive=True)
    x = sp.Symbol("x", real=True)
    m = sp.Symbol("M", real=True, positive=True)
    a = sp.Symbol("a", real=True)
    locals_map: dict[str, Any] = {"r": r, "x": x, "M": m, "a": a}
    locals_map.update(UNARY_OPS)
    return locals_map


def exact_zero(expr: sp.Basic) -> bool:
    simplified = sp.factor(sp.cancel(sp.together(expr)))
    return simplified == 0 or sp.simplify(simplified) == 0


def kerr_surrogate_lhs(u: sp.Basic) -> sp.Basic:
    locals_map = sympify_locals()
    r = locals_map["r"]
    x = locals_map["x"]
    m = locals_map["M"]
    a = locals_map["a"]
    u = u.subs(
        [
            (s, {"r": r, "x": x, "M": m, "a": a}[str(s)])
            for s in u.free_symbols
            if str(s) in {"r", "x", "M", "a"}
        ]
    )
    delta = r**2 - 2 * m * r + a**2
    g = 1 - (2 * m * r) / (r**2 + a**2 * x**2)
    return sp.diff(g / (1 - x**2) * sp.diff(u, r), r) + sp.diff(g / delta * sp.diff(u, x), x)


def full_kerr_split_monopole_residual(psi: sp.Basic) -> sp.Basic:
    """Mahlmann et al. Eq. GSLightCylinder in x = cos(theta) coordinates.

    The target is closed by the split-monopole potential functions used as the
    paper's setup anchor:

        omega(Psi) = (1/2) * a / (r_+**2 + a**2)
        I(Psi) = -(1/2) * omega * Psi * (2 - Psi)

    The GSE uses II' = I(Psi) dI/dPsi.  Since omega is fixed for this closed
    gate, omega_,r and omega_,theta vanish.
    """
    locals_map = sympify_locals()
    r = locals_map["r"]
    x = locals_map["x"]
    m = locals_map["M"]
    a = locals_map["a"]
    psi = psi.subs(
        [
            (s, {"r": r, "x": x, "M": m, "a": a}[str(s)])
            for s in psi.free_symbols
            if str(s) in {"r", "x", "M", "a"}
        ]
    )

    s2 = 1 - x**2
    sigma = r**2 + a**2 * x**2
    delta = r**2 - 2 * m * r + a**2
    big_a = (r**2 + a**2) ** 2 - delta * a**2 * s2
    r_plus = m + sp.sqrt(m**2 - a**2)
    omega = a / (2 * (r_plus**2 + a**2))
    ii_prime = sp.Rational(1, 2) * omega**2 * psi * (2 - psi) * (1 - psi)

    psi_r = sp.diff(psi, r)
    psi_x = sp.diff(psi, x)
    psi_rr = sp.diff(psi, r, 2)
    psi_xx = sp.diff(psi, x, 2)

    sigma_r = sp.diff(sigma, r)
    big_a_r = sp.diff(big_a, r)

    second_order_block = (
        psi_rr
        + (s2 / delta) * psi_xx
        + (big_a_r / big_a - sigma_r / sigma) * psi_r
    )
    light_surface_factor = (
        omega**2 * big_a * s2 / sigma
        - 4 * m * a * r * omega * s2 / sigma
        - 1
        + 2 * m * r / sigma
    )

    term_r_metric = (big_a_r / big_a - sigma_r / sigma) * psi_r
    term_a_theta = 8 * m * a**3 * r * omega * x * s2**2 / (sigma * big_a) * psi_x
    term_sigma_theta = -4 * m * r * a**2 * x * s2 / (delta * sigma**2) * psi_x
    theta_factor_times_psi_theta = (
        -2 * x
        + 2 * delta * a**2 * x * s2 / big_a
        - 2 * a**2 * x * s2 / sigma
    ) * psi_x
    term_light_theta = (
        theta_factor_times_psi_theta
        * big_a
        * omega
        * (omega - 4 * m * a * r / big_a)
        * s2
        / (delta * sigma)
    )
    term_r_drag = -(
        2 * m * r / sigma - 4 * m * a * r * omega * s2 / sigma
    ) * (big_a_r / big_a - 1 / r) * psi_r

    rhs = (
        second_order_block * light_surface_factor
        + term_r_metric
        + term_a_theta
        + term_sigma_theta
        + term_light_theta
        + term_r_drag
    )
    lhs = 4 * sigma / delta * ii_prime
    return rhs - lhs


def full_target_validation(expr: sp.Basic) -> dict[str, Any]:
    try:
        residual = full_kerr_split_monopole_residual(expr)
    except Exception as exc:
        return {
            "residual_defined": False,
            "exact_zero": False,
            "point_checks": [],
            "residual_simplified": f"<residual-error: {exc}>",
            "reason": f"full target residual error: {exc}",
        }

    point_checks = []
    all_points_zero = True
    for point in SAFE_POINTS:
        try:
            value = sp.factor(sp.cancel(sp.together(residual.subs(_point_subs(point)))))
            is_zero = value == 0 or sp.simplify(value) == 0
            value_repr = sp.sstr(value)
        except Exception as exc:
            is_zero = False
            value_repr = f"<point-error: {exc}>"
        if not is_zero:
            all_points_zero = False
        point_checks.append(
            {
                "point": {key: sp.sstr(value) for key, value in point.items()},
                "value": value_repr[:2000],
                "zero": bool(is_zero),
            }
        )

    exact = False
    residual_repr = "<skipped; point checks nonzero>"
    if all_points_zero:
        try:
            simplified = sp.factor(sp.cancel(sp.together(residual)))
            exact = simplified == 0 or sp.simplify(simplified) == 0
            residual_repr = sp.sstr(simplified)[:2000]
        except Exception as exc:
            residual_repr = f"<exact-simplify-error: {exc}>"

    return {
        "residual_defined": True,
        "exact_zero": bool(exact),
        "point_checks": point_checks,
        "residual_simplified": residual_repr,
        "reason": "valid" if exact else FULL_TARGET_REJECTION,
    }


def _point_subs(point: dict[str, sp.Basic]) -> dict[sp.Symbol, sp.Basic]:
    locals_map = sympify_locals()
    return {
        locals_map["M"]: point["M"],
        locals_map["a"]: point["a"],
        locals_map["r"]: point["r"],
        locals_map["x"]: point["x"],
    }


def _has_bad_atom(expr: sp.Basic) -> bool:
    try:
        return bool(expr.has(sp.zoo, sp.oo, -sp.oo, sp.nan))
    except Exception:
        return True


def denominator_rejections(expr: sp.Basic) -> list[str]:
    rejections: list[str] = []
    if _has_bad_atom(expr):
        rejections.append("non-finite symbolic expression")
        return rejections

    try:
        denominator = sp.together(expr).as_numer_denom()[1]
    except Exception as exc:
        return [f"could not extract denominator: {exc}"]

    if _has_bad_atom(denominator):
        rejections.append("non-finite denominator")

    for point in SAFE_POINTS:
        try:
            value = sp.factor(sp.cancel(sp.together(denominator.subs(_point_subs(point)))))
        except Exception as exc:
            rejections.append(f"denominator check failed at {point}: {exc}")
            continue
        if value == 0 or _has_bad_atom(value):
            rejections.append(f"singular denominator at {point}")
    return rejections


def finite_point_rejections(expr: sp.Basic) -> list[str]:
    rejections: list[str] = []
    if _has_bad_atom(expr):
        return ["non-finite symbolic expression"]
    for point in SAFE_POINTS:
        try:
            value = sp.simplify(expr.subs(_point_subs(point)))
        except Exception as exc:
            rejections.append(f"finite-value check failed at {point}: {exc}")
            continue
        if _has_bad_atom(value):
            rejections.append(f"non-finite value at {point}")
            continue
        try:
            numeric = sp.N(value, 40)
            if numeric.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
                rejections.append(f"non-finite numeric value at {point}")
        except Exception as exc:
            rejections.append(f"numeric finite-value check failed at {point}: {exc}")
    return rejections


def small_spin_limit(expr: sp.Basic) -> tuple[str | None, bool]:
    locals_map = sympify_locals()
    a = locals_map["a"]
    x = locals_map["x"]
    targets = [1 - x, x]
    try:
        limit_expr = sp.simplify(sp.limit(expr, a, 0))
    except Exception:
        try:
            limit_expr = sp.simplify(expr.subs(a, 0))
        except Exception:
            return None, False

    matches = False
    for target in targets:
        try:
            if sp.simplify(limit_expr - target) == 0:
                matches = True
                break
        except Exception:
            continue
    return sp.sstr(limit_expr), matches


def equivalent_to_known_anchor(expr: sp.Basic) -> bool:
    locals_map = sympify_locals()
    x = locals_map["x"]
    for target in (1 - x, x):
        try:
            if sp.simplify(expr - target) == 0:
                return True
        except Exception:
            continue
    return False


def surrogate_validation(expr: sp.Basic) -> dict[str, Any]:
    try:
        lhs = kerr_surrogate_lhs(expr)
        point_checks = []
        all_points_zero = True
        for point in SAFE_POINTS:
            value = sp.factor(sp.cancel(sp.together(lhs.subs(_point_subs(point)))))
            is_zero = value == 0 or sp.simplify(value) == 0
            if not is_zero:
                all_points_zero = False
            point_checks.append(
                {
                    "point": {key: sp.sstr(value) for key, value in point.items()},
                    "value": sp.sstr(value)[:500],
                    "zero": bool(is_zero),
                }
            )
    except Exception as exc:
        point_checks = []
        all_points_zero = False
        reason = f"linear surrogate diagnostic error: {exc}"
    else:
        reason = "linear surrogate point checks zero" if all_points_zero else "linear surrogate point checks nonzero"

    return {
        "independent_linear_surrogate_point_zero": bool(all_points_zero),
        "linear_surrogate_point_checks": point_checks,
        "linear_surrogate_reason": reason,
    }


def assess_expression(expression: str, row: CandidateRow | None = None) -> dict[str, Any]:
    locals_map = sympify_locals()
    rejections: list[str] = []
    expr: sp.Basic | None = None
    parse_error = None

    try:
        expr = sp.sympify(expression, locals=locals_map)
    except Exception as exc:
        parse_error = str(exc)
        rejections.append(f"parse error: {exc}")

    variables: list[str] = []
    limit_repr: str | None = None
    limit_matches_anchor = False
    anchor_equivalent = False
    surrogate = {
        "independent_linear_surrogate_point_zero": False,
        "linear_surrogate_point_checks": [],
        "linear_surrogate_reason": "not evaluated",
    }
    full_target = {
        "residual_defined": False,
        "exact_zero": False,
        "point_checks": [],
        "residual_simplified": "<not-evaluated>",
        "reason": "not evaluated",
    }

    if expr is not None:
        variables = sorted(str(sym) for sym in expr.free_symbols)
        missing = sorted(REQUIRED_VARIABLES - set(variables))
        if missing:
            rejections.append("missing required variables: " + ", ".join(missing))

        rejections.extend(denominator_rejections(expr))
        rejections.extend(finite_point_rejections(expr))

        limit_repr, limit_matches_anchor = small_spin_limit(expr)
        if not limit_matches_anchor:
            rejections.append("fails small-spin anchor limit to 1 - x or x")

        anchor_equivalent = equivalent_to_known_anchor(expr)
        if anchor_equivalent:
            rejections.append("equivalent to known small-spin anchor")

        surrogate = surrogate_validation(expr)

        pre_target_rejections = list(rejections)
        if not pre_target_rejections:
            full_target = full_target_validation(expr)
            if not full_target["exact_zero"]:
                rejections.append(full_target["reason"])
        else:
            full_target["reason"] = "skipped because strict prechecks failed"

    deduped_rejections = list(dict.fromkeys(rejections))
    strict_prechecks_pass = not [
        reason for reason in deduped_rejections if reason != FULL_TARGET_REJECTION
    ]
    paper_admissible = (
        PAPER_TARGET_VALIDATOR_IMPLEMENTED
        and strict_prechecks_pass
        and bool(full_target["exact_zero"])
    )

    return {
        "row_id": row.row_id if row else None,
        "source": row.source if row else "manual",
        "expression": expression,
        "parse_error": parse_error,
        "variables": variables,
        "required_variables": sorted(REQUIRED_VARIABLES),
        "small_spin_limit": limit_repr,
        "small_spin_anchor_matches": limit_matches_anchor,
        "equivalent_to_known_anchor": anchor_equivalent,
        "strict_prechecks": {
            "passed_before_full_target_residual": strict_prechecks_pass,
            "rejections_before_full_target_residual": [
                reason for reason in deduped_rejections if reason != FULL_TARGET_REJECTION
            ],
        },
        "linear_surrogate": surrogate,
        "full_target": full_target,
        "paper_admissible": paper_admissible,
        "paper_rejection_reasons": deduped_rejections,
        "engine_depth": row.depth if row else None,
        "engine_validation_reason": row.validation_reason if row else None,
    }


def run_engine(max_depth: int, validators: int, timeout_s: float, validation_timeout_s: float) -> DbRun:
    command = [
        sys.executable,
        "general_method_paper_reproduction.py",
        "--problem",
        "kerr_magnetosphere",
        "--max-depth",
        str(max_depth),
        "--validators",
        str(validators),
    ]
    display_command = ["python3", *command[1:]]
    env = os.environ.copy()
    env["PDE_ENGINE_VALIDATION_TIMEOUT_S"] = str(validation_timeout_s)
    proc = subprocess.run(
        command,
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        timeout=timeout_s,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "engine run failed with return code "
            f"{proc.returncode}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )

    run_ids = re.findall(r"(paper_repro_[0-9_]+_[0-9a-f]+)", proc.stdout)
    if not run_ids:
        raise RuntimeError(f"could not parse run id from engine output:\n{proc.stdout}")
    run_id = run_ids[-1]
    db_matches = re.findall(r"Database:\s*(problems/kerr_magnetosphere/outputs/parallel_runs_[^\s]+\.db)", proc.stdout)
    db_path = ROOT / db_matches[-1] if db_matches else ROOT / "problems" / "kerr_magnetosphere" / "outputs" / f"parallel_runs_{run_id}.db"
    return DbRun(
        db_path=db_path,
        table_name=f"expressions_{run_id.replace('-', '_')}",
        run_id=run_id,
        command=display_command,
        max_depth=max_depth,
        validators=validators,
        wall_timeout_s=timeout_s,
        validation_timeout_s=validation_timeout_s,
    )


def load_rows(run: DbRun, max_rows: int) -> tuple[list[CandidateRow], dict[str, Any]]:
    if not run.db_path.exists():
        raise FileNotFoundError(run.db_path)
    with sqlite3.connect(run.db_path) as conn:
        cur = conn.cursor()
        cur.execute(f"SELECT COUNT(*) FROM {run.table_name}")
        total_generated = int(cur.fetchone()[0])
        cur.execute(f"SELECT COUNT(*) FROM {run.table_name} WHERE validation_status = 'completed'")
        total_completed = int(cur.fetchone()[0])
        cur.execute(f"SELECT COUNT(*) FROM {run.table_name} WHERE is_valid = 1")
        total_valid = int(cur.fetchone()[0])
        cur.execute(f"SELECT COUNT(*) FROM {run.table_name} WHERE is_paper_solution = 1")
        known_rows = int(cur.fetchone()[0])
        cur.execute(
            f"""
            SELECT id, expression, depth, validation_reason
            FROM {run.table_name}
            ORDER BY id
            LIMIT ?
            """,
            (max_rows,),
        )
        rows = [
            CandidateRow(
                row_id=int(row[0]),
                expression=str(row[1]),
                depth=int(row[2]) if row[2] is not None else None,
                validation_reason=row[3],
                source="engine",
            )
            for row in cur.fetchall()
        ]

    return rows, {
        "total_generated": total_generated,
        "total_completed": total_completed,
        "total_valid_rows": total_valid,
        "known_solution_rows": known_rows,
        "rows_loaded_for_full_target_gate": len(rows),
    }


def build_artifact(
    run: DbRun | None,
    rows: list[CandidateRow],
    run_summary: dict[str, Any],
    include_probes: bool,
    include_corrections: bool,
    include_pair_corrections: bool,
    include_coefficient_solve: bool,
) -> dict[str, Any]:
    candidate_rows = list(rows)
    correction_rows: list[CandidateRow] = []
    if include_corrections:
        correction_rows = generate_anchor_correction_rows()
        candidate_rows.extend(correction_rows)
    pair_correction_rows: list[CandidateRow] = []
    if include_pair_corrections:
        pair_correction_rows = generate_anchor_pair_correction_rows()
        candidate_rows.extend(pair_correction_rows)
    coefficient_solve_rows: list[CandidateRow] = []
    coefficient_solve_metadata: dict[str, Any] = {
        "enabled": False,
        "generated_candidates": 0,
    }
    if include_coefficient_solve:
        coefficient_solve_rows, coefficient_solve_metadata = solve_leading_order_coefficient_rows()
        candidate_rows.extend(coefficient_solve_rows)
    if include_probes:
        candidate_rows.extend(
            CandidateRow(None, expr, None, "manual strict-gate probe", "probe")
            for expr in DEFAULT_PROBES
        )

    assessments = [assess_expression(row.expression, row) for row in candidate_rows]
    admitted = [item for item in assessments if item["paper_admissible"]]
    status = "candidate_found" if admitted else "no_candidate_yet"

    return {
        "schema": "kerr_paper_target_gate_v1",
        "status": status,
        "paper_target": PAPER_TARGET,
        "literature_sources": LITERATURE_SOURCES,
        "gate": {
            "required_variables": sorted(REQUIRED_VARIABLES),
            "safe_points": [{key: sp.sstr(value) for key, value in point.items()} for point in SAFE_POINTS],
            "requires_small_spin_anchor": "limit a -> 0 must equal 1 - x or x",
            "rejects_known_anchors": ["1 - x", "x"],
            "rejects_singular_denominators": True,
            "full_target_residual": {
                "source_equation": "Mahlmann et al. 2018 Eq. GSLightCylinder",
                "coordinate_change": "x = cos(theta)",
                "closed_potential_functions": (
                    "split-monopole setup: omega = a/(2*(r_+**2 + a**2)), "
                    "I(Psi)=-(omega/2)*Psi*(2-Psi), II'=I*dI/dPsi"
                ),
                "implemented": PAPER_TARGET_VALIDATOR_IMPLEMENTED,
            },
            "targeted_search_grammar": {
                "enabled": include_corrections,
                "form": "Psi = 1 - x + a**2 * c * basis(r,x)",
                "coefficients": [sp.sstr(item) for item in CORRECTION_COEFFICIENTS],
                "angular_factors": CORRECTION_ANGULAR_FACTORS,
                "radial_factors": CORRECTION_RADIAL_FACTORS,
                "basis": [name for name, _ in correction_basis()],
                "basis_construction": "cartesian product of angular_factors and radial_factors",
                "generated_candidates": len(correction_rows),
            },
            "two_term_search_grammar": {
                "enabled": include_pair_corrections,
                "form": "Psi = 1 - x + a**2 * (c1*basis_i(r,x) + c2*basis_j(r,x))",
                "coefficients": [sp.sstr(item) for item in PAIR_CORRECTION_COEFFICIENTS],
                "angular_factors": PAIR_CORRECTION_ANGULAR_FACTORS,
                "radial_factors": PAIR_CORRECTION_RADIAL_FACTORS,
                "basis": [name for name, _ in pair_correction_basis()],
                "basis_construction": "unordered pairs from the cartesian product basis",
                "generated_candidates": len(pair_correction_rows),
            },
            "coefficient_solve_screen": coefficient_solve_metadata,
        },
        "source_engine_run": None
        if run is None
        else {
            "run_id": run.run_id,
            "db_path": str(run.db_path.relative_to(ROOT)),
            "table_name": run.table_name,
            "command": run.command,
            "search_bounds": {
                "max_depth": run.max_depth,
                "validators": run.validators,
                "wall_timeout_s": run.wall_timeout_s,
                "per_expression_validation_timeout_s": run.validation_timeout_s,
            },
        },
        "run_summary": run_summary,
        "candidate_inputs_scanned": len(assessments),
        "admitted_count": len(admitted),
        "assessments": assessments,
    }


def markdown_from_artifact(artifact: dict[str, Any]) -> str:
    run = artifact.get("source_engine_run")
    if run:
        source_run = f"""- Run id: `{run["run_id"]}`
- Database: `{run["db_path"]}`
- Table: `{run["table_name"]}`
- Command: `{" ".join(run.get("command") or [])}`
- Bounds: `max_depth={run["search_bounds"]["max_depth"]}`, `validators={run["search_bounds"]["validators"]}`, `per_expression_validation_timeout_s={run["search_bounds"]["per_expression_validation_timeout_s"]}`
"""
    else:
        source_run = "- No engine run attached; artifact contains strict-gate probes only.\n"

    grammar = artifact["gate"]["targeted_search_grammar"]
    pair_grammar = artifact["gate"]["two_term_search_grammar"]
    coefficient_solve = artifact["gate"]["coefficient_solve_screen"]
    correction_assessments = [
        item for item in artifact["assessments"] if item["source"] == "anchor_correction_grammar"
    ]
    pair_correction_assessments = [
        item
        for item in artifact["assessments"]
        if item["source"] == "anchor_pair_correction_grammar"
    ]
    correction_prechecks = [
        item
        for item in correction_assessments
        if item["strict_prechecks"]["passed_before_full_target_residual"]
    ]
    pair_correction_prechecks = [
        item
        for item in pair_correction_assessments
        if item["strict_prechecks"]["passed_before_full_target_residual"]
    ]
    correction_exact = [item for item in correction_assessments if item["full_target"]["exact_zero"]]
    pair_correction_exact = [
        item for item in pair_correction_assessments if item["full_target"]["exact_zero"]
    ]
    coeffs = ", ".join(grammar["coefficients"])
    pair_coeffs = ", ".join(pair_grammar["coefficients"])
    angular_lines = "\n".join(f"- `{item}`" for item in grammar["angular_factors"])
    radial_lines = "\n".join(f"- `{item}`" for item in grammar["radial_factors"])
    pair_angular_lines = "\n".join(f"- `{item}`" for item in pair_grammar["angular_factors"])
    pair_radial_lines = "\n".join(f"- `{item}`" for item in pair_grammar["radial_factors"])
    if coefficient_solve.get("enabled"):
        coefficient_section = f"""
The gate then runs a leading-order coefficient solve over the full one-term
basis instead of only trying fixed scalar coefficients:

```text
{coefficient_solve["ansatz"]}
{coefficient_solve["series"]}
{coefficient_solve["mass_normalization"]}
```

This screen produced `{coefficient_solve["equation_count"]}` polynomial
equations in `{coefficient_solve["unknown_count"]}` unknown coefficients.  The
linear system matrix shape was `{coefficient_solve["matrix_shape"]}`, and SymPy
returned `{coefficient_solve["linsolve_result"]}`.  It therefore generated
`{coefficient_solve["generated_candidates"]}` additional candidates.
"""
    else:
        coefficient_section = "\nThe leading-order coefficient solve screen was disabled for this artifact.\n"
    grammar_section = f"""## Targeted finite-spin correction grammar

The gate also appends an expanded bounded correction grammar around the
split-monopole anchor:

```text
{grammar["form"]}
c in {{{coeffs}}}
basis = angular_factor * radial_factor
```

Angular factors:

{angular_lines}

Radial factors:

{radial_lines}

This generated `{grammar["generated_candidates"]}` correction candidates. In
this run, `{len(correction_prechecks)}` passed the strict prechecks before
the full residual, and `{len(correction_exact)}` had exact-zero full residual.
The full Cartesian-product basis list is stored in the JSON artifact.

The gate then appends a bounded two-term correction screen:

```text
{pair_grammar["form"]}
c1, c2 in {{{pair_coeffs}}}
```

Two-term angular factors:

{pair_angular_lines}

Two-term radial factors:

{pair_radial_lines}

This generated `{pair_grammar["generated_candidates"]}` two-term correction
candidates. In this run, `{len(pair_correction_prechecks)}` passed strict
prechecks before the full residual, and `{len(pair_correction_exact)}` had
exact-zero full residual.
{coefficient_section}
"""

    rows = []
    for item in artifact["assessments"][:20]:
        reasons = "; ".join(item["paper_rejection_reasons"][:3])
        rows.append(
            "| `{}` | {} | {} | {} | {} |".format(
                item["expression"],
                ", ".join(item["variables"]) or "-",
                item["strict_prechecks"]["passed_before_full_target_residual"],
                item["full_target"]["exact_zero"],
                reasons,
            )
        )
    table = "\n".join(rows)

    return f"""# Kerr paper-target gate

Status: **{artifact["status"]}**.

This artifact is the strict admission gate for the next paper target. The
current `kerr_magnetosphere` linear surrogate is used only to generate an
expression table; admission is decided by the nonlinear Kerr force-free
Grad-Shafranov residual gate below.

## Literature target

- Mahlmann et al. use the relativistic Grad-Shafranov equation as the central
  target for static, axisymmetric, force-free Kerr magnetospheres and emphasize
  numerical solution methodology for standard setups.
  Source: https://arxiv.org/abs/1802.00815
- Camilloni et al. state that for extreme Kerr there is no known exact analytic
  force-free solution that is stationary, axisymmetric, and magnetically
  dominated.
  Source: https://arxiv.org/abs/2007.15665

## Gate boundary

A row is not a paper candidate unless it passes the closed split-monopole form
of the nonlinear Kerr force-free Grad-Shafranov equation:

- depends on `r`, `x`, and `a`
- has no singular denominator on the rational safe points
- is finite on the rational safe points
- has small-spin limit `1 - x` or `x`
- is not exactly the known small-spin anchor
- passes the full nonlinear Kerr force-free Grad-Shafranov residual from
  Mahlmann et al. Eq. `GSLightCylinder`, rewritten with `x = cos(theta)`
- uses the split-monopole potential functions
  `omega = a/(2*(r_+**2 + a**2))` and
  `I(Psi)=-(omega/2)*Psi*(2-Psi)`

The finite-spin paper target remains intentionally hard: the Schwarzschild
anchor is used only as an `a -> 0` limit, not as a finite-spin solution.

{grammar_section}

## Source run

{source_run}
## Candidate assessments

Only the first 20 assessments are shown here; the JSON contains the full list.

| Expression | Variables | Strict prechecks before full residual | Full residual exact zero | First rejection reasons |
| --- | --- | --- | --- | --- |
{table}
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-engine", action="store_true", help="Run a fresh bounded kerr_magnetosphere search first.")
    parser.add_argument("--max-depth", type=int, default=2)
    parser.add_argument("--validators", type=int, default=1)
    parser.add_argument("--timeout-s", type=float, default=90)
    parser.add_argument("--validation-timeout-s", type=float, default=3)
    parser.add_argument("--db-path", type=Path)
    parser.add_argument("--table-name")
    parser.add_argument("--run-id", default="manual")
    parser.add_argument("--max-rows", type=int, default=1000)
    parser.add_argument("--include-probes", action="store_true", default=True)
    parser.add_argument("--no-corrections", action="store_true", help="Do not add the bounded anchor-preserving correction grammar.")
    parser.add_argument(
        "--no-pair-corrections",
        action="store_true",
        help="Do not add the bounded two-term anchor-preserving correction grammar.",
    )
    parser.add_argument(
        "--no-coefficient-solve",
        action="store_true",
        help="Do not run the leading-order coefficient solve screen.",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--markdown", type=Path, default=DEFAULT_MD)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run: DbRun | None = None
    rows: list[CandidateRow] = []
    summary: dict[str, Any] = {
        "total_generated": 0,
        "total_completed": 0,
        "total_valid_rows": 0,
        "known_solution_rows": 0,
        "rows_loaded_for_full_target_gate": 0,
    }

    if args.run_engine:
        run = run_engine(args.max_depth, args.validators, args.timeout_s, args.validation_timeout_s)
        rows, summary = load_rows(run, args.max_rows)
    elif args.db_path or args.table_name:
        if not args.db_path or not args.table_name:
            raise SystemExit("--db-path and --table-name must be provided together")
        db_path = args.db_path if args.db_path.is_absolute() else ROOT / args.db_path
        run = DbRun(
            db_path=db_path,
            table_name=args.table_name,
            run_id=args.run_id,
            command=[
                "python3",
                "general_method_paper_reproduction.py",
                "--problem",
                "kerr_magnetosphere",
                "--max-depth",
                str(args.max_depth),
                "--validators",
                str(args.validators),
            ],
            max_depth=args.max_depth,
            validators=args.validators,
            wall_timeout_s=args.timeout_s,
            validation_timeout_s=args.validation_timeout_s,
        )
        rows, summary = load_rows(run, args.max_rows)

    artifact = build_artifact(
        run,
        rows,
        summary,
        include_probes=args.include_probes,
        include_corrections=not args.no_corrections,
        include_pair_corrections=not args.no_pair_corrections,
        include_coefficient_solve=not args.no_coefficient_solve,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.markdown.write_text(markdown_from_artifact(artifact), encoding="utf-8")
    print(args.output)
    print(args.markdown)
    print(f"status={artifact['status']} admitted={artifact['admitted_count']} scanned={artifact['candidate_inputs_scanned']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
