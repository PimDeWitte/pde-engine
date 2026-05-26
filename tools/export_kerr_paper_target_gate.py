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
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
DEFAULT_JSON = DOCS / "kerr-paper-target-gate.json"
DEFAULT_MD = DOCS / "kerr-paper-target-gate.md"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


PAPER_TARGET_VALIDATOR_IMPLEMENTED = False
TARGET_BLOCKER = "full nonlinear Kerr force-free Grad-Shafranov validator is not implemented"
REQUIRED_VARIABLES = {"r", "x", "a"}

DEFAULT_PROBES = [
    "1 - x",
    "x",
    "1/(1 - 1)",
    "1 - x + a**2*r*x",
]

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
    "status": "not_implemented_in_repo",
    "validator_implemented": PAPER_TARGET_VALIDATOR_IMPLEMENTED,
    "claim": (
        "No pde-engine row is a paper candidate until the full nonlinear Kerr "
        "force-free target equation, physical domain, and regularity gates are encoded."
    ),
    "current_repo_problem": "kerr_magnetosphere",
    "current_repo_problem_boundary": (
        "linear surrogate only; useful as a harness, not sufficient for the "
        "paper target"
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
    from problems.kerr_magnetosphere.validator import KerrMagnetosphereValidator

    try:
        lhs = kerr_surrogate_lhs(expr)
        lhs_simplified = sp.factor(sp.cancel(sp.together(lhs)))
        independent_exact_zero = exact_zero(lhs)
    except Exception as exc:
        lhs_simplified = f"<lhs-error: {exc}>"
        independent_exact_zero = False

    try:
        locals_map = sympify_locals()
        validator = KerrMagnetosphereValidator(
            locals_map["r"],
            locals_map["x"],
            locals_map["M"],
            locals_map["a"],
            M_value=sp.Integer(1),
            a_value=sp.Rational(1, 10),
            use_lean=False,
        )
        is_valid, reason = validator.validate(
            expr,
            check_regularity=True,
            fast_point_only=False,
            lean_first=False,
            defer_heavy_checks=False,
            enforce_anchor=True,
        )
    except Exception as exc:
        is_valid = False
        reason = f"surrogate validator error: {exc}"

    return {
        "independent_linear_surrogate_exact_zero": bool(independent_exact_zero),
        "linear_surrogate_lhs_simplified": sp.sstr(lhs_simplified)[:2000],
        "repo_linear_surrogate_valid": bool(is_valid),
        "repo_linear_surrogate_reason": reason,
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
        "independent_linear_surrogate_exact_zero": False,
        "linear_surrogate_lhs_simplified": "<not-evaluated>",
        "repo_linear_surrogate_valid": False,
        "repo_linear_surrogate_reason": "not evaluated",
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
        if not surrogate["repo_linear_surrogate_valid"]:
            rejections.append("does not pass current repo linear surrogate heavy validation")

    if not PAPER_TARGET_VALIDATOR_IMPLEMENTED:
        rejections.append(TARGET_BLOCKER)

    deduped_rejections = list(dict.fromkeys(rejections))
    strict_prechecks_pass = not [
        reason
        for reason in deduped_rejections
        if reason != TARGET_BLOCKER
    ]

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
            "passed_before_full_target_blocker": strict_prechecks_pass,
            "rejections_before_full_target_blocker": [
                reason for reason in deduped_rejections if reason != TARGET_BLOCKER
            ],
        },
        "linear_surrogate": surrogate,
        "paper_admissible": False if not PAPER_TARGET_VALIDATOR_IMPLEMENTED else strict_prechecks_pass,
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
            WHERE is_valid = 1
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
        "valid_rows_loaded_for_strict_gate": len(rows),
    }


def build_artifact(
    run: DbRun | None,
    rows: list[CandidateRow],
    run_summary: dict[str, Any],
    include_probes: bool,
) -> dict[str, Any]:
    candidate_rows = list(rows)
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
            "requires_full_target_validator": True,
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

    rows = []
    for item in artifact["assessments"][:20]:
        reasons = "; ".join(item["paper_rejection_reasons"][:3])
        rows.append(
            "| `{}` | {} | {} | {} | {} |".format(
                item["expression"],
                ", ".join(item["variables"]) or "-",
                item["strict_prechecks"]["passed_before_full_target_blocker"],
                item["linear_surrogate"]["repo_linear_surrogate_valid"],
                reasons,
            )
        )
    table = "\n".join(rows)

    return f"""# Kerr paper-target gate

Status: **{artifact["status"]}**.

This artifact is the strict admission gate for the next paper target. It does
not upgrade the current `kerr_magnetosphere` linear surrogate into a nonlinear
Kerr force-free / Grad-Shafranov result.

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

The current repo problem is explicitly a linear surrogate. A row is not a paper
candidate unless the full nonlinear target validator exists and the expression
passes the following checks:

- depends on `r`, `x`, and `a`
- has no singular denominator on the rational safe points
- is finite on the rational safe points
- has small-spin limit `1 - x` or `x`
- is not exactly the known small-spin anchor
- passes axis, horizon, and target-specific regularity checks
- passes the full nonlinear Kerr force-free Grad-Shafranov residual

Current full-target blocker: `{TARGET_BLOCKER}`.

## Source run

{source_run}
## Candidate assessments

Only the first 20 assessments are shown here; the JSON contains the full list.

| Expression | Variables | Strict prechecks before target blocker | Linear surrogate valid | First rejection reasons |
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
    parser.add_argument("--max-rows", type=int, default=200)
    parser.add_argument("--include-probes", action="store_true", default=True)
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
        "valid_rows_loaded_for_strict_gate": 0,
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

    artifact = build_artifact(run, rows, summary, include_probes=args.include_probes)
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
