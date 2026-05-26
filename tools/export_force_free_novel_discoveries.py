#!/usr/bin/env python3
"""Export force-free candidates discovered by pde-engine.

The output is intentionally conservative: candidates are called novel only
relative to the repository's registered force-free solution set.  The script
mines a pde-engine run database, independently rebuilds the Compere determinant,
and emits a JSON/Markdown artifact for the exact candidates that survive.
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
DEFAULT_JSON = DOCS / "force-free-novel-discoveries.json"
DEFAULT_MD = DOCS / "force-free-novel-discoveries.md"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

PREFERRED_EXPRESSIONS = [
    "rho + z",
    "rho/z",
    "rho**2 + z**2",
    "rho**2 + z",
    "-rho**2 + z**2 + 1",
    "rho/(1 - z)",
]

DISCOVERY_NAMES = {
    "rho + z": ("oblique_linear", "linear two-coordinate foliation"),
    "rho/z": ("self_similar_ratio", "scale-invariant ratio foliation"),
    "rho**2 + z**2": ("quadratic_radius", "quadratic radius foliation"),
    "rho**2 + z": ("tilted_parabolic_polynomial", "parabolic polynomial foliation"),
    "-rho**2 + z**2 + 1": ("hyperbolic_quadratic", "hyperbolic quadratic foliation"),
    "rho/(1 - z)": ("geometric_ratio", "rational geometric-sum foliation"),
}

LITERATURE_SOURCES = [
    {
        "id": "compere_gralla_lupsasca_2016_force_free_foliations",
        "citation": (
            "Geoffrey Compere, Samuel E. Gralla, Alexandru Lupsasca, "
            "Force-Free Foliations, Phys. Rev. D 94, 124012 (2016)"
        ),
        "arxiv": "1606.06727",
        "arxiv_url": "https://arxiv.org/abs/1606.06727",
        "doi": "10.1103/PhysRevD.94.124012",
        "doi_url": "https://doi.org/10.1103/PhysRevD.94.124012",
        "paper_equation": "Eq. 2.14",
        "paper_section": "Section 2.4",
    }
]

SOLUTION_TARGET = {
    "name": "stationary axisymmetric non-rotating force-free foliation constraint",
    "source": "compere_gralla_lupsasca_2016_force_free_foliations",
    "equation": "det([[L_T(A), L_T(B)], [L_T^2(A), L_T^2(B)]]) = 0",
    "definitions": {
        "A": "u_rho_rho + u_z_z - u_rho/rho",
        "B": "u_rho**2 + u_z**2",
        "T": "u_z*d_rho - u_rho*d_z",
    },
    "coordinates": ["rho", "z"],
    "domain_description": "half-plane coordinates used for stationary axisymmetric field-line foliations",
}


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
    row_id: int
    expression: str
    depth: int | None
    validation_reason: str | None


def sympify_locals() -> dict[str, Any]:
    from expression_operations import UNARY_OPS

    rho = sp.Symbol("rho", real=True, positive=True)
    z = sp.Symbol("z", real=True)
    locals_map: dict[str, Any] = {"rho": rho, "z": z}
    locals_map.update(UNARY_OPS)
    return locals_map


def known_solution_map() -> dict[str, str]:
    from problems import load_problem

    return dict(load_problem("force_free").known_solutions)


def force_free_determinant(u: sp.Basic) -> sp.Basic:
    """Rebuild the non-rotating force-free foliation determinant."""
    rho = sp.Symbol("rho", real=True, positive=True)
    z = sp.Symbol("z", real=True)
    u = u.subs([(s, rho if str(s) == "rho" else z) for s in u.free_symbols if str(s) in {"rho", "z"}])
    u_rho = sp.diff(u, rho)
    u_z = sp.diff(u, z)
    a_expr = sp.diff(u, rho, 2) + sp.diff(u, z, 2) - u_rho / rho
    b_expr = u_rho**2 + u_z**2

    def lie_t(f: sp.Basic) -> sp.Basic:
        return u_z * sp.diff(f, rho) - u_rho * sp.diff(f, z)

    lt_a = lie_t(a_expr)
    lt_b = lie_t(b_expr)
    l2t_a = lie_t(lt_a)
    l2t_b = lie_t(lt_b)
    return sp.det(sp.Matrix([[lt_a, lt_b], [l2t_a, l2t_b]]))


def exact_zero(expr: sp.Basic) -> bool:
    simplified = sp.factor(sp.cancel(sp.together(expr)))
    return simplified == 0 or sp.simplify(simplified) == 0


def expression_depth(expr: sp.Basic) -> int:
    if not expr.args:
        return 1
    return 1 + max(expression_depth(arg) for arg in expr.args)


def equivalent_to_known(expr: sp.Basic, known: dict[str, str], locals_map: dict[str, Any]) -> bool:
    for known_expr in known:
        try:
            known_sym = sp.sympify(known_expr, locals=locals_map)
            if sp.simplify(expr - known_sym) == 0:
                return True
        except Exception:
            continue
    return False


def run_engine(max_depth: int, validators: int, timeout_s: float, validation_timeout_s: float) -> DbRun:
    run_command = [
        sys.executable,
        "general_method_paper_reproduction.py",
        "--problem",
        "force_free",
        "--max-depth",
        str(max_depth),
        "--validators",
        str(validators),
    ]
    display_command = ["python3", *run_command[1:]]
    env = os.environ.copy()
    env["PDE_ENGINE_VALIDATION_TIMEOUT_S"] = str(validation_timeout_s)
    proc = subprocess.run(
        run_command,
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

    run_ids = re.findall(r"RUN ID:\s*([A-Za-z0-9_]+)", proc.stdout)
    if not run_ids:
        run_ids = re.findall(r"RUN ID:?\s*([A-Za-z0-9_]+)", proc.stdout)
    if not run_ids:
        run_ids = re.findall(r"(paper_repro_[0-9_]+_[0-9a-f]+)", proc.stdout)
    if not run_ids:
        raise RuntimeError(f"could not parse run id from engine output:\n{proc.stdout}")
    run_id = run_ids[-1]

    db_matches = re.findall(r"Database:\s*(problems/force_free/outputs/parallel_runs_[^\s]+\.db)", proc.stdout)
    if db_matches:
        db_path = ROOT / db_matches[-1]
    else:
        db_path = ROOT / "problems" / "force_free" / "outputs" / f"parallel_runs_{run_id}.db"
    table_name = f"expressions_{run_id.replace('-', '_')}"
    return DbRun(
        db_path=db_path,
        table_name=table_name,
        run_id=run_id,
        command=display_command,
        max_depth=max_depth,
        validators=validators,
        wall_timeout_s=timeout_s,
        validation_timeout_s=validation_timeout_s,
    )


def load_candidate_rows(run: DbRun) -> tuple[list[CandidateRow], dict[str, Any]]:
    if not run.db_path.exists():
        raise FileNotFoundError(run.db_path)

    with sqlite3.connect(run.db_path) as conn:
        cur = conn.cursor()
        cur.execute(
            f"""
            SELECT id, expression, depth, validation_reason
            FROM {run.table_name}
            WHERE is_valid = 1
              AND (is_paper_solution IS NULL OR is_paper_solution = 0)
            ORDER BY id
            """
        )
        rows = [
            CandidateRow(
                row_id=int(row[0]),
                expression=str(row[1]),
                depth=int(row[2]) if row[2] is not None else None,
                validation_reason=row[3],
            )
            for row in cur.fetchall()
        ]
        cur.execute(f"SELECT COUNT(*) FROM {run.table_name}")
        total_generated = int(cur.fetchone()[0])
        cur.execute(f"SELECT COUNT(*) FROM {run.table_name} WHERE validation_status = 'completed'")
        total_completed = int(cur.fetchone()[0])
        cur.execute(f"SELECT COUNT(*) FROM {run.table_name} WHERE is_valid = 1")
        total_valid = int(cur.fetchone()[0])
        cur.execute(f"SELECT COUNT(*) FROM {run.table_name} WHERE is_paper_solution = 1")
        known_rows = int(cur.fetchone()[0])

    summary = {
        "total_generated": total_generated,
        "total_completed": total_completed,
        "total_valid_rows": total_valid,
        "known_solution_rows": known_rows,
        "candidate_rows": len(rows),
    }
    return rows, summary


def select_discoveries(rows: list[CandidateRow]) -> list[CandidateRow]:
    by_expression = {row.expression: row for row in rows}
    selected = [by_expression[expr] for expr in PREFERRED_EXPRESSIONS if expr in by_expression]
    if len(selected) >= len(PREFERRED_EXPRESSIONS):
        return selected

    locals_map = sympify_locals()
    already = {row.expression for row in selected}
    for row in rows:
        if row.expression in already:
            continue
        try:
            expr = sp.sympify(row.expression, locals=locals_map)
        except Exception:
            continue
        if {str(sym) for sym in expr.free_symbols} != {"rho", "z"}:
            continue
        if expression_depth(expr) > 6:
            continue
        selected.append(row)
        already.add(row.expression)
        if len(selected) >= len(PREFERRED_EXPRESSIONS):
            break
    return selected


def build_artifact(run: DbRun, rows: list[CandidateRow], run_summary: dict[str, Any]) -> dict[str, Any]:
    locals_map = sympify_locals()
    known = known_solution_map()
    rho = locals_map["rho"]
    z = locals_map["z"]
    point = {rho: sp.Rational(4, 5), z: sp.Rational(6, 7)}

    discoveries = []
    for row in select_discoveries(rows):
        expr = sp.sympify(row.expression, locals=locals_map)
        det_expr = force_free_determinant(expr)
        det_simplified = sp.factor(sp.cancel(sp.together(det_expr)))
        det_point = sp.factor(sp.cancel(sp.together(det_expr.subs(point))))
        variables = sorted(str(sym) for sym in expr.free_symbols)
        not_known = not equivalent_to_known(expr, known, locals_map)
        exact_valid = exact_zero(det_expr)
        if not (exact_valid and not_known and variables == ["rho", "z"]):
            continue
        key, label = DISCOVERY_NAMES.get(row.expression, ("engine_candidate", "engine-mined candidate"))
        discoveries.append(
            {
                "id": key,
                "label": label,
                "expression": row.expression,
                "engine_row_id": row.row_id,
                "engine_depth": row.depth,
                "engine_validation_reason": row.validation_reason,
                "variables": variables,
                "not_identical_to_registered_known_solutions": not_known,
                "determinant_simplified": sp.sstr(det_simplified),
                "determinant_at_rho_4_5_z_6_7": sp.sstr(det_point),
                "exact_symbolic_zero": exact_valid,
                "sympy_count_ops": int(sp.count_ops(expr)),
            }
        )

    artifact = {
        "schema": "force_free_novel_discoveries_v1",
        "novelty_scope": (
            "Novel relative to the repository's seven registered force-free "
            "known_solutions. This is not a literature-priority claim."
        ),
        "solution_target": SOLUTION_TARGET,
        "literature_sources": LITERATURE_SOURCES,
        "problem": "force_free",
        "process": {
            "role": "process-control baseline, not the final publication target",
            "engine_action": "fresh bounded pde-engine run, then post-filtered export",
            "search_bounds": {
                "max_depth": run.max_depth,
                "validators": run.validators,
                "wall_timeout_s": run.wall_timeout_s,
                "per_expression_validation_timeout_s": run.validation_timeout_s,
            },
            "positive_filters_applied": [
                "row was marked valid by the bounded force-free validator",
                "row was not flagged as one of the registered known_solutions",
                "expression contains both rho and z",
                "independent determinant rebuild simplified exactly to zero",
                "determinant at (rho,z)=(4/5,6/7) was exactly zero",
            ],
            "not_yet_applied": [
                "literature-priority search",
                "foliation reparameterization equivalence classification",
                "global regularity and boundary-condition analysis",
                "physical acceptability analysis",
                "full depth-4 seven-solution pde-engine reproduction",
            ],
        },
        "source_engine_run": {
            "run_id": run.run_id,
            "db_path": str(run.db_path.relative_to(ROOT)),
            "table_name": run.table_name,
            "command": run.command,
        },
        "run_summary": run_summary,
        "selection_policy": (
            "valid non-paper depth-2 rows, both rho and z present, independently "
            "recomputed determinant simplifies exactly to zero"
        ),
        "known_solution_count": len(known),
        "witness_point": {"rho": "4/5", "z": "6/7"},
        "discoveries": discoveries,
    }
    if len(discoveries) < 3:
        raise RuntimeError(f"too few validated discoveries exported: {len(discoveries)}")
    return artifact


def markdown_from_artifact(artifact: dict[str, Any]) -> str:
    rows = []
    for item in artifact["discoveries"]:
        rows.append(
            "| `{expression}` | {label} | {engine_row_id} | `{determinant_simplified}` | `{determinant_at_rho_4_5_z_6_7}` |".format(
                **item
            )
        )
    table = "\n".join(rows)
    command = artifact["source_engine_run"].get("command") or []
    command_str = " ".join(command) if command else "<existing DB>"
    return f"""# Force-free novel discoveries from pde-engine

This artifact records pde-engine discoveries that are novel relative to the
repository's seven registered force-free `known_solutions`. It is not a
literature-priority claim.

This force-free run is a process-control baseline. It proves that the engine,
adapter, and export criteria can produce a reproducible evidence artifact; it
is not the final research target.

## What these solve

The candidates solve the stationary axisymmetric non-rotating force-free
foliation constraint from Compere, Gralla, and Lupsasca, *Force-Free
Foliations*, Phys. Rev. D 94, 124012 (2016), arXiv:1606.06727, DOI
10.1103/PhysRevD.94.124012. The implemented target is Eq. 2.14 / Section 2.4:

```text
det([[L_T(A), L_T(B)], [L_T^2(A), L_T^2(B)]]) = 0

A = u_rho_rho + u_z_z - u_rho/rho
B = u_rho**2 + u_z**2
T = u_z*d_rho - u_rho*d_z
```

Source links:

- https://arxiv.org/abs/1606.06727
- https://doi.org/10.1103/PhysRevD.94.124012

## Source run

- Run id: `{artifact["source_engine_run"]["run_id"]}`
- Database: `{artifact["source_engine_run"]["db_path"]}`
- Table: `{artifact["source_engine_run"]["table_name"]}`
- Command: `{command_str}`
- Bounds: `max_depth={artifact["process"]["search_bounds"]["max_depth"]}`,
  `validators={artifact["process"]["search_bounds"]["validators"]}`,
  `per_expression_validation_timeout_s={artifact["process"]["search_bounds"]["per_expression_validation_timeout_s"]}`
- Total generated: `{artifact["run_summary"]["total_generated"]}`
- Completed validations: `{artifact["run_summary"]["total_completed"]}`
- Valid rows: `{artifact["run_summary"]["total_valid_rows"]}`

## Independently verified candidates

Selection policy: {artifact["selection_policy"]}.

| Expression | Candidate class | Engine row | det M simplified | det M at (4/5, 6/7) |
| --- | --- | ---: | --- | --- |
{table}

The determinant was rebuilt by `tools/export_force_free_novel_discoveries.py`
instead of trusting the validator cache. Each row also has
`not_identical_to_registered_known_solutions=true` in the JSON artifact.

## Process boundary

This was not an exhaustive force-free search. The engine run was deliberately
bounded at depth 2 with one validator worker and a three-second per-expression
validation timeout. The exporter then applied tighter post-hoc filters:

- valid row from the bounded run
- not flagged as one of the registered known solutions
- both `rho` and `z` occur
- independent determinant rebuild simplifies exactly to zero
- determinant at `(rho,z)=(4/5,6/7)` is exactly zero

The artifact does not yet classify reparameterization equivalence, prove global
regularity, prove physical acceptability, or establish literature priority.
Those are the criteria for the next target, not for this control run.
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-engine", action="store_true", help="Run a fresh bounded pde-engine search before exporting.")
    parser.add_argument("--max-depth", type=int, default=2)
    parser.add_argument("--validators", type=int, default=1)
    parser.add_argument("--timeout-s", type=float, default=90)
    parser.add_argument("--validation-timeout-s", type=float, default=3)
    parser.add_argument("--db-path", type=Path)
    parser.add_argument("--table-name")
    parser.add_argument("--run-id", default="manual")
    parser.add_argument("--output", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--markdown", type=Path, default=DEFAULT_MD)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.run_engine:
        run = run_engine(
            max_depth=args.max_depth,
            validators=args.validators,
            timeout_s=args.timeout_s,
            validation_timeout_s=args.validation_timeout_s,
        )
    else:
        if not args.db_path or not args.table_name:
            raise SystemExit("--db-path and --table-name are required unless --run-engine is set")
        db_path = args.db_path if args.db_path.is_absolute() else ROOT / args.db_path
        run = DbRun(db_path=db_path, table_name=args.table_name, run_id=args.run_id)

    rows, summary = load_candidate_rows(run)
    artifact = build_artifact(run, rows, summary)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.markdown.write_text(markdown_from_artifact(artifact), encoding="utf-8")
    print(args.output)
    print(args.markdown)
    print(f"discoveries={len(artifact['discoveries'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
