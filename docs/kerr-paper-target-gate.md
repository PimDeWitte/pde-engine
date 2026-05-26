# Kerr paper-target gate

Status: **no_candidate_yet**.

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

Current full-target blocker: `full nonlinear Kerr force-free Grad-Shafranov validator is not implemented`.

## Source run

- Run id: `paper_repro_20260526_091857_d29966e9`
- Database: `problems/kerr_magnetosphere/outputs/parallel_runs_paper_repro_20260526_091857_d29966e9.db`
- Table: `expressions_paper_repro_20260526_091857_d29966e9`
- Command: `python3 general_method_paper_reproduction.py --problem kerr_magnetosphere --max-depth 2 --validators 1`
- Bounds: `max_depth=2`, `validators=1`, `per_expression_validation_timeout_s=3.0`

## Candidate assessments

Only the first 20 assessments are shown here; the JSON contains the full list.

| Expression | Variables | Strict prechecks before target blocker | Linear surrogate valid | First rejection reasons |
| --- | --- | --- | --- | --- |
| `1 - x` | x | False | False | missing required variables: a, r; equivalent to known small-spin anchor; does not pass current repo linear surrogate heavy validation |
| `x` | x | False | False | missing required variables: a, r; equivalent to known small-spin anchor; does not pass current repo linear surrogate heavy validation |
| `1/(1 - 1)` | - | False | False | missing required variables: a, r, x; non-finite symbolic expression; fails small-spin anchor limit to 1 - x or x |
| `1 - x + a**2*r*x` | a, r, x | False | False | does not pass current repo linear surrogate heavy validation; full nonlinear Kerr force-free Grad-Shafranov validator is not implemented |
