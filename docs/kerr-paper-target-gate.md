# Kerr paper-target gate

Status: **no_candidate_yet**.

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

## Targeted finite-spin correction grammar

The gate also appends a bounded correction grammar around the split-monopole
anchor:

```text
Psi = 1 - x + a**2 * c * basis(r,x)
c in {-1, -1/2, 1/2, 1}
```

Basis functions:

- `x*(1-x**2)/r`
- `x*(1-x**2)/r**2`
- `x*(1-x**2)/(r - 2*M)`
- `(1-x**2)/r`
- `(1-x**2)/r**2`
- `x/r`
- `x/r**2`

This generated `28` correction candidates. In
this run, `28` passed the strict prechecks before the
full residual, and `0` had exact-zero full residual.


## Source run

- Run id: `paper_repro_20260526_093939_80027ff6`
- Database: `problems/kerr_magnetosphere/outputs/parallel_runs_paper_repro_20260526_093939_80027ff6.db`
- Table: `expressions_paper_repro_20260526_093939_80027ff6`
- Command: `python3 general_method_paper_reproduction.py --problem kerr_magnetosphere --max-depth 2 --validators 1`
- Bounds: `max_depth=2`, `validators=1`, `per_expression_validation_timeout_s=3.0`

## Candidate assessments

Only the first 20 assessments are shown here; the JSON contains the full list.

| Expression | Variables | Strict prechecks before full residual | Full residual exact zero | First rejection reasons |
| --- | --- | --- | --- | --- |
| `r` | r | False | False | missing required variables: a, x; fails small-spin anchor limit to 1 - x or x |
| `x` | x | False | False | missing required variables: a, r; equivalent to known small-spin anchor |
| `1` | - | False | False | missing required variables: a, r, x; fails small-spin anchor limit to 1 - x or x |
| `1/3` | - | False | False | missing required variables: a, r, x; fails small-spin anchor limit to 1 - x or x |
| `1 - x` | x | False | False | missing required variables: a, r; equivalent to known small-spin anchor |
| `a**2` | a | False | False | missing required variables: r, x; fails small-spin anchor limit to 1 - x or x |
| `a**2*x**2 + r**2` | a, r, x | False | False | fails small-spin anchor limit to 1 - x or x |
| `-2*M*r + a**2 + r**2` | M, a, r | False | False | missing required variables: x; fails small-spin anchor limit to 1 - x or x |
| `-2*M*r/(a**2*x**2 + r**2) + 1` | M, a, r, x | False | False | fails small-spin anchor limit to 1 - x or x |
| `neg(r)` | r | False | False | missing required variables: a, x; fails small-spin anchor limit to 1 - x or x |
| `inv(r)` | r | False | False | missing required variables: a, x; fails small-spin anchor limit to 1 - x or x |
| `sqrt(r)` | r | False | False | missing required variables: a, x; fails small-spin anchor limit to 1 - x or x |
| `square(r)` | r | False | False | missing required variables: a, x; fails small-spin anchor limit to 1 - x or x |
| `pow_3_2(r)` | r | False | False | missing required variables: a, x; fails small-spin anchor limit to 1 - x or x |
| `pow_neg_3_2(r)` | r | False | False | missing required variables: a, x; fails small-spin anchor limit to 1 - x or x |
| `exp(r)` | r | False | False | missing required variables: a, x; fails small-spin anchor limit to 1 - x or x |
| `exp_neg(r)` | r | False | False | missing required variables: a, x; fails small-spin anchor limit to 1 - x or x |
| `neg(x)` | x | False | False | missing required variables: a, r; fails small-spin anchor limit to 1 - x or x |
| `inv(x)` | x | False | False | missing required variables: a, r; fails small-spin anchor limit to 1 - x or x |
| `sqrt(x)` | x | False | False | missing required variables: a, r; fails small-spin anchor limit to 1 - x or x |
