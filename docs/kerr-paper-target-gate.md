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

The gate also appends an expanded bounded correction grammar around the
split-monopole anchor:

```text
Psi = 1 - x + a**2 * c * basis(r,x)
c in {-1, -1/2, 1/2, 1}
basis = angular_factor * radial_factor
```

Angular factors:

- `1`
- `x`
- `x**2`
- `1-x**2`
- `x*(1-x**2)`
- `(1-x**2)**2`
- `x**2*(1-x**2)`
- `x*(1-x**2)**2`

Radial factors:

- `1/r`
- `1/r**2`
- `1/r**3`
- `1/(r - 2*M)`
- `1/(r - 2*M)**2`
- `1/(r*(r - 2*M))`

This generated `192` correction candidates. In
this run, `192` passed the strict prechecks before
the full residual, and `0` had exact-zero full residual.
The full Cartesian-product basis list is stored in the JSON artifact.

The gate then appends a bounded two-term correction screen:

```text
Psi = 1 - x + a**2 * (c1*basis_i(r,x) + c2*basis_j(r,x))
c1, c2 in {-1, 1}
```

Two-term angular factors:

- `x`
- `1-x**2`
- `x*(1-x**2)`
- `(1-x**2)**2`

Two-term radial factors:

- `1/r`
- `1/r**2`
- `1/(r - 2*M)`

This generated `264` two-term correction
candidates. In this run, `264` passed strict
prechecks before the full residual, and `0` had
exact-zero full residual.

The gate then runs a leading-order coefficient solve over the full one-term
basis instead of only trying fixed scalar coefficients:

```text
Psi = 1 - x + a**2 * sum_i c_i*basis_i(r,x)
coefficient of a**2 in full residual series through O(a**4)
M = 1
```

This screen produced `66` polynomial
equations in `48` unknown coefficients.  The
linear system matrix shape was `[66, 48]`, and SymPy
returned `EmptySet`.  It therefore generated
`0` additional candidates.


The gate also calibrates itself against the literature slow-rotation
split-monopole perturbative correction:

```text
Psi = 1 - x + a**2*x*(1-x**2)*R(r)
coefficient of a**2 in full residual series through O(a**4)
M = 1
```

Using `polylog(1, z) = -log(1 - z)`, the leading-order
residual simplifies to `0`.
The screen status is `passes_leading_order_anchor`. This is not added as a
finite-spin exact candidate: This is a perturbative O(a**2) literature anchor, not an exact finite-spin Kerr paper candidate.



## Criteria status matrix

This matrix prevents the negative artifact from being overread. It states
which parts of the paper criteria are implemented in this gate and which remain
future work before a positive solution claim.

| Criterion id | Criterion | Current status | Evidence |
| --- | --- | --- | --- |
| `literature_target` | Target equation and no-known-exact-solution boundary are cited from the literature. | implemented | Mahlmann et al. 2018 supplies the Kerr GSE target; Camilloni et al. 2020 supplies the no-known-exact-analytic extreme-Kerr boundary. |
| `engine_generation` | The engine supplies a bounded expression table rather than hand-picked paper rows. | implemented | 306 generated expressions loaded for the full target gate. |
| `strict_prechecks` | Candidate rows must depend on r, x, and a, be finite on rational safe points, meet the small-spin anchor, and not equal the known anchor. | implemented | Every assessment records strict_prechecks before the full nonlinear residual is evaluated. |
| `full_residual` | Candidate rows must pass the closed split-monopole nonlinear Kerr Grad-Shafranov residual. | implemented_negative | 0 candidates have exact-zero full residual; admitted_count=0. |
| `bounded_correction_screens` | The split-monopole anchor is tested with bounded one-term and two-term finite-spin correction grammars. | implemented_negative | 192 one-term and 264 two-term rows generated; 456 pass strict prechecks; 0 exact-zero residual rows. |
| `coefficient_solve` | The one-term basis is tested with arbitrary leading-order coefficients, not only sampled constants. | no_leading_order_solution | 66 equations, 48 unknowns, matrix [66, 48], linsolve=EmptySet. |
| `literature_perturbative_anchor` | The gate is calibrated against the known O(a**2) Blandford-Znajek split-monopole perturbative correction. | passes_leading_order_anchor | leading_residual_exact_zero=True; This is a perturbative O(a**2) literature anchor, not an exact finite-spin Kerr paper candidate. |
| `equivalence_filters` | Candidate rows must not be equivalent to known solutions under broader gauge, scaling, coordinate, or reparameterization transformations. | partial_not_sufficient_for_positive_claim | The current gate rejects exact known anchors and trivial equivalents only; broader equivalence filters remain future work. |
| `global_regularities` | Candidate rows must pass horizon, axis, and light-surface regularity checks for a positive paper claim. | not_implemented_for_positive_claim | The current gate checks symbolic finiteness and denominator safety on rational safe points; it does not prove global horizon/axis/light-surface regularity. |

## Source run

- Run id: `paper_repro_20260526_104304_2b096bdc`
- Database: `problems/kerr_magnetosphere/outputs/parallel_runs_paper_repro_20260526_104304_2b096bdc.db`
- Table: `expressions_paper_repro_20260526_104304_2b096bdc`
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
