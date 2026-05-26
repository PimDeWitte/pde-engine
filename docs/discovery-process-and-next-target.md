# Discovery process record and next target

This document records what was actually done in this PR and why the force-free
candidate list should not be mistaken for the final research objective.

## What I actually did

The force-free work was a control problem. I did not perform a broad literature
search and then claim new physics. I used the existing Compere-Gralla-Lupsasca
force-free foliation problem to debug the pde-engine/Lagra interface and to
make the evidence path reproducible.

The exact source equation was the non-rotating stationary axisymmetric
force-free foliation determinant from Compere, Gralla, and Lupsasca,
"Force-Free Foliations," Phys. Rev. D 94, 124012 (2016), arXiv:1606.06727,
DOI:10.1103/PhysRevD.94.124012.

The implementation target was:

```text
det([[L_T(A), L_T(B)], [L_T^2(A), L_T^2(B)]]) = 0
A = u_rho_rho + u_z_z - u_rho/rho
B = u_rho**2 + u_z**2
T = u_z*d_rho - u_rho*d_z
```

The work proceeded in these steps:

1. Fixed the force-free validator cache key so validation context is included:
   expression, schema, `Omega`, `check_regularity`, `fast_point_only`, and
   `use_lean`.
2. Removed the broken fast-point placeholder derivative shortcut. Fast point
   mode now still builds the exact determinant and only skips the expensive
   full-plane simplification.
3. Fixed the parallel worker import path so spawned workers load
   `problems.load_problem` from this repo instead of a stale package namespace.
4. Stabilized the bounded reproduction smoke with atomic queue claims, writer
   flush fixes, child-process cleanup, and a per-expression validation timeout.
5. Removed the unnecessary Lean `mathlib` dependency for the tiny normalizer so
   `lake build` is a small local check instead of a ProofWidgets/mathlib build.
6. Ran a bounded engine pass through:

```sh
python3 tools/export_force_free_novel_discoveries.py \
  --run-engine --max-depth 2 --validators 1 \
  --timeout-s 90 --validation-timeout-s 3
```

7. Exported only rows that passed these post-run filters:

- row marked valid by the bounded force-free validator
- not flagged as one of the seven registered force-free known solutions
- expression contains both `rho` and `z`
- independently rebuilt determinant simplifies exactly to zero
- determinant at `(rho,z)=(4/5,6/7)` is exactly zero

## What I did not do

The force-free list is not a publication result by itself.

I did not establish any of the following:

- literature priority
- global physical regularity
- boundary-condition acceptability
- equivalence under foliation reparameterization
- inequivalence to known solution families
- a full depth-4 seven-solution pde-engine reproduction

The correct label for the six rows in `docs/force-free-novel-discoveries.*` is:
**process-control candidate foliations, novel only relative to the repo
registry**.

## Same criteria for a real paper target

For a paper, a candidate must satisfy stronger criteria than the force-free
control run:

1. The target equation and boundary conditions are cited from the literature.
2. The literature states, or strongly implies, that no exact analytic solution
   meeting the chosen criteria is known.
3. The engine generates the candidate from an explicit grammar, not manual
   algebra.
4. The candidate residual is independently rebuilt and simplified to exact zero.
5. The candidate survives rational point checks away from singular sets.
6. The candidate is nontrivial: it depends on the required variables and
   parameters.
7. The candidate is finite on the physical domain and passes the relevant
   horizon/axis/light-surface regularity checks.
8. The candidate meets the literature boundary/anchor condition.
9. The candidate is not equivalent to known solutions under allowed gauge,
   scaling, coordinate, or reparameterization transformations.
10. The full evidence is exported as JSON plus a human-readable paper artifact.

## Proposed next target

The next target should be the Kerr force-free magnetosphere / relativistic
Grad-Shafranov problem, not another registry-rich force-free foliation problem.

Reason:

- Mahlmann et al. describe the static, axisymmetric, force-free Kerr
  magnetosphere problem as relying on solutions of the relativistic
  Grad-Shafranov equation, and emphasize numerical solution methodologies for
  established setups including split-monopole, paraboloidal, black-hole-disk,
  and uniform configurations.
  Source: https://arxiv.org/abs/1802.00815
- The same paper notes that general solutions for arbitrarily large black-hole
  spin require numerical evaluation of the GSE, while small-spin analytic
  solutions are a special case.
  Source: https://academic.oup.com/mnras/article/477/3/3927/4963751
- Camilloni, Grignani, Harmark, Oliveri, and Orselli state that for extreme
  Kerr there is no known exact analytic force-free solution that is stationary,
  axisymmetric, and magnetically dominated; they proceed perturbatively away
  from the NHEK attractor.
  Source: https://arxiv.org/abs/2007.15665

This target matches the real goal: find an exact symbolic candidate where a
paper-level solution meeting the criteria is not already sitting in the known
solution registry.

## Current repo status for that target

The existing `kerr_magnetosphere` problem is a **linear surrogate**:

```text
partial_r[(G/(1-x^2)) partial_r Psi] +
partial_x[(G/Delta) partial_x Psi] = 0
Delta = r^2 - 2 M r + a^2
G = 1 - 2 M r/(r^2 + a^2 x^2)
```

It is useful as a harness, but it is not yet the full nonlinear Kerr
force-free Grad-Shafranov problem. Historical local outputs also contain many
degenerate expressions such as `1/(1 - 1)`, so the next step cannot be "run it
and trust the valid rows." The next step must tighten the target before search.

## Next implementation step

Before running a serious search, add a new target manifest and validator gate
for the Kerr paper target:

- source equation and citation metadata
- physical-domain assumptions (`r > r_+`, `|x| < 1`, admissible spin)
- denominator/singularity rejection before validation
- nontriviality requirements: dependence on `r`, `x`, and `a`
- small-spin or NHEK anchor condition from the target paper
- axis/horizon/light-surface regularity checks
- equivalence filters for constant shifts, scalings, and known anchors
- JSON exporter with a `no_candidate_yet` result if no row meets all criteria

Only after that should we let the engine search. A negative result under the
strict gate is useful; a degenerate "valid" row is not.

## First gate now added

This PR now includes the first machine-readable version of that gate:

- `tools/export_kerr_paper_target_gate.py`
- `tests/test_kerr_paper_target_gate.py`
- `docs/kerr-paper-target-gate.json`
- `docs/kerr-paper-target-gate.md`

The gate ran the current `kerr_magnetosphere` harness with:

```sh
python3 tools/export_kerr_paper_target_gate.py \
  --run-engine --max-depth 2 --validators 1 \
  --timeout-s 90 --validation-timeout-s 3
```

It generated `306` rows, completed `306` validations, found `0` valid rows in
the current linear surrogate run, and admitted `0` paper candidates. It also
probed `1 - x`, `x`, `1/(1 - 1)`, and `1 - x + a**2*r*x` to prove that the
gate rejects known anchors, undefined expressions, and anchor-like expressions
until the full nonlinear Kerr force-free Grad-Shafranov validator is
implemented.

The status is therefore **`no_candidate_yet`**. That is intentional. It is the
right artifact for the next paper program until the real target equation,
regularity conditions, and equivalence filters are encoded.
