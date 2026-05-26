# Force-free novel discoveries from pde-engine

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

- Run id: `paper_repro_20260526_091146_fa5af2b9`
- Database: `problems/force_free/outputs/parallel_runs_paper_repro_20260526_091146_fa5af2b9.db`
- Table: `expressions_paper_repro_20260526_091146_fa5af2b9`
- Command: `python3 general_method_paper_reproduction.py --problem force_free --max-depth 2 --validators 1`
- Bounds: `max_depth=2`,
  `validators=1`,
  `per_expression_validation_timeout_s=3.0`
- Total generated: `112`
- Completed validations: `112`
- Valid rows: `76`

## Independently verified candidates

Selection policy: valid non-paper depth-2 rows, both rho and z present, independently recomputed determinant simplifies exactly to zero.

| Expression | Candidate class | Engine row | det M simplified | det M at (4/5, 6/7) |
| --- | --- | ---: | --- | --- |
| `rho + z` | linear two-coordinate foliation | 41 | `0` | `0` |
| `rho/z` | scale-invariant ratio foliation | 4 | `0` | `0` |
| `rho**2 + z**2` | quadratic radius foliation | 3 | `0` | `0` |
| `rho**2 + z` | parabolic polynomial foliation | 77 | `0` | `0` |
| `-rho**2 + z**2 + 1` | hyperbolic quadratic foliation | 108 | `0` | `0` |
| `rho/(1 - z)` | rational geometric-sum foliation | 44 | `0` | `0` |

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
