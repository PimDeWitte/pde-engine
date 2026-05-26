# Force-free novel discoveries from pde-engine

This artifact records pde-engine discoveries that are novel relative to the
repository's seven registered force-free `known_solutions`. It is not a
literature-priority claim.

## Source run

- Run id: `paper_repro_20260526_071520_c7dcbaa2`
- Database: `problems/force_free/outputs/parallel_runs_paper_repro_20260526_071520_c7dcbaa2.db`
- Table: `expressions_paper_repro_20260526_071520_c7dcbaa2`
- Command: `python3 general_method_paper_reproduction.py --problem force_free --max-depth 2 --validators 1`
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
