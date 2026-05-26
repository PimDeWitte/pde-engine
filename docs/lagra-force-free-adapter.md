# Lagra force-free adapter note

This change records the pde-engine side of Lagra's external force-free
adapter gate.

## What changed in pde-engine

- `problems/force_free/validator.py` now keys validator-cache entries by the
  expression and the validation context: schema, `Omega`, `check_regularity`,
  `fast_point_only`, and `use_lean`.
- The fast point mode no longer uses a second-derivative placeholder shortcut.
  It builds the exact symbolic determinant and stops after the exact paper
  point check. The determinant uses `L_T(A)`, `L_T(B)`, `L_T^2(A)`, and
  `L_T^2(B)`, so the shortcut was under-specified.
- `general_method_paper_reproduction.py` now imports `problems.load_problem`
  from the repository package namespace inside spawned worker processes, with
  a local fallback import for `PreciseFoliationValidator`.
- `tests/test_force_free_validator_cache.py` proves the cache separates
  rotation context: `rho**2*exp(-2*z)` is valid at `Omega=0`, invalid at
  `Omega=1`, and those two validations receive different cache hashes even
  when they share the same cache database.
- `tools/export_force_free_novel_discoveries.py` runs a bounded force-free
  engine pass, mines valid non-registered rows, independently rebuilds the
  determinant, and emits JSON/Markdown process-control artifacts.
- `tests/test_force_free_novel_discoveries.py` checks the independent
  determinant builder and the registry-novel candidate filter.
- `docs/discovery-process-and-next-target.md` records the full process,
  including the search bounds, post-hoc filters, missing criteria, and the next
  Kerr-force-free target.

## How Lagra used pde-engine

Lagra used pde-engine as an external operator and registry boundary, not as an
internal source of truth.

The literature target is Geoffrey Compere, Samuel E. Gralla, and Alexandru
Lupsasca, "Force-Free Foliations," Phys. Rev. D 94, 124012 (2016),
arXiv:1606.06727, DOI:10.1103/PhysRevD.94.124012. The implemented equation is
the stationary axisymmetric non-rotating force-free foliation determinant from
Eq. 2.14 / Section 2.4 of that paper.

The Lagra gate
`experiments/proof-search/pde_engine_force_free_point_gate.py` imports the
local pde-engine force-free symbols and loads the seven known force-free
foliations from the pde-engine problem registry. Lagra then evaluates the
force-free determinant

```text
det([[L_T(A), L_T(B)], [L_T^2(A), L_T^2(B)]])
```

with

```text
A = u_rho_rho + u_z_z - u_rho/rho
B = u_rho**2 + u_z**2
T = u_z*d_rho - u_rho*d_z
```

with its own exact SymPy adapter and records the result into its proof-search
artifact catalog. The current gate result is:

- `7/7` registry solutions simplify to exact symbolic determinant `0`.
- `2/2` false controls simplify to exact nonzero determinants:
  `rho*z -> 16*rho*z` and `exp(rho*z) -> 16*rho*z*exp(6*rho*z)`.
- The same expressions pass three high-precision rational checkpoint tests:
  known solutions remain below `1e-60`, while false controls exceed `1e-6`.
- The upstream validator accepts three bounded non-rotating controls.
- The rotation-context regression proves `rho**2*exp(-2*z)` is valid at
  `Omega=0`, invalid at `Omega=1`, and has rotating determinant
  `2048*rho**9*exp(-12*z)`.

Lagra binds those results into its `.lagra` proof-search constants and checks
them through publication-readiness and integrity gates before the result is
allowed into the rendered paper.

## Engine process-control candidates

The bounded depth-2 engine run now produces a committed process-control
artifact:

- `docs/force-free-novel-discoveries.json`
- `docs/force-free-novel-discoveries.md`

The exporter reruns:

```sh
python3 tools/export_force_free_novel_discoveries.py \
  --run-engine --max-depth 2 --validators 1
```

It then filters for valid rows not registered as one of the seven force-free
known solutions, requires both `rho` and `z` to appear, and rebuilds the
non-rotating force-free determinant independently of the validator cache.

The current exported candidates are:

| Expression | Candidate class | Independent witness |
| --- | --- | --- |
| `rho + z` | linear two-coordinate foliation | `det M = 0` |
| `rho/z` | scale-invariant ratio foliation | `det M = 0` |
| `rho**2 + z**2` | quadratic radius foliation | `det M = 0` |
| `rho**2 + z` | parabolic polynomial foliation | `det M = 0` |
| `-rho**2 + z**2 + 1` | hyperbolic quadratic foliation | `det M = 0` |
| `rho/(1 - z)` | rational geometric-sum foliation | `det M = 0` |

Novelty is scoped to the repository registry: these rows are not identical to
the seven registered `known_solutions`. This is not a literature-priority
claim.

The full process record is in
`docs/discovery-process-and-next-target.md`. The important correction is that
this force-free list is not the final research target. It is the control
problem used to debug and document the pde-engine/Lagra evidence path.

## Literature sources

- Geoffrey Compere, Samuel E. Gralla, Alexandru Lupsasca, "Force-Free
  Foliations," Phys. Rev. D 94, 124012 (2016), arXiv:1606.06727,
  DOI:10.1103/PhysRevD.94.124012.
- arXiv: <https://arxiv.org/abs/1606.06727>
- DOI: <https://doi.org/10.1103/PhysRevD.94.124012>

## Claim boundary

This is a positive pde-engine adapter, exact symbolic determinant check, point
check, and validator-cache regression. It is not a Lean proof. It is not a
full pde-engine reproduction claim. The new candidate artifact is a bounded
engine process-control result, not a claim of literature priority.

The Lagra health gate still records the broader pde-engine reproduction path
as not fully reproduced on the local checkout used for this paper. After the
second patch, the bounded depth-2 smoke no longer hits the old
`PreciseFoliationValidator` worker `NameError`, no longer fails the Lean build,
and no longer times out under the 20 second health wrapper. It exits with
return code `0`, generates `112` expressions, reports `76` valid rows, and finds
`2` known vertical canonical forms. It also exports six independently rechecked
non-registered process-control candidates. That is a process-clean smoke and a
bounded evidence-path artifact, not a full seven-solution Compere reproduction
and not the paper target.

## Included paper artifacts

- `docs/lagra-force-free-adapter-paper.pdf`
- `docs/lagra-force-free-adapter-paper.html`
- `docs/force-free-novel-discoveries.json`
- `docs/force-free-novel-discoveries.md`
- `docs/discovery-process-and-next-target.md`
- `tools/render_lagra_force_free_adapter_paper.py`
- `tools/export_force_free_novel_discoveries.py`

These are pde-engine-facing paper artifacts generated on 2026-05-26 in the
visual style of `/Users/p/Downloads/LagraPaperV1 (1).pdf`. The paper maps the
diagram style in that reference paper onto the pde-engine force-free adapter:
pipeline, vocabulary chips, graph anatomy, data flow, exploded kernel, verifier
triangle, derivation-style cache regression, and engine-discovery export. Each
diagram states what enters Lagra and what measurement comes out.

Artifact hashes:

```text
2da86b8029a0ad252c4ba8907b7e10019df42b8b949721a281b0f1db5cde5995  docs/lagra-force-free-adapter-paper.html
2555cb05126b691127c3d5db1900a94409ce02d86a78c23f61e81e86c8c0bf52  docs/lagra-force-free-adapter-paper.pdf
7a9c0cc5219b697843e6c80d9d263cc55537b3fefad34d4cd553c5b8d1887d2d  tools/render_lagra_force_free_adapter_paper.py
f67202596b8fe94c85b6ca9ebba02267a7891ef6929add75bf0a0b5f73f67d6c  docs/force-free-novel-discoveries.json
b6073e74d079d19d7183101f5a21770064e0f3ea450fd257c803afc9e2cd301c  docs/force-free-novel-discoveries.md
cd5999772df7b0e2a06d10831f9863ddb51ea5ecca5276e314eedc6d76517b41  tools/export_force_free_novel_discoveries.py
53a95f937062e0817358e5fbb0cb00ba4183b7eafd376fdd7175c6d7e3a0d334  docs/discovery-process-and-next-target.md
7feafdfd3e9dee42cbeed306ca2bf62a55d49a1659c23d662de360a83da6f44b  pde_engine_force_free_point_gate.json
e5a4c2b7ce63ba4044dcdfa5f0aa58758b3e2e5db92e12b8936e01770d3f99a1  pde_engine_reproduction_health_gate.json
```

## Reproduction commands

From this pde-engine checkout:

```sh
python3 -m py_compile \
  general_method_paper_reproduction.py \
  problems/force_free/validator.py \
  tests/test_force_free_validator_cache.py \
  tools/export_force_free_novel_discoveries.py \
  tests/test_force_free_novel_discoveries.py

python3 -m unittest \
  tests/test_force_free_validator_cache.py \
  tests/test_force_free_novel_discoveries.py

python3 tools/export_force_free_novel_discoveries.py \
  --run-engine --max-depth 2 --validators 1

python3 tools/render_lagra_force_free_adapter_paper.py
```

From the Lagra checkout that produced the included paper:

```sh
python3 experiments/proof-search/pde_engine_force_free_point_gate.py \
  --pde-engine-root /Users/p/dev/pde-engine

python3 experiments/proof-search/pde_engine_reproduction_health_gate.py \
  --pde-engine-root /Users/p/dev/pde-engine \
  --timeout-seconds 20

python3 experiments/proof-search/render_yukawa_finite_range_paper.py
python3 experiments/proof-search/publication_readiness_audit.py
python3 experiments/proof-search/proof_search_integrity_gate.py
```
