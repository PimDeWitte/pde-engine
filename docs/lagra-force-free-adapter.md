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
  determinant, and emits JSON/Markdown discovery artifacts.
- `tests/test_force_free_novel_discoveries.py` checks the independent
  determinant builder and the registry-novel candidate filter.

## How Lagra used pde-engine

Lagra used pde-engine as an external operator and registry boundary, not as an
internal source of truth.

The Lagra gate
`experiments/proof-search/pde_engine_force_free_point_gate.py` imports the
local pde-engine force-free symbols and loads the seven known force-free
foliations from the pde-engine problem registry. Lagra then evaluates the
force-free determinant

```text
det([[L_T(A), L_T(B)], [L_T^2(A), L_T^2(B)]])
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

## Engine-discovered candidates

The bounded depth-2 engine run now produces a committed discovery artifact:

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

## Claim boundary

This is a positive pde-engine adapter, exact symbolic determinant check, point
check, and validator-cache regression. It is not a Lean proof. It is not a
full pde-engine reproduction claim. The new candidate artifact is a bounded
engine-discovery result, not a claim of literature priority.

The Lagra health gate still records the broader pde-engine reproduction path
as not fully reproduced on the local checkout used for this paper. After the
second patch, the bounded depth-2 smoke no longer hits the old
`PreciseFoliationValidator` worker `NameError`, no longer fails the Lean build,
and no longer times out under the 20 second health wrapper. It exits with
return code `0`, generates `112` expressions, reports `76` valid rows, and finds
`2` known vertical canonical forms. It also exports six independently rechecked
non-registered candidates. That is a process-clean smoke and a bounded discovery
artifact, not a full seven-solution Compere reproduction.

## Included paper artifacts

- `docs/lagra-force-free-adapter-paper.pdf`
- `docs/lagra-force-free-adapter-paper.html`
- `docs/force-free-novel-discoveries.json`
- `docs/force-free-novel-discoveries.md`
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
276d62d04c25163713d93a6cf2cdee6411d75c8650f617a4d55c19cd7e208369  docs/lagra-force-free-adapter-paper.html
41189d9b2a4000dafdf765cafa34e18264161642499b44765d7e5bdfde8b1a02  docs/lagra-force-free-adapter-paper.pdf
80988342b7e719212599941e199a1a5d010937ab363b36d2fb2d240bbc991965  tools/render_lagra_force_free_adapter_paper.py
7554ceef41776dfab13151ad4f891b8d8fee464aceec673479481d0db4a21488  docs/force-free-novel-discoveries.json
b038e4da247b6d13c2aaa0565157bd6f0f967821a49c54585879088b4afed426  docs/force-free-novel-discoveries.md
70975317bf0606f66c9e9b2c6c97f8dac3ff65623a491e08095174a1ecc6e7a7  tools/export_force_free_novel_discoveries.py
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
