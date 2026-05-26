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
- `tools/export_kerr_paper_target_gate.py` runs the stricter next-target gate
  for the closed split-monopole Kerr force-free / Grad-Shafranov paper program,
  appends bounded one-term and two-term finite-spin correction grammars, and
  runs a leading-order coefficient solve plus a literature slow-rotation anchor
  check before emitting a `no_candidate_yet` artifact instead of promoting the
  current linear surrogate.
- `tests/test_force_free_novel_discoveries.py` checks the independent
  determinant builder and the registry-novel candidate filter.
- `tests/test_kerr_paper_target_gate.py` checks that undefined expressions,
  known small-spin anchors, the Schwarzschild monopole limit, and an
  anchor-like nontrivial probe are handled by the full residual gate.
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

## Kerr paper-target gate

The next-target artifact is:

- `docs/kerr-paper-target-gate.json`
- `docs/kerr-paper-target-gate.md`

It ran the current `kerr_magnetosphere` harness with:

```sh
python3 tools/export_kerr_paper_target_gate.py \
  --run-engine --max-depth 2 --validators 1 \
  --timeout-s 90 --validation-timeout-s 3
```

The result is `no_candidate_yet`: `306` generated rows, `306` completed
validations, all `306` generated expressions scanned by the full target gate,
`192` bounded one-term finite-spin correction candidates, `264` bounded
two-term correction candidates, four additional probes, a leading-order
coefficient solve over the full one-term basis, and `0` admitted paper
candidates across `766` total candidate inputs. The one-term targeted grammar
is:

```text
Psi = 1 - x + a**2*c*basis(r,x)
c in {-1, -1/2, 1/2, 1}
basis = angular_factor * radial_factor
```

The angular factors are `1`, `x`, `x**2`, `1-x**2`, `x*(1-x**2)`,
`(1-x**2)**2`, `x**2*(1-x**2)`, and `x*(1-x**2)**2`. The radial factors are
`1/r`, `1/r**2`, `1/r**3`, `1/(r - 2*M)`, `1/(r - 2*M)**2`, and
`1/(r*(r - 2*M))`.

All `192` correction candidates pass the strict prechecks before the full
residual and then fail exact-zero residual. The two-term screen uses
`Psi = 1 - x + a**2*(c1*basis_i(r,x) + c2*basis_j(r,x))` with
`c1,c2 in {-1,1}`, a 12-function sub-basis, and produces `264` candidates.
All `264` pass the strict prechecks before the full residual and then fail
exact-zero residual. The gate also probes `1 - x`, `x`,
`1/(1 - 1)`, and `1 - x + a**2*r*x`. These are rejected respectively as known
anchors, undefined expressions, or finite-spin anchor-like expressions whose
full nonlinear Kerr split-monopole Grad-Shafranov residual is nonzero at
rational safe points.

The coefficient solve is stricter than the fixed correction screens:
`Psi = 1 - x + a**2*sum_i c_i*basis_i(r,x)` with `M = 1`, using the coefficient
of `a**2` in the full residual series through `O(a**4)`. It produces `66`
polynomial equations in `48` unknown coefficients, matrix shape `[66, 48]`,
and SymPy `linsolve` returns `EmptySet`. That rules out the entire one-term
basis as a leading-order correction family rather than only the sampled
coefficients.

The gate now also checks the known Blandford-Znajek slow-rotation perturbative
anchor from Tanabe-Nagataki and Pan-Yu:
`Psi = 1 - x + a**2*x*(1-x**2)*R(r)`, where `R(r)` contains the literature
logarithm and dilogarithm radial correction. After applying
`polylog(1,z) = -log(1-z)`, the coefficient of `a**2` in the closed residual
simplifies exactly to `0`, and all rational point checks are zero. The artifact
marks this as a perturbative `O(a**2)` anchor, not an exact finite-spin paper
candidate.

The artifact also carries a criteria-status matrix. It marks the literature
target, bounded engine generation, strict prechecks, nonlinear residual screen,
finite-spin correction screens, coefficient solve, and literature perturbative
anchor calibration as implemented. It marks broader equivalence filtering and
global horizon/axis/light-surface regularity as not sufficient for a positive
paper claim. This is deliberate: the current artifact is a negative gate and a
target packet, not a claim that all positive paper criteria have been
implemented.

## Literature sources

- Geoffrey Compere, Samuel E. Gralla, Alexandru Lupsasca, "Force-Free
  Foliations," Phys. Rev. D 94, 124012 (2016), arXiv:1606.06727,
  DOI:10.1103/PhysRevD.94.124012.
- arXiv: <https://arxiv.org/abs/1606.06727>
- DOI: <https://doi.org/10.1103/PhysRevD.94.124012>
- J. F. Mahlmann, P. Cerda-Duran, M. A. Aloy et al., "Numerically solving the
  relativistic Grad-Shafranov equation in Kerr spacetimes: Numerical
  techniques," MNRAS 477, 3927-3946 (2018), arXiv:1802.00815,
  DOI:10.1093/mnras/sty858.
- F. Camilloni, G. Grignani, T. Harmark, R. Oliveri, M. Orselli, "Moving away
  from the Near-Horizon Attractor of the Extreme Kerr Force-Free
  Magnetosphere," JCAP 10, 048 (2020), arXiv:2007.15665,
  DOI:10.1088/1475-7516/2020/10/048.

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
- `docs/kerr-paper-target-gate.json`
- `docs/kerr-paper-target-gate.md`
- `tools/render_lagra_force_free_adapter_paper.py`
- `tools/export_force_free_novel_discoveries.py`
- `tools/export_kerr_paper_target_gate.py`

These are pde-engine-facing paper artifacts generated on 2026-05-26 in the
visual style of `/Users/p/Downloads/LagraPaperV1 (1).pdf`. The paper maps the
diagram style in that reference paper onto the pde-engine force-free adapter:
pipeline, vocabulary chips, graph anatomy, data flow, exploded kernel, verifier
triangle, derivation-style cache regression, engine-discovery export, and the
Kerr paper-target gate. Each diagram states what enters Lagra and what
measurement comes out.

Artifact hashes:

```text
d992821aae6d3bd0d56cf565ef972fbf1d2327217d8ad3e149db1bc1024e4fba  docs/lagra-force-free-adapter-paper.html
9c6587e617de2e5e189d4fc3f174ffa46fa47f26b0a342c274b9276d8c364cb1  docs/lagra-force-free-adapter-paper.pdf
dee917615ba99060e29970219c4f8726929f3c296ab6d5acde0da819b993c89a  tools/render_lagra_force_free_adapter_paper.py
f67202596b8fe94c85b6ca9ebba02267a7891ef6929add75bf0a0b5f73f67d6c  docs/force-free-novel-discoveries.json
b6073e74d079d19d7183101f5a21770064e0f3ea450fd257c803afc9e2cd301c  docs/force-free-novel-discoveries.md
cd5999772df7b0e2a06d10831f9863ddb51ea5ecca5276e314eedc6d76517b41  tools/export_force_free_novel_discoveries.py
024919b59a15aca96c17da99f1522097633715e072a88f9bd4e494ef463bce0b  docs/discovery-process-and-next-target.md
3c03c4114ca86268bd213aa71ee0e41be023bdd73201986a5b4bcfb1154679f8  docs/kerr-paper-target-gate.json
ba389f6999d42c79e6db15956d3a25ab6a5e4dba3324714cbe6a9a2a4134fb3f  docs/kerr-paper-target-gate.md
b06927bdb7165c5b076d5440189b3fb8dadfd35b807f5bee0827804b429dd8be  tools/export_kerr_paper_target_gate.py
```

## Reproduction commands

From this pde-engine checkout:

```sh
python3 -m py_compile \
  general_method_paper_reproduction.py \
  problems/force_free/validator.py \
  tests/test_force_free_validator_cache.py \
  tools/export_force_free_novel_discoveries.py \
  tools/export_kerr_paper_target_gate.py \
  tests/test_force_free_novel_discoveries.py \
  tests/test_kerr_paper_target_gate.py

python3 -m unittest \
  tests/test_force_free_validator_cache.py \
  tests/test_force_free_novel_discoveries.py \
  tests/test_kerr_paper_target_gate.py

python3 tools/export_force_free_novel_discoveries.py \
  --run-engine --max-depth 2 --validators 1 \
  --timeout-s 90 --validation-timeout-s 3

python3 tools/export_kerr_paper_target_gate.py \
  --run-engine --max-depth 2 --validators 1 \
  --timeout-s 90 --validation-timeout-s 3

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
