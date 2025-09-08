"""
Kerr magnetosphere linear surrogate PDE validator.

PDE (divergence form):
  ∂r[ (G/(1-x^2)) ∂r u ] + ∂x[ (G/Δ) ∂x u ] = 0,
  with Δ = r^2 - 2Mr + a^2,
       G = 1 - (2Mr)/(r^2 + a^2 x^2).

This validator avoids heavy simplification; optionally uses Lean normalizer for
symbolic zero detection on small expressions and otherwise falls back to fast
point checks.
"""

from __future__ import annotations

from typing import Any, Dict, Tuple, Optional
import sympy as sp

try:
    # Optional Lean integration for lightweight symbolic zero detection
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    from lean_normalizer.lean_bridge_fixed import LeanNormalizer
except Exception:  # pragma: no cover
    LeanNormalizer = None  # type: ignore


class KerrMagnetosphereValidator:
    def __init__(
        self,
        r: sp.Symbol,
        x: sp.Symbol,
        M: sp.Symbol,
        a: sp.Symbol,
        M_value: Any = sp.Integer(1),
        a_value: Any = sp.Rational(1, 10),
        use_lean: bool = True,
        lean_det_str_max_len: int = 12000,
        require_monopole_extension: bool = True,
        monopole_target: str = '1-x',  # '1-x' | 'x' | 'either'
        allow_normalization: bool = False,
        strict_sympy_check: bool = True,
        exclude_constants: bool = True,
    ) -> None:
        self.r = r
        self.x = x
        self.M = M
        self.a = a
        self.M_value = M_value
        self.a_value = a_value
        self.use_lean = use_lean and (LeanNormalizer is not None)
        self.lean_det_str_max_len = lean_det_str_max_len
        self.require_monopole_extension = require_monopole_extension
        self.monopole_target = monopole_target
        self.allow_normalization = allow_normalization
        self.strict_sympy_check = strict_sympy_check
        self.exclude_constants = exclude_constants
        self._lean = None  # type: ignore[var-annotated]
        if self.use_lean:
            try:
                self._lean = LeanNormalizer()
            except Exception:
                self._lean = None
                self.use_lean = False
        # Simple cache for residual exact-zero verdicts by expression string
        self._residual_zero_cache: dict[str, bool] = {}

    def _delta(self) -> sp.Basic:
        r, M, a = self.r, self.M, self.a
        return r**2 - 2*M*r + a**2

    def _G(self) -> sp.Basic:
        r, x, M, a = self.r, self.x, self.M, self.a
        return 1 - (2*M*r) / (r**2 + a**2 * x**2)

    def _lhs(self, u: sp.Basic) -> sp.Basic:
        r, x = self.r, self.x
        Delta = self._delta()
        G = self._G()
        try:
            ur = sp.diff(u, r)
            ux = sp.diff(u, x)
        except Exception:
            # If u is not in terms of (r,x), try to coerce by replacing symbols named 'r','x'
            ur = sp.diff(u.subs({sp.Symbol('r'): r, sp.Symbol('x'): x}), r)
            ux = sp.diff(u.subs({sp.Symbol('r'): r, sp.Symbol('x'): x}), x)
        # Avoid simplify/expand; keep as direct derivatives
        term_r = sp.diff(G / (1 - x**2) * ur, r)
        term_x = sp.diff(G / Delta * ux, x)
        return term_r + term_x

    def _finite_classical(self, expr: sp.Basic) -> bool:
        """Check that the expression is finite (no infinities/NaNs) and evaluates
        to finite values at a few safe points away from singular sets (Δ=0, x=±1).
        """
        try:
            e = sp.simplify(expr)
        except Exception:
            e = expr
        try:
            if e.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
                return False
        except Exception:
            return False
        tests = [
            {self.M: sp.Integer(1), self.a: sp.Rational(3, 5), self.r: sp.Rational(7, 3), self.x: sp.Rational(1, 3)},
            {self.M: sp.Integer(1), self.a: sp.Rational(4, 5), self.r: sp.Integer(3), self.x: -sp.Rational(2, 5)},
        ]
        for s in tests:
            try:
                val = sp.simplify(e.subs(s))
                if val.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
                    return False
            except Exception:
                return False
        return True

    def _is_nonconstant(self, u: sp.Basic) -> bool:
        try:
            ur = sp.simplify(sp.diff(u, self.r))
            ux = sp.simplify(sp.diff(u, self.x))
            return not (ur == 0 and ux == 0)
        except Exception:
            return True

    def _is_monopole_extension(self, Psi: sp.Basic) -> bool:
        """Require Ψ → 1 - x (default) as a → 0. If allow_normalization=True, also
        accept a finite constant offset at a=0. If monopole_target='either', also
        accept x.
        """
        targets: list[sp.Basic] = []
        if self.monopole_target in ('1-x', 'either'):
            targets.append(1 - self.x)
        if self.monopole_target in ('x', 'either'):
            targets.append(self.x)
        for tgt in targets:
            try:
                diff = sp.simplify(Psi - tgt)
            except Exception:
                diff = Psi - tgt
            try:
                L = sp.simplify(sp.limit(diff, self.a, 0))
            except Exception:
                try:
                    L = sp.simplify(diff.subs(self.a, 0))
                except Exception:
                    continue
            try:
                if L == 0:
                    return True
                if self.allow_normalization and not L.has(sp.oo, sp.zoo, sp.nan):
                    # Accept finite constants; reject symbolic non-constants
                    if L.free_symbols.issubset({self.M}):
                        # permit dependence on M only as an overall constant
                        return True
                    if L.is_number:
                        return True
            except Exception:
                continue
        return False

    def _fast_point_check(self, expr: sp.Basic) -> Tuple[bool, str]:
        Ms = self.M_value
        as_ = self.a_value
        subs_base = {self.M: Ms, self.a: as_}
        test_points = [
            {self.r: sp.Rational(5, 2), self.x: sp.Rational(3, 5)},
            {self.r: sp.Rational(7, 3), self.x: sp.Rational(1, 3)},
            {self.r: sp.Rational(5, 1), self.x: sp.Rational(-2, 5)},
        ]
        max_abs = 0.0
        successes = 0
        for tp in test_points:
            try:
                val = expr.subs({**subs_base, **tp})
                # Evaluate numerically with decent precision if symbolic
                val_num = sp.N(val, 40)
                if val_num.is_real is False and val_num.is_real is not None:
                    return False, "Invalid (non-real at test point)"
                fv = float(val_num)
                if fv != fv:  # NaN check
                    return False, "Invalid (NaN at test point)"
                max_abs = max(max_abs, abs(fv))
                successes += 1
            except Exception:
                continue
        if successes == 0:
            return False, "Indeterminate (no evaluable test points)"
        if max_abs < 1e-10:
            return True, "Valid (point checks ≈ 0)"
        return False, f"Invalid (point checks ≈ {max_abs:.2e})"

    def _lean_normalized(self, expr: sp.Basic) -> Optional[str]:
        if not self.use_lean or self._lean is None:
            return None
        s = str(expr)
        if len(s) > self.lean_det_str_max_len:
            return None
        try:
            res = self._lean.normalize_batch([(s, 0)])
            if res and isinstance(res, list):
                normalized = res[0].get('normalized')
                if isinstance(normalized, str):
                    return normalized
        except Exception:
            return None
        return None

    def validate(
        self,
        u: sp.Basic,
        check_regularity: bool = True,
        fast_point_only: bool = False,
        *,
        lean_first: bool = True,
        defer_heavy_checks: bool = True,
        enforce_anchor: Optional[bool] = None,
    ) -> Tuple[bool, str]:
        """Validator for exact solutions to the PDE.

        Modes:
        - lean_first=True, defer_heavy_checks=True (fast path): prove lhs==0 via Lean (or SymPy) first;
          skip finiteness/regularity/anchor checks unless exact zero is established. This is ideal for
          discovery since almost all candidates are rejected quickly.
        - defer_heavy_checks=False: after exact zero is established, also check finiteness, regularity,
          and small-spin anchor (controlled by enforce_anchor or self.require_monopole_extension).
        """
        try:
            # Quick constant elimination (very cheap): drop if u simplifies to a constant
            if self.exclude_constants:
                try:
                    us = sp.simplify(u)
                except Exception:
                    us = u
                try:
                    if not (us.has(self.r) or us.has(self.x)):
                        return False, "Trivial constant solution excluded"
                except Exception:
                    pass

            lhs = self._lhs(u)
            # Fast exact check: Lean first (if available), else SymPy simplify
            lean_zero = False
            sympy_zero = False
            normalized = None

            # Helper: short residual representation without heavy simplification
            def _short_residual_repr(expr: sp.Basic) -> str:
                try:
                    # Replace abstract derivatives wrt placeholders with their symbolic names for compactness
                    s_expr = expr.replace(lambda e: isinstance(e, sp.Derivative), lambda e: sp.Symbol('d'))
                    num, den = sp.as_numer_denom(s_expr)
                    s = f"{sp.sstr(num)}/{sp.sstr(den)}"
                    return s[:240]
                except Exception:
                    try:
                        s = sp.sstr(expr)
                        return s[:240]
                    except Exception:
                        return '<residual-unavailable>'

            # Quick numeric pre-filter on lhs to avoid expensive Lean/SymPy for obvious non-solutions
            try:
                ok_fast, _ = self._fast_point_check(lhs)
                if not ok_fast:
                    # Attach a short residual representation (no simplify)
                    residual_repr = _short_residual_repr(lhs)
                    return False, f"PDE residual != 0 (fast point check) | residual: {residual_repr[:240]}"
            except Exception:
                pass

            # Cache check
            try:
                key = str(u)
                if key in self._residual_zero_cache:
                    cached = self._residual_zero_cache[key]
                    if not cached:
                        return False, "PDE residual != 0 (cached)"
            except Exception:
                pass

            if lean_first and self.use_lean and self._lean is not None:
                normalized = self._lean_normalized(lhs)
                if isinstance(normalized, str) and normalized.strip() == '0':
                    lean_zero = True
            
            if not lean_zero and self.strict_sympy_check:
                try:
                    # Use a cheaper check first: expand+together+cancel before full simplify
                    lhs_q = sp.together(sp.cancel(lhs))
                    sympy_zero = (lhs_q == 0) or (sp.simplify(lhs_q) == 0)
                except Exception:
                    sympy_zero = False

            # Record evidence
            try:
                lhs_str = str(lhs)
            except Exception:
                lhs_str = '<lhs-string-error>'
            self._last_evidence = {
                'lhs_string': lhs_str if len(lhs_str) <= 4000 else lhs_str[:4000] + '...truncated...',
                'lean_normalized': normalized,
                'sympy_simplified_is_zero': bool(sympy_zero),
                'params': {'M': str(self.M_value), 'a': str(self.a_value)},
            }

            if not (lean_zero or sympy_zero):
                # Provide short residual representation to aid debugging (no simplify)
                residual_repr = _short_residual_repr(lhs)
                try:
                    self._residual_zero_cache[str(u)] = False
                except Exception:
                    pass
                return False, f"PDE residual != 0 | residual: {residual_repr[:240]}"

            # Exact zero established. Optionally defer heavy checks for speed.
            if defer_heavy_checks:
                try:
                    self._residual_zero_cache[str(u)] = True
                except Exception:
                    pass
                return True, "Valid (exact zero; heavy checks deferred)"

            # Heavy checks: constants, finiteness, regularity, anchor
            if self.exclude_constants and not self._is_nonconstant(u):
                return False, "Trivial constant solution excluded"

            if not self._finite_classical(u):
                return False, "non-finite"

            if not self._finite_classical(lhs):
                return False, "residual non-finite"

            if check_regularity and not self._check_regularity(u):
                return False, "Symbolic zero but fails regularity checks"

            must_anchor = self.require_monopole_extension if enforce_anchor is None else bool(enforce_anchor)
            if must_anchor and not self._is_monopole_extension(u):
                return False, "fails a->0 monopole anchor"

            return True, "valid"

        except Exception as e:
            return False, f"Validation error: {e}"

    def _check_regularity(self, u: sp.Basic) -> bool:
        r, x = self.r, self.x
        Delta = self._delta()
        G = self._G()
        try:
            lim1 = sp.limit(G / (1 - x**2) * sp.diff(u, r), x, 1)
            lim2 = sp.limit(G / (1 - x**2) * sp.diff(u, r), x, -1)
            if any(val in (sp.oo, -sp.oo, sp.zoo) for val in (lim1, lim2)):
                return False
        except Exception:
            return False
        try:
            Ms = self.M_value
            as_ = self.a_value
            r_plus = Ms + sp.sqrt(Ms**2 - as_**2)
            lim_h = sp.limit((G / Delta).subs({self.M: Ms, self.a: as_}) * sp.diff(u, x), r, r_plus)
            if lim_h in (sp.oo, -sp.oo, sp.zoo):
                return False
        except Exception:
            return False
        return True

    def describe(self) -> Dict[str, str]:
        u = sp.Function('u')(self.r, self.x)
        G = self._G()
        Delta = self._delta()
        lhs = sp.Derivative(G / (1 - self.x**2) * sp.Derivative(u, self.r), self.r) \
            + sp.Derivative(G / Delta * sp.Derivative(u, self.x), self.x)
        return {
            'method_name': f"{self.__class__.__module__}.{self.__class__.__name__}.validate",
            'math_definition': str(lhs),
        }

    def last_evidence(self) -> Dict[str, Any]:  # type: ignore[name-defined]
        return getattr(self, '_last_evidence', {})


