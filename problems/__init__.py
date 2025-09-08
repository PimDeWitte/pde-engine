"""
Problem plugin interface and built-in problems.

This module defines a minimal plugin system to make the discovery/validation
engine agnostic to the specific PDE being solved. Each problem provides:

- Symbols and constants
- Primitives and operation sets
- A validator with a uniform interface
- Known solutions (optional)
- An output directory root where reports/databases are stored
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Tuple, Any, Optional, Callable

import sympy as sp

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from expression_operations import (
    UNARY_OPS,
    BINARY_OPS,
    SPECIAL_OPS,
    ALL_BINARY_OPS,
)


@dataclass
class ProblemSpec:
    """Specification container for a PDE discovery problem."""

    name: str
    slug: str

    # Symbols (coordinates) and constants
    symbols: Dict[str, sp.Symbol]
    constants: Dict[str, sp.Symbol]

    # Expression generation building blocks
    primitives: List[sp.Basic]
    unary_ops: Dict[str, Callable]
    binary_ops: Dict[str, Callable]
    special_ops: Dict[str, Callable]
    all_binary_ops: Dict[str, Callable]

    # Validator must expose: validate(u: sp.Basic, check_regularity: bool = True, fast_point_only: bool = False) -> Tuple[bool, str]
    validator: Any

    # Optional mapping of known solutions (stringified expression -> human name)
    known_solutions: Dict[str, str]

    # Root directory for outputs of this problem (reports, DBs, etc.)
    output_root: str

    def get_output_dir(self) -> str:
        os.makedirs(self.output_root, exist_ok=True)
        return self.output_root


def _create_force_free_problem() -> ProblemSpec:
    """Create force-free foliation problem specification."""
    from problems.force_free.validator import PreciseFoliationValidator
    
    rho = sp.Symbol('rho', real=True, positive=True)
    z = sp.Symbol('z', real=True)

    primitives: List[sp.Basic] = [
        rho,
        z,
        rho**2 + z**2,
        rho / z,
        sp.Integer(1),
    ]

    # Use problem-specific cache database
    cache_path = os.path.join('problems', 'force_free', 'outputs', 'validator_cache.db')
    validator = PreciseFoliationValidator(cache_db=cache_path, use_lean=True)

    known_solutions = {
        'rho**2': 'Vertical field',
        'rho**2*z': 'X-point',
        '1 - z/sqrt(rho**2 + z**2)': 'Radial',
        'rho**2/(rho**2 + z**2)**(3/2)': 'Dipolar',
        'sqrt(rho**2 + z**2) - z': 'Parabolic',
        'sqrt(z**2 + (rho - 1)**2) - sqrt(z**2 + (rho + 1)**2)': 'Hyperbolic',
        'rho**2*exp(-2*z)': 'Bent',
    }

    return ProblemSpec(
        name="Force-Free Foliations",
        slug="force_free",
        symbols={'rho': rho, 'z': z},
        constants={},
        primitives=primitives,
        unary_ops=UNARY_OPS,
        binary_ops=BINARY_OPS,
        special_ops=SPECIAL_OPS,
        all_binary_ops=ALL_BINARY_OPS,
        validator=validator,
        known_solutions=known_solutions,
        output_root=os.path.join('problems', 'force_free', 'outputs'),
    )


class _LegacyKerrValidatorToRemove:
    """
    Validator for the linear surrogate PDE modeling a Kerr magnetosphere:

        ∂r [ (G / (1 - x^2)) ∂r Ψ ] + ∂x [ (G / Δ) ∂x Ψ ] = 0

    with Δ = r^2 - 2 M r + a^2,
         G = 1 - 2 M r / (r^2 + a^2 x^2).

    Supports a fast point-evaluation mode and optional regularity checks at the
    axis (x = ±1) and at the horizon (r = r_+).
    """

    def __init__(self, r: sp.Symbol, x: sp.Symbol, M: sp.Symbol, a: sp.Symbol,
                 M_value: sp.Rational | int | float = sp.Integer(1),
                 a_value: sp.Rational | int | float = sp.Rational(1, 10)):
        self.r = r
        self.x = x
        self.M = M
        self.a = a
        self.M_value = M_value
        self.a_value = a_value

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
        ur = sp.diff(u, r)
        ux = sp.diff(u, x)
        term_r = sp.diff(G / (1 - x**2) * ur, r)
        term_x = sp.diff(G / Delta * ux, x)
        return sp.simplify(term_r + term_x)

    def _safe_eval(self, expr: sp.Basic, subs_dict: Dict[sp.Symbol, Any]) -> Optional[sp.Basic]:
        try:
            val = sp.simplify(expr.subs(subs_dict))
            if val.is_number:
                return val
            return sp.simplify(val)
        except Exception:
            return None

    def validate(self, u: sp.Basic, check_regularity: bool = True, fast_point_only: bool = False) -> Tuple[bool, str]:
        try:
            lhs = self._lhs(u)

            if not fast_point_only:
                if sp.simplify(lhs) == 0:
                    # Optional regularity checks
                    if check_regularity and not self._check_regularity(u):
                        return False, "Symbolically zero but fails regularity checks"
                    return True, "Valid solution (symbolically zero)"

            # Fast numeric/rational point checks
            Ms = self.M_value
            as_ = self.a_value
            subs_base = {self.M: Ms, self.a: as_}
            test_points = [
                {self.r: sp.Rational(5, 2), self.x: sp.Rational(3, 5)},
                {self.r: sp.Rational(7, 3), self.x: sp.Rational(1, 3)},
                {self.r: sp.Rational(5, 1), self.x: sp.Rational(-2, 5)},
            ]

            max_abs = 0.0
            for tp in test_points:
                val = self._safe_eval(lhs, {**subs_base, **tp})
                if val is None:
                    continue
                try:
                    fv = float(val)
                    max_abs = max(max_abs, abs(fv))
                    if val == 0:
                        continue
                except Exception:
                    # Non-float but exact zero check
                    if sp.simplify(val) != 0:
                        return False, f"Constraint not satisfied at test point {tp}"

            if max_abs > 1e-10:
                return False, f"Constraint not satisfied (max |LHS| ≈ {max_abs:.2e})"

            if check_regularity and not fast_point_only:
                if not self._check_regularity(u):
                    return False, "Satisfies PDE at points but fails regularity checks"

            return True, "Valid up to fast checks"

        except Exception as e:
            return False, f"Validation error: {e}"

    def _check_regularity(self, u: sp.Basic) -> bool:
        # Axis regularity: flux term in r should be finite as x -> ±1
        r, x = self.r, self.x
        Delta = self._delta()
        G = self._G()
        ur = sp.diff(u, r)
        ux = sp.diff(u, x)

        flux_r = sp.simplify(G / (1 - x**2) * ur)
        flux_x = sp.simplify(G / Delta * ux)

        try:
            lim1 = sp.limit(flux_r, x, 1)
            lim2 = sp.limit(flux_r, x, -1)
            if any(val in (sp.oo, -sp.oo, sp.zoo) for val in (lim1, lim2)):
                return False
        except Exception:
            # If limit computation fails, be conservative
            return False

        # Horizon regularity: flux_x finite as r -> r_+
        try:
            Ms = self.M_value
            as_ = self.a_value
            r_plus = Ms + sp.sqrt(Ms**2 - as_**2)
            lim_h = sp.limit(flux_x.subs({self.M: Ms, self.a: as_}), self.r, r_plus)
            if lim_h in (sp.oo, -sp.oo, sp.zoo):
                return False
        except Exception:
            return False

        return True

    def describe(self) -> Dict[str, str]:
        """Descriptor of the validator and PDE operator without heavy simplification.

        Builds a symbolic operator skeleton using SymPy Derivative nodes only (cheap to construct),
        with the exact coefficient functions G and Δ injected from this validator.
        """
        u = sp.Function('u')(self.r, self.x)
        G = self._G()
        Delta = self._delta()
        lhs = sp.Derivative(G / (1 - self.x**2) * sp.Derivative(u, self.r), self.r) \
            + sp.Derivative(G / Delta * sp.Derivative(u, self.x), self.x)
        return {
            'method_name': f"{self.__class__.__module__}.{self.__class__.__name__}.validate",
            'math_definition': str(lhs),
        }


def _create_kerr_magnetosphere_problem() -> ProblemSpec:
    """Create Kerr magnetosphere problem specification."""
    from problems.kerr_magnetosphere.validator import KerrMagnetosphereValidator
    
    r = sp.Symbol('r', real=True, positive=True)
    x = sp.Symbol('x', real=True)
    M = sp.Symbol('M', real=True, positive=True)
    a = sp.Symbol('a', real=True)

    Delta = r**2 - 2*M*r + a**2
    G = 1 - (2*M*r) / (r**2 + a**2 * x**2)

    primitives: List[sp.Basic] = [
        r,
        x,
        sp.Integer(1),
        sp.Rational(1, 3),  # <reason>chain: Expose exact one-third coefficient to allow Kerr small-spin families with 1/3 factors</reason>
        (1 - x),
        a**2,  # <reason>chain: Direct access to spin-squared factor a^2 for constructing analytic Kerr ansätze</reason>
        r**2 + a**2 * x**2,
        Delta,
        G,
    ]

    validator = KerrMagnetosphereValidator(r, x, M, a, M_value=sp.Integer(1), a_value=sp.Rational(1, 10))

    known_solutions = {
        '1 - x': 'Monopole (a -> 0 limit)'
    }

    return ProblemSpec(
        name="Kerr Magnetosphere (linear surrogate)",
        slug="kerr_magnetosphere",
        symbols={'r': r, 'x': x},
        constants={'M': M, 'a': a},
        primitives=primitives,
        unary_ops=UNARY_OPS,
        binary_ops=BINARY_OPS,
        special_ops=SPECIAL_OPS,
        all_binary_ops=ALL_BINARY_OPS,
        validator=validator,
        known_solutions=known_solutions,
        output_root=os.path.join('problems', 'kerr_magnetosphere', 'outputs'),
    )


def derive_small_spin_odes(M_value: int | float = 1):
    """
    Derive the O(a^2) small-spin correction system (★) projected onto P1 and P3.

    Returns a tuple (odes, context) where:
      - odes is a dict with keys 'f1', 'f3' mapping to sympy Eq for the ODEs
      - context contains the symbols and functions used
    """
    r = sp.Symbol('r', real=True, positive=True)
    x = sp.Symbol('x', real=True)
    M = sp.Integer(M_value) if isinstance(M_value, int) else sp.nsimplify(M_value)

    # Legendre polynomials
    P1 = x
    P3 = sp.Rational(1,2) * (5*x**3 - 3*x)

    f1 = sp.Function('f1')(r)
    f3 = sp.Function('f3')(r)

    U = f1*P1 + f3*P3

    # Left-hand side operator in (★)
    L_U = (1 - x**2) * sp.diff(U, x, 2) - r * (r - 2*M) * sp.diff(U, r, 2) - 2*M * sp.diff(U, r)

    # Right-hand side forcing in (★)
    rhs = 4*M * r**2 * (r - 2*M) * (x**3 - x)

    # Project onto P1 and P3 using L2 inner product over x in [-1,1]
    eq1 = sp.integrate((L_U - rhs) * P1, (x, -1, 1))
    eq3 = sp.integrate((L_U - rhs) * P3, (x, -1, 1))

    # Simplify and factor common pieces
    eq1_s = sp.simplify(sp.together(eq1))
    eq3_s = sp.simplify(sp.together(eq3))

    odes = {
        'f1': sp.Eq(eq1_s, 0),
        'f3': sp.Eq(eq3_s, 0),
    }

    context = {
        'r': r, 'x': x, 'M': M,
        'P1': P1, 'P3': P3,
        'f1': f1, 'f3': f3,
        'U': U,
    }

    return odes, context


def load_problem(name: str) -> ProblemSpec:
    key = (name or '').strip().lower()
    if key in ('force_free', 'forcefree', 'foliation', 'foliations'):
        return _create_force_free_problem()
    if key in ('kerr', 'kerr_magnetosphere', 'kerr-magnetosphere'):
        return _create_kerr_magnetosphere_problem()
    raise ValueError(f"Unknown problem '{name}'. Available: 'force_free', 'kerr_magnetosphere'")


__all__ = [
    'ProblemSpec',
    'load_problem',
]


