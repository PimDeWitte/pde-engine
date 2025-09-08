"""
Precise implementation of the force-free foliation constraint from Compère et al.

This module implements the exact constraint from equation 2.14 of the paper:
det[L_T(A) L_T(B); L²_T(A) L²_T(B)] = 0

where:
- A = u_ρρ + u_zz - u_ρ/ρ  
- B = u_ρ² + u_z²
- T = u_z ∂_ρ - u_ρ ∂_z (tangent vector field)
- L_T is the Lie derivative along T
"""

import sympy as sp
from sympy import symbols, simplify, Matrix, det, sqrt, exp, log, expand
from typing import Tuple, Dict, Optional, Any
import sqlite3
import hashlib
import logging
import time

logger = logging.getLogger(__name__)
import json
import os

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from lean_normalizer.lean_bridge_fixed import LeanNormalizer


class PreciseFoliationValidator:
    """
    Validates expressions against the exact force-free foliation constraint.
    Includes SQLite caching to avoid redundant computations.
    """
    
    def __init__(self, cache_db: str = None, use_lean: bool = True, Omega: Any = 0):
        """
        Initialize validator with symbols and cache.
        
        Args:
            cache_db: Path to SQLite database for caching results
            use_lean: Whether to use Lean normalizer for simplification
            Omega: Field line angular velocity (can be a constant or a function of u)
        """
        self.rho = symbols('rho', real=True, positive=True)
        self.z = symbols('z', real=True)
        self.Omega = Omega
        # Use problem-specific cache path if not provided
        if cache_db is None:
            cache_db = os.path.join(os.path.dirname(__file__), 'outputs', 'validator_cache.db')
        self.cache_db = cache_db
        self.use_lean = use_lean
        self._init_cache()
        
        # Try to initialize Lean normalizer if requested
        if use_lean:
            try:
                self.lean_normalizer = LeanNormalizer()
                logger.info("Using Lean normalizer for simplification")
            except Exception as e:
                logger.warning(f"Could not load Lean normalizer: {e}")
                self.use_lean = False
                self.lean_normalizer = None
        else:
            self.lean_normalizer = None

    # -------------------- Fast exact pointwise AD (up to 2nd order) --------------------
    def _ad_eval_point(self, expr: sp.Basic, rho0: sp.Basic, z0: sp.Basic) -> Tuple[sp.Basic, sp.Basic, sp.Basic, sp.Basic, sp.Basic, sp.Basic]:
        """
        Forward-mode automatic differentiation at a point using exact SymPy numbers.
        Returns (v, dr, dz, drr, dzz, drz) for expr at (rho0, z0).
        """
        rho, z = self.rho, self.z

        # Base cases
        if expr == rho:
            return rho0, sp.Integer(1), sp.Integer(0), sp.Integer(0), sp.Integer(0), sp.Integer(0)
        if expr == z:
            return z0, sp.Integer(0), sp.Integer(1), sp.Integer(0), sp.Integer(0), sp.Integer(0)
        if expr.is_Number:
            return expr, sp.Integer(0), sp.Integer(0), sp.Integer(0), sp.Integer(0), sp.Integer(0)

        # Helper for unary composition y=f(b)
        def compose_unary(bv, br, bz, brr, bzz, brz, f, fprime, fsecond):
            y = f(bv)
            fp = fprime(bv)
            fs = fsecond(bv)
            dr = fp * br
            dz = fp * bz
            drr = fs * br * br + fp * brr
            dzz = fs * bz * bz + fp * bzz
            drz = fs * br * bz + fp * brz
            return y, dr, dz, drr, dzz, drz

        # Operations
        func = expr.func
        args = expr.args

        # Addition and subtraction
        if func in (sp.Add,):
            vs = [self._ad_eval_point(arg, rho0, z0) for arg in args]
            v = sum(val[0] for val in vs)
            dr = sum(val[1] for val in vs)
            dz = sum(val[2] for val in vs)
            drr = sum(val[3] for val in vs)
            dzz = sum(val[4] for val in vs)
            drz = sum(val[5] for val in vs)
            return v, dr, dz, drr, dzz, drz

        # Multiplication
        if func in (sp.Mul,):
            # Reduce pairwise to keep formulas manageable
            def mul_pair(a, b):
                av, adr, adz, adrr, adzz, adrz = a
                bv, bdr, bdz, bdrr, bdzz, bdrz = b
                v = av * bv
                dr = adr * bv + av * bdr
                dz = adz * bv + av * bdz
                drr = adrr * bv + 2 * adr * bdr + av * bdrr
                dzz = adzz * bv + 2 * adz * bdz + av * bdzz
                drz = adrz * bv + adr * bdz + adz * bdr + av * bdrz
                return (v, dr, dz, drr, dzz, drz)
            vals = [self._ad_eval_point(arg, rho0, z0) for arg in args]
            acc = vals[0]
            for nxt in vals[1:]:
                acc = mul_pair(acc, nxt)
            return acc

        # Power
        if func in (sp.Pow,):
            base, expo = args
            bv, br, bz, brr, bzz, brz = self._ad_eval_point(base, rho0, z0)
            # Handle common exponents exactly (rational or integer)
            if expo.is_Number:
                k = sp.Rational(expo) if expo.is_Rational else sp.Integer(int(expo))
                f = lambda x: x**k
                fprime = lambda x: k * x**(k - 1)
                fsecond = lambda x: k * (k - 1) * x**(k - 2)
                return compose_unary(bv, br, bz, brr, bzz, brz, f, fprime, fsecond)
            # Fallback: y = exp(expo * log(base))
            # y = g(h), g=exp, h=expo*log(base)
            hv, hdr, hdz, hdrr, hdzz, hdrz = self._ad_eval_point(expo * sp.log(base), rho0, z0)
            f = lambda x: sp.exp(x)
            fprime = lambda x: sp.exp(x)
            fsecond = lambda x: sp.exp(x)
            return compose_unary(hv, hdr, hdz, hdrr, hdzz, hdrz, f, fprime, fsecond)

        # Unary functions: sqrt, exp, log
        if func in (sp.sqrt,):
            bv, br, bz, brr, bzz, brz = self._ad_eval_point(args[0], rho0, z0)
            f = lambda x: sp.sqrt(x)
            fprime = lambda x: sp.Rational(1, 2) * x**sp.Rational(-1, 2)
            fsecond = lambda x: -sp.Rational(1, 4) * x**sp.Rational(-3, 2)
            return compose_unary(bv, br, bz, brr, bzz, brz, f, fprime, fsecond)

        if func in (sp.exp,):
            bv, br, bz, brr, bzz, brz = self._ad_eval_point(args[0], rho0, z0)
            f = lambda x: sp.exp(x)
            fprime = lambda x: sp.exp(x)
            fsecond = lambda x: sp.exp(x)
            return compose_unary(bv, br, bz, brr, bzz, brz, f, fprime, fsecond)

        if func in (sp.log,):
            bv, br, bz, brr, bzz, brz = self._ad_eval_point(args[0], rho0, z0)
            f = lambda x: sp.log(x)
            fprime = lambda x: 1 / x
            fsecond = lambda x: -1 / (x**2)
            return compose_unary(bv, br, bz, brr, bzz, brz, f, fprime, fsecond)

        # Generic function not covered: fallback to sympy diff+subs (slower but rare)
        # Compute derivatives symbolically at the point for this node
        u_r = sp.diff(expr, rho).subs({rho: rho0, z: z0})
        u_z = sp.diff(expr, z).subs({rho: rho0, z: z0})
        u_rr = sp.diff(expr, rho, rho).subs({rho: rho0, z: z0})
        u_zz = sp.diff(expr, z, z).subs({rho: rho0, z: z0})
        u_rz = sp.diff(expr, rho, z).subs({rho: rho0, z: z0})
        u_v = expr.subs({rho: rho0, z: z0})
        return u_v, u_r, u_z, u_rr, u_zz, u_rz
        
    def _init_cache(self):
        """Initialize SQLite cache database."""
        self.conn = sqlite3.connect(self.cache_db)
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS validation_cache (
                expr_hash TEXT PRIMARY KEY,
                expr_str TEXT,
                is_valid INTEGER,
                constraint_value TEXT,
                reason TEXT,
                timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
            )
        """)
        self.conn.commit()
        
    def _get_expr_hash(self, expr: sp.Basic) -> str:
        """Get a hash for an expression for caching."""
        # Use string representation for hashing
        expr_str = str(expr)
        return hashlib.sha256(expr_str.encode()).hexdigest()
        
    def _check_cache(self, expr_hash: str) -> Optional[Tuple[bool, str]]:
        """Check if result is in cache."""
        cursor = self.conn.execute(
            "SELECT is_valid, reason FROM validation_cache WHERE expr_hash = ?",
            (expr_hash,)
        )
        result = cursor.fetchone()
        if result:
            return bool(result[0]), result[1]
        return None
        
    def _save_to_cache(self, expr_hash: str, expr_str: str, is_valid: bool, 
                       constraint_value: str, reason: str):
        """Save validation result to cache."""
        self.conn.execute("""
            INSERT OR REPLACE INTO validation_cache 
            (expr_hash, expr_str, is_valid, constraint_value, reason)
            VALUES (?, ?, ?, ?, ?)
        """, (expr_hash, expr_str, int(is_valid), constraint_value, reason))
        self.conn.commit()
        
    def _simplify_with_lean(self, expr: sp.Basic, timeout: float = 5.0) -> sp.Basic:
        """Use Lean normalizer to simplify expression if available."""
        if not self.use_lean or not self.lean_normalizer:
            # Use limited simplification for complex expressions
            if len(str(expr)) > 1000:
                return expand(expr)
            return simplify(expr)
            
        try:
            # Convert to string and normalize with Lean
            expr_str = str(expr)
            
            # Skip Lean for very complex expressions (prevent huge IPC)
            if len(expr_str) > 10000:
                logger.info("Expression too large for Lean call; returning expanded form")
                return expand(expr)
            
            start_time = time.time()
            batch_res = self.lean_normalizer.normalize_batch([(expr_str, 0)])
            if time.time() - start_time > timeout:
                logger.warning(f"Lean simplification timed out after {timeout}s")
                return expand(expr)
            
            if batch_res and isinstance(batch_res, list) and 'normalized' in batch_res[0]:
                normalized_str = batch_res[0]['normalized']
                # Fast-path: exact zero
                if normalized_str.strip() == '0':
                    return sp.Integer(0)
                return sp.sympify(normalized_str)
            
            # Fallback to a cheap expand if Lean returned nothing
            return expand(expr)
        except Exception as e:
            logger.warning(f"Lean simplification failed, using expand: {e}")
            return expand(expr)
        
    def validate(self, u: sp.Basic, check_regularity: bool = True, fast_point_only: bool = False) -> Tuple[bool, str]:
        """
        Validate if expression u satisfies the foliation constraint.
        
        Paper-faithful procedure (Section 2.4):
        - First evaluate the foliation constraint at (ρ,z) = (4/5, 6/7) using exact arithmetic.
        - If it is exactly 0 at that point, then test the constraint on the entire plane
          by asking Lean to simplify the symbolic determinant to 0.
        
        Args:
            u: Expression to validate
            check_regularity: If True, also check axis regularity
            
        Returns:
            (is_valid, reason) tuple
        """
        # Check cache first
        expr_hash = self._get_expr_hash(u)
        cached_result = self._check_cache(expr_hash)
        if cached_result is not None:
            return cached_result
        
        # Substitute symbols to ensure consistency
        u = u.subs([(s, self.rho if str(s) == 'rho' else self.z) 
                    for s in u.free_symbols if str(s) in ['rho', 'z']])
        
        try:
            # Check regularity on axis if requested
            if check_regularity:
                axis_value = u.subs(self.rho, 0)
                if axis_value.has(sp.oo, sp.zoo, sp.nan):
                    result = (False, "Singular on axis")
                    self._save_to_cache(expr_hash, str(u), False, "N/A", result[1])
                    return result
            
            # Compute derivatives quickly at the point using AD if fast mode
            rho_pt = sp.Rational(4, 5)
            z_pt = sp.Rational(6, 7)
            if fast_point_only:
                u_val, u_rho_pt, u_z_pt, u_rr_pt, u_zz_pt, u_rz_pt = self._ad_eval_point(u, rho_pt, z_pt)
                u_rho = sp.Function('u_rho')
                u_z = sp.Function('u_z')
                u_rho_rho = sp.Function('u_rr')
                u_z_z = sp.Function('u_zz')
            else:
                u_rho = u.diff(self.rho)
                u_z = u.diff(self.z)
            
            # Check if gradient is zero (trivial case)
            if u_rho == 0 and u_z == 0:
                result = (False, "Zero gradient (constant expression)")
                self._save_to_cache(expr_hash, str(u), False, "N/A", result[1])
                return result
            
            # Second derivatives
            if fast_point_only:
                u_rho_rho = sp.Function('u_rr')
                u_z_z = sp.Function('u_zz')
            else:
                u_rho_rho = u_rho.diff(self.rho)
                u_z_z = u_z.diff(self.z)
            
            # A and B (Eq. 2.10), with optional rotation
            A_non_rotating = u_rho_rho + u_z_z - u_rho/self.rho
            B_non_rotating = u_rho**2 + u_z**2

            if self.Omega != 0:
                A = (1 - self.rho**2 * self.Omega**2) * (u_rho_rho + u_z_z) - \
                    (1 + self.rho**2 * self.Omega**2)/self.rho * u_rho
                B = (1 - self.rho**2 * self.Omega**2) * (u_rho**2 + u_z**2)
            else:
                A = A_non_rotating
                B = B_non_rotating
            
            # Lie derivative along T = u_z ∂_ρ - u_ρ ∂_z
            def lie_derivative_T(f):
                if fast_point_only:
                    # In fast mode we only need determinant at the point; we will substitute numeric derivatives directly
                    return u_z * f.diff(self.rho) - u_rho * f.diff(self.z)
                return u_z * f.diff(self.rho) - u_rho * f.diff(self.z)
            
            LT_A = lie_derivative_T(A)
            LT_B = lie_derivative_T(B)
            L2T_A = lie_derivative_T(LT_A)
            L2T_B = lie_derivative_T(LT_B)
            
            # Determinant of [[LT_A, LT_B], [L2T_A, L2T_B]]
            det_M = det(Matrix([[LT_A, LT_B], [L2T_A, L2T_B]]))
            
            # Step 1 (paper): check at (ρ, z) = (4/5, 6/7)
            test_point = {self.rho: rho_pt, self.z: z_pt}
            if fast_point_only:
                # Substitute first- and second-order derivatives numerically at the point to avoid symbolic growth
                subs_map = {
                    self.rho: rho_pt,
                    self.z: z_pt,
                    u_rho: u_rho_pt,
                    u_z: u_z_pt,
                    u_rho_rho: u_rr_pt,
                    u_z_z: u_zz_pt,
                }
                det_at_point = det_M.subs(subs_map)
            else:
                det_at_point = det_M.subs(test_point)
            det_at_point = det_M.subs(test_point)
            
            # Fast path: avoid Lean, prefer cheap simplifications
            det_at_point_simple = det_at_point
            try:
                # Try very cheap canonical reductions first
                det_at_point_simple = sp.cancel(sp.together(det_at_point_simple))
                if det_at_point_simple.is_Number:
                    if det_at_point_simple == 0:
                        if fast_point_only:
                            result = (True, "Valid foliation (point check = 0)")
                            self._save_to_cache(expr_hash, str(u), True, "point_only", result[1])
                            return result
                    else:
                        result = (False, "Invalid (point check != 0)")
                        self._save_to_cache(expr_hash, str(u), False, "point_only", result[1])
                        return result
                # Try a small simplify budget
                det_at_point_simple = sp.simplify(det_at_point_simple)
            except Exception:
                pass
            
            # If still not a number, fall back to high-precision numeric check at the point
            try:
                det_val = complex(det_at_point_simple.evalf(50))
                if abs(det_val) < 1e-20:
                    if fast_point_only:
                        result = (True, "Valid foliation (point check ≈ 0)")
                        self._save_to_cache(expr_hash, str(u), True, "point_only", result[1])
                        return result
                else:
                    result = (False, f"Invalid (point check ≈ {abs(det_val):.2e})")
                    self._save_to_cache(expr_hash, str(u), False, "point_only", result[1])
                    return result
            except Exception:
                # If evaluation fails, mark as indeterminate
                result = (False, "Could not evaluate point check")
                self._save_to_cache(expr_hash, str(u), False, "point_only", result[1])
                return result
            
            # Optional full-plane check (slow). Only run when not in fast mode.
            if not fast_point_only:
                # Prefer Lean only for small expressions
                det_str = str(det_M)
                if self.use_lean and self.lean_normalizer and len(det_str) < 3000:
                    det_symbolic_simpl = self._simplify_with_lean(det_M)
                    if det_symbolic_simpl == 0:
                        result = (True, "Valid foliation (Lean: det = 0 symbolically)")
                        self._save_to_cache(expr_hash, str(u), True, "lean_symbolic", result[1])
                        return result
                    result = (False, "Invalid (Lean could not simplify det to 0 symbolically)")
                    self._save_to_cache(expr_hash, str(u), False, "lean_symbolic", result[1])
                    return result
                else:
                    # Last resort: expanded check
                    try:
                        if expand(det_M) == 0:
                            result = (True, "Valid foliation (expanded det = 0)")
                        else:
                            result = (False, "Invalid (expanded det != 0)")
                    except Exception:
                        result = (False, "Could not simplify det symbolically")
                    self._save_to_cache(expr_hash, str(u), result[0], "symbolic", result[1])
                    return result
            
            # If we reached here, treat as inconclusive
            result = (False, "Inconclusive")
            self._save_to_cache(expr_hash, str(u), False, "inconclusive", result[1])
            return result
        
        except Exception as e:
            result = (False, f"Error: {str(e)}")
            self._save_to_cache(expr_hash, str(u), False, "Error", result[1])
            return result
            
    def validate_known_solutions(self) -> Dict[str, bool]:
        """Test the validator on the 7 known solutions from the paper."""
        known_solutions = {
            'Vertical': self.rho**2,
            'X-point': self.rho**2 * self.z,
            'Radial': 1 - self.z/sqrt(self.rho**2 + self.z**2),
            'Dipolar': self.rho**2/(self.rho**2 + self.z**2)**(sp.Rational(3,2)),
            'Parabolic': sqrt(self.rho**2 + self.z**2) - self.z,
            'Hyperbolic': sqrt(self.z**2 + (self.rho - 1)**2) - sqrt(self.z**2 + (self.rho + 1)**2),
            'Bent': self.rho**2 * exp(-2*self.z)
        }
        
        results = {}
        print("Validating known solutions from the paper:")
        print("=" * 60)
        
        for name, expr in known_solutions.items():
            is_valid, reason = self.validate(expr)
            results[name] = is_valid
            status = "✓" if is_valid else "✗"
            print(f"{status} {name:12s}: {reason}")
            
        return results
        
    def get_cache_stats(self) -> Dict[str, int]:
        """Get statistics about the cache."""
        cursor = self.conn.execute("""
            SELECT 
                COUNT(*) as total,
                SUM(is_valid) as valid,
                COUNT(*) - SUM(is_valid) as invalid
            FROM validation_cache
        """)
        row = cursor.fetchone()
        return {
            'total': row[0],
            'valid': row[1],
            'invalid': row[2]
        }
        
    def clear_cache(self):
        """Clear the validation cache."""
        self.conn.execute("DELETE FROM validation_cache")
        self.conn.commit()
        
    def __del__(self):
        """Close database connection."""
        if hasattr(self, 'conn'):
            self.conn.close()


def test_validator():
    """Test the precise validator."""
    
    # Test non-rotating solutions
    print("\nTesting non-rotating solutions (Omega = 0):")
    validator_no_rotation = PreciseFoliationValidator(Omega=0)
    results_no_rotation = validator_no_rotation.validate_known_solutions()

    # Test rotating solutions (with constant Omega)
    print("\nTesting rigidly rotating solutions (Omega = 1):")
    Omega_test = 1 
    validator_rotation = PreciseFoliationValidator(Omega=Omega_test)
    
    # Re-using non-rotating solutions as they should still be valid with rotation
    results_rotation = validator_rotation.validate_known_solutions()

    # Test some expressions that should fail
    print("\n" + "="*60)
    print("Testing expressions that should fail (with and without rotation):")
    
    test_fail = [
        (sp.symbols('rho'), 'Linear rho'),
        (sp.symbols('z'), 'Linear z'),
        (sp.symbols('rho')*sp.symbols('z'), 'Product rho*z'),
        (sp.symbols('rho')**3, 'Cubic rho'),
        (sp.exp(sp.symbols('z')), 'Exponential z')
    ]
    
    for expr, desc in test_fail:
        is_valid_no_rot, reason_no_rot = validator_no_rotation.validate(expr)
        status_no_rot = "✗" if not is_valid_no_rot else "✓"
        print(f"{status_no_rot} (No-rot) {desc:20s}: {reason_no_rot}")

        is_valid_rot, reason_rot = validator_rotation.validate(expr)
        status_rot = "✗" if not is_valid_rot else "✓"
        print(f"{status_rot} (Rot)    {desc:20s}: {reason_rot}")

    # Show cache stats
    print("\n" + "="*60)
    stats = validator_no_rotation.get_cache_stats()
    print(f"Cache statistics: {stats}")
    
    return {
        "non_rotating": results_no_rotation,
        "rotating": results_rotation
    }


if __name__ == "__main__":
    test_validator()
