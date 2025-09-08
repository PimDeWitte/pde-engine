"""
Bridge between Python and Lean for expression normalization.
Implements the mathematical functions from the Force-Free Foliations paper.
"""

import subprocess
import json
import os
import tempfile
from typing import List, Dict, Set, Tuple, Any, Optional
import sympy as sp
from pathlib import Path
import hashlib


class LeanNormalizer:
    """Normalizes expressions using Lean for canonical forms."""
    
    def __init__(self):
        self.lean_path = Path(__file__).parent
        self.lake_exe = self._find_lake()
        self._build_lean_project()
        
    def _find_lake(self) -> str:
        """Find the lake executable."""
        # Try common locations
        locations = [
            "lake",  # In PATH
            "~/.elan/bin/lake",
            "/usr/local/bin/lake",
        ]
        
        for loc in locations:
            try:
                expanded = os.path.expanduser(loc)
                result = subprocess.run([expanded, "--version"], 
                                      capture_output=True, text=True)
                if result.returncode == 0:
                    return expanded
            except:
                continue
                
        # If not found, assume it's in PATH
        return "lake"
    
    def _build_lean_project(self):
        """Build the Lean project if needed."""
        try:
            # Check if already built
            olean_file = self.lean_path / ".lake" / "build" / "lib" / "PhysicsExpr.olean"
            if olean_file.exists():
                return
                
            # Build the project
            result = subprocess.run(
                [self.lake_exe, "build"],
                cwd=self.lean_path,
                capture_output=True,
                text=True
            )
            
            if result.returncode != 0:
                print(f"Warning: Failed to build Lean project: {result.stderr}")
        except Exception as e:
            print(f"Warning: Could not build Lean project: {e}")
    
    def normalize(self, expr_str: str) -> str:
        """
        Normalize a single expression to canonical form.
        For now, uses SymPy as fallback since full Lean integration is complex.
        """
        try:
            expr = sp.sympify(expr_str)
            # Apply canonical ordering and simplification
            normalized = self._canonical_form(expr)
            return str(normalized)
        except:
            return expr_str
    
    def _canonical_form(self, expr):
        """Convert expression to canonical form."""
        # Expand and collect terms
        expr = sp.expand(expr)
        
        # Collect like terms
        if expr.has(sp.Symbol('rho')) and expr.has(sp.Symbol('z')):
            expr = sp.collect(expr, [sp.Symbol('rho'), sp.Symbol('z')])
        
        # Apply specific simplifications
        expr = self._apply_rules(expr)
        
        return expr
    
    def _apply_rules(self, expr):
        """Apply simplification rules from the paper."""
        # Handle special cases
        rho = sp.Symbol('rho', positive=True)
        z = sp.Symbol('z')
        
        # Common substitutions
        substitutions = [
            (sp.exp(sp.log(rho)), rho),
            (sp.log(sp.exp(z)), z),
            (sp.sqrt(rho**2), rho),  # Since rho is positive
            (rho/rho, 1),
            (z - z, 0),
        ]
        
        for pattern, replacement in substitutions:
            expr = expr.subs(pattern, replacement)
        
        return expr
    
    def compute_signature(self, expr_str: str) -> int:
        """Compute a hash signature for the expression."""
        normalized = self.normalize(expr_str)
        return hash(normalized) % (2**32)
    
    def batch_normalize(self, expressions: List[str]) -> List[str]:
        """Normalize a batch of expressions."""
        return [self.normalize(expr) for expr in expressions]


class FastExpressionGenerator:
    """
    Fast expression generator using Lean-inspired canonical forms.
    Implements the mathematical operations from Section 2.4 of the paper.
    """
    
    def __init__(self, normalizer: LeanNormalizer):
        self.normalizer = normalizer
        self.seen_signatures = set()
        
    def generate_expressions(self, primitives: List, unary_ops: Dict, 
                           binary_ops: Dict, max_depth: int) -> Dict[int, List[str]]:
        """
        Generate expressions up to max_depth using the paper's method.
        
        This implements:
        - The 4 primitives from Section 2.4
        - Binary operations including special ones from the paper
        - Unary operations for transformations
        """
        import sys
        expressions_by_depth = {1: []}
        all_seen = set()
        
        # Add primitives at depth 1
        for p in primitives:
            p_str = str(p)
            normalized = self.normalizer.normalize(p_str)
            sig = self.normalizer.compute_signature(normalized)
            
            if sig not in self.seen_signatures:
                expressions_by_depth[1].append(normalized)
                self.seen_signatures.add(sig)
                all_seen.add(normalized)
        
        print(f"Depth 1: {len(expressions_by_depth[1])} expressions")
        
        # Generate higher depths
        for depth in range(2, max_depth + 1):
            expressions_by_depth[depth] = []
            new_count = 0
            print(f"\nGenerating depth {depth}...")
            sys.stdout.flush()
            
            # Apply unary operations
            unary_count = 0
            total_unary = len(expressions_by_depth[depth - 1]) * len(unary_ops)
            print(f"  Processing {total_unary} unary combinations...")
            sys.stdout.flush()
            
            for expr_str in expressions_by_depth[depth - 1]:
                try:
                    expr = sp.sympify(expr_str)
                    for op_name, op_func in unary_ops.items():
                        unary_count += 1
                        if unary_count % 100 == 0:
                            print(f"    Unary: {unary_count}/{total_unary}")
                            sys.stdout.flush()
                        try:
                            result = op_func(expr)
                            if result is None:
                                continue
                            
                            # Normalize immediately
                            result_str = str(result)
                            normalized = self.normalizer.normalize(result_str)
                            
                            if normalized not in all_seen:
                                expressions_by_depth[depth].append(normalized)
                                all_seen.add(normalized)
                                new_count += 1
                        except:
                            continue
                except:
                    continue
            
            # Apply binary operations
            binary_count = 0
            print(f"  Processing binary combinations...")
            sys.stdout.flush()
            
            for d1 in range(1, depth):
                d2 = depth - d1
                if d2 < 1 or d2 >= depth:
                    continue
                
                for expr1_str in expressions_by_depth[d1]:
                    for expr2_str in expressions_by_depth[d2]:
                        try:
                            expr1 = sp.sympify(expr1_str)
                            expr2 = sp.sympify(expr2_str)
                            
                            for op_name, op_func in binary_ops.items():
                                binary_count += 1
                                if binary_count % 1000 == 0:
                                    print(f"    Binary: {binary_count} processed, {new_count} new found")
                                    sys.stdout.flush()
                                # Handle special operations from the paper
                                if op_name in ['sqrt_shift_neg', 'sqrt_shift_pos']:
                                    # These are specific to building hyperbolic solution
                                    # Only apply to coordinate-like expressions
                                    if not self._is_coordinate_expr(expr1_str):
                                        continue
                                
                                try:
                                    result = op_func(expr1, expr2)
                                    if result is None:
                                        continue
                                    
                                    # Normalize
                                    result_str = str(result)
                                    normalized = self.normalizer.normalize(result_str)
                                    
                                    if normalized not in all_seen:
                                        expressions_by_depth[depth].append(normalized)
                                        all_seen.add(normalized)
                                        new_count += 1
                                        
                                    # For non-commutative ops, try reversed order
                                    if op_name in ['sub', 'div', 'geom_sum']:
                                        result = op_func(expr2, expr1)
                                        if result is not None:
                                            result_str = str(result)
                                            normalized = self.normalizer.normalize(result_str)
                                            
                                            if normalized not in all_seen:
                                                expressions_by_depth[depth].append(normalized)
                                                all_seen.add(normalized)
                                                new_count += 1
                                except:
                                    continue
                        except:
                            continue
            
            print(f"Depth {depth}: {new_count} new expressions")
        
        return expressions_by_depth
    
    def _is_coordinate_expr(self, expr_str: str) -> bool:
        """Check if expression is coordinate-like (for special operations)."""
        # Simple coordinate expressions used in the paper
        return expr_str in ['rho', 'z', 'rho + z', 'rho - z', 'z - rho']


class FoliationConstraintValidator:
    """
    Validates the foliation constraint from equation 2.14 of the paper.
    
    The constraint is: det([[LT A, LT B], [LT^2 A, LT^2 B]]) = 0
    where T = u_z ∂_ρ - u_ρ ∂_z is the tangent vector field.
    """
    
    def __init__(self):
        self.rho = sp.Symbol('rho', real=True, positive=True)
        self.z = sp.Symbol('z', real=True)
    
    def validate(self, u_expr: sp.Basic) -> Tuple[bool, str]:
        """
        Check if expression satisfies the foliation constraint.
        
        Returns (is_valid, reason_string)
        """
        try:
            # Compute derivatives
            u_rho = sp.diff(u_expr, self.rho)
            u_z = sp.diff(u_expr, self.z)
            u_rho_rho = sp.diff(u_rho, self.rho)
            u_z_z = sp.diff(u_z, self.z)
            u_rho_z = sp.diff(u_rho, self.z)
            
            # Compute A and B from equation 2.10
            A = u_rho_rho + u_z_z - u_rho / self.rho
            B = u_rho**2 + u_z**2
            
            # Tangent operator T = u_z ∂_ρ - u_ρ ∂_z
            # LT A means applying T to A
            LT_A = u_z * sp.diff(A, self.rho) - u_rho * sp.diff(A, self.z)
            LT_B = u_z * sp.diff(B, self.rho) - u_rho * sp.diff(B, self.z)
            
            # Second application of T
            LT2_A = u_z * sp.diff(LT_A, self.rho) - u_rho * sp.diff(LT_A, self.z)
            LT2_B = u_z * sp.diff(LT_B, self.rho) - u_rho * sp.diff(LT_B, self.z)
            
            # The determinant constraint
            det = LT_A * LT2_B - LT_B * LT2_A
            
            # Simplify the determinant
            det_simplified = sp.simplify(det)
            
            # Check if it's zero
            is_zero = det_simplified == 0
            
            if is_zero:
                return True, "Satisfies foliation constraint"
            else:
                # Check numerically at a test point
                test_val = det_simplified.subs({self.rho: sp.Rational(4, 5), 
                                               self.z: sp.Rational(6, 7)})
                if abs(float(test_val)) < 1e-10:
                    return True, "Satisfies constraint (numerically)"
                else:
                    return False, f"Constraint not satisfied: det = {det_simplified}"
                    
        except Exception as e:
            return False, f"Error checking constraint: {e}"
    
    def check_regularity(self, u_expr: sp.Basic) -> bool:
        """
        Check if expression is regular on the axis (rho = 0).
        The expression must be constant along the axis.
        """
        try:
            # Substitute rho = 0
            on_axis = u_expr.subs(self.rho, 0)
            
            # Check if it depends on z
            if on_axis.has(self.z):
                # If it still has z, it's not constant on axis
                return False
            
            return True
        except:
            # If substitution fails, assume not regular
            return False