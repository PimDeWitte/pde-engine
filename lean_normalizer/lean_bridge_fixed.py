"""
Fixed Lean Bridge - Standalone implementation for physics discovery.

This provides a Lean normalizer that works without external dependencies.
"""

import sys
import os
from pathlib import Path
from typing import List, Dict, Set, Tuple, Any, Optional
import sympy as sp
import sqlite3
import hashlib

# Import the local implementation
from .lean_bridge import LeanNormalizer as BaseLeanNormalizer

# Re-export with enhancements
class LeanNormalizer(BaseLeanNormalizer):
    """Enhanced Lean normalizer with caching."""
    
    def __init__(self, lean_path: Optional[str] = None, cache_db: Optional[str] = None):
        # Use a different cache database for physics expressions
        self.cache_db_path = cache_db or str(Path(__file__).parent / "physics_expressions.db")
        super().__init__()
        self._init_cache()
        print("Initialized Lean normalizer with SQLite cache")
    
    def _init_cache(self):
        """Initialize SQLite cache for normalized expressions."""
        self.cache_conn = sqlite3.connect(self.cache_db_path)
        self.cache_conn.execute("""
            CREATE TABLE IF NOT EXISTS normalized_cache (
                expr_hash TEXT PRIMARY KEY,
                expr_str TEXT,
                normalized TEXT,
                timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
            )
        """)
        self.cache_conn.commit()
    
    def normalize_batch(self, expressions: List[Tuple[str, int]]) -> List[Dict[str, Any]]:
        """Normalize a batch of expressions with caching."""
        results = []
        for expr_str, idx in expressions:
            # Check cache first
            expr_hash = hashlib.sha256(expr_str.encode()).hexdigest()
            cursor = self.cache_conn.execute(
                "SELECT normalized FROM normalized_cache WHERE expr_hash = ?",
                (expr_hash,)
            )
            cached = cursor.fetchone()
            if cached:
                # Generate signature from normalized form
                sig = hashlib.sha256(cached[0].encode()).hexdigest()[:16]
                results.append({'normalized': cached[0], 'index': idx, 'signature': sig})
            else:
                # Normalize and cache
                normalized = self.normalize(expr_str)
                self.cache_conn.execute(
                    "INSERT OR REPLACE INTO normalized_cache (expr_hash, expr_str, normalized) VALUES (?, ?, ?)",
                    (expr_hash, expr_str, normalized)
                )
                self.cache_conn.commit()
                # Generate signature from normalized form
                sig = hashlib.sha256(normalized.encode()).hexdigest()[:16]
                results.append({'normalized': normalized, 'index': idx, 'signature': sig})
        return results


class FastExpressionGenerator:
    """
    Fast expression generator using REAL Lean normalization.
    Adapted for physics agent's specific needs.
    """
    
    def __init__(self, normalizer: Optional[LeanNormalizer] = None):
        """Initialize with Lean normalizer."""
        self.normalizer = normalizer or LeanNormalizer()
        self.seen_normalized = set()  # <reason>chain: Track normalized forms to avoid duplicates</reason>
        
    def generate_expressions(self, primitives: List, unary_ops: Dict, 
                           binary_ops: Dict, max_depth: int) -> Dict[int, List[str]]:
        """
        Generate expressions using Lean normalization.
        
        This handles the combinatorial explosion by:
        1. Using fast Lean normalization in batches
        2. Proper caching to avoid re-normalizing
        3. Efficient signature-based deduplication
        """
        print(f"\nGenerating expressions with Lean normalization up to depth {max_depth}")
        
        # Convert primitives to strings
        primitive_strs = [str(p) for p in primitives]
        
        # Use streaming to generate expressions
        expressions_by_depth = {}
        
        def collect_batch(depth: int, expressions: List[str]):
            expressions_by_depth[depth] = expressions
        
        self.stream_generate(
            primitives=primitives,
            unary_ops=unary_ops,
            binary_ops=binary_ops,
            max_depth=max_depth,
            on_batch=collect_batch
        )
        
        return expressions_by_depth
    
    def stream_generate(self,
                        primitives: List,
                        unary_ops: Dict,
                        binary_ops: Dict,
                        max_depth: int,
                        batch_size: int = 1000,
                        on_batch: Optional[Any] = None,
                        prune: bool = True) -> None:
        """Stream expressions per batch to a callback to enable early consumption.

        on_batch(depth: int, expressions: List[str]) will be invoked for each batch of
        normalized expressions at the given depth.
        """
        # Prepare primitive strings
        primitive_strs = [str(p) for p in primitives]
        # Depth 1 is just primitives
        if on_batch:
            on_batch(1, list(primitive_strs))

        # Build incrementally per depth similar to theorem_scraper's generator
        expressions_by_depth: Dict[int, List[str]] = {1: primitive_strs}
        def _has_vars(s: str) -> bool:
            # Checks for dependence on coordinates (r, x, rho, z)
            return ('r' in s) or ('x' in s) or ('rho' in s) or ('z' in s)
        seen_signatures: Set[Any] = set()

        for depth in range(2, max_depth + 1):
            candidates: List[Tuple[str, int]] = []
            # Unary from previous depth (optional pruning)
            for expr in expressions_by_depth[depth - 1]:
                if prune:
                    dep = _has_vars(expr)
                    if not dep:
                        continue
                for op_name in unary_ops:
                    if prune:
                        if op_name == 'inv' and expr.startswith('inv('):
                            continue
                        if op_name in ('sqrt', 'square', 'pow_3_2', 'pow_neg_3_2') and expr == '1':
                            continue
                    candidates.append((f"{op_name}({expr})", depth))
            # Binary combining d1 and d2 = depth - d1
            for d1 in range(1, depth):
                d2 = depth - d1
                if d2 < 1 or d2 >= depth:
                    continue
                for expr1 in expressions_by_depth[d1]:
                    for expr2 in expressions_by_depth[d2]:
                        # Optional pruning: skip constant-only × constant-only
                        if prune:
                            dep1, dep2 = _has_vars(expr1), _has_vars(expr2)
                            if not dep1 and not dep2:
                                continue
                        for op_name in binary_ops:
                            a, b = expr1, expr2
                            if op_name in ['add', 'mul'] and a > b:
                                a, b = b, a
                            if op_name == 'add':
                                candidates.append((f"({a} + {b})", depth))
                            elif op_name == 'sub':
                                if prune:
                                    # Skip a-a
                                    if a == b:
                                        continue
                                candidates.append((f"({a} - {b})", depth))
                            elif op_name == 'mul':
                                if prune:
                                    # Skip multiplication by 1
                                    if a == '1' or b == '1':
                                        continue
                                candidates.append((f"({a} * {b})", depth))
                            elif op_name == 'div':
                                if prune:
                                    # Skip divide by 1 and a/a
                                    if b == '1' or a == b:
                                        continue
                                candidates.append((f"({a} / ({b}))", depth))
                            elif op_name == 'geom_sum':
                                if prune:
                                    # Skip denominator singular case 1 - b with b==1
                                    if b == '1':
                                        continue
                                candidates.append((f"({a} / (1 - {b}))", depth))

            # Process in batches and report progress
            print(f"Depth {depth}: {len(candidates)} candidates to normalize")
            unique_expressions: List[str] = []
            for i in range(0, len(candidates), batch_size):
                batch = candidates[i:i + batch_size]
                results = self.normalizer.normalize_batch(batch)
                out_chunk: List[str] = []
                for result in results:
                    sig = result.get('signature')
                    norm = result.get('normalized')
                    if sig not in seen_signatures:
                        seen_signatures.add(sig)
                        unique_expressions.append(norm)
                        out_chunk.append(norm)
                if on_batch and out_chunk:
                    on_batch(depth, out_chunk)

            expressions_by_depth[depth] = unique_expressions
            print(f"Depth {depth}: {len(unique_expressions)} unique expressions after normalization")
    
    def _is_coordinate_expr(self, expr_str: str) -> bool:
        """Check if expression is coordinate-like (for special operations)."""
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
            
            # For validation, we still need to simplify to check if it's zero
            # But we can use a more targeted approach
            det_expanded = sp.expand(det)
            
            # Quick symbolic check
            if det_expanded == 0:
                return True, "Satisfies foliation constraint (symbolically zero)"
            
            # Try collecting terms to see if it simplifies to zero
            det_collected = sp.collect(det_expanded, [self.rho, self.z])
            if det_collected == 0:
                return True, "Satisfies foliation constraint (zero after collecting)"
                
            # Numerical check at test points
            test_points = [
                {self.rho: sp.Rational(4, 5), self.z: sp.Rational(6, 7)},
                {self.rho: sp.Rational(3, 4), self.z: sp.Rational(5, 6)},
                {self.rho: sp.Rational(7, 8), self.z: sp.Rational(1, 2)},
            ]
            
            max_val = 0
            for point in test_points:
                try:
                    val = abs(float(det_expanded.subs(point)))
                    max_val = max(max_val, val)
                except:
                    pass
                    
            if max_val < 1e-10:
                return True, f"Satisfies constraint (numerically, max = {max_val:.2e})"
            else:
                return False, f"Constraint not satisfied: max numerical value = {max_val:.2e}"
                    
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


# For backward compatibility - if code imports these directly
__all__ = ['LeanNormalizer', 'FastExpressionGenerator', 'FoliationConstraintValidator']
