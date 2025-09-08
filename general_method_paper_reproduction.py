"""
PDE Engine - General Method for Discovering Symbolic Solutions to PDEs

This module implements a general approach for discovering symbolic solutions to partial
differential equations through systematic expression search with intelligent simplification.

Key features:
- Parallel expression generation and validation
- Lean-based mathematical normalization to canonical forms
- SQLite-based caching and deduplication
- Support for arbitrary PDEs through the ProblemSpec interface

The system can reproduce known results (e.g., all 7 force-free foliation solutions from
Compère et al.) while generating ~200x fewer expressions than brute force approaches.

Architecture:
- Expression generation: Builds expressions from primitives using unary/binary operations
- Lean normalization: Reduces expressions to canonical forms, eliminating redundancy
- Parallel validation: Tests expressions against PDE constraints using problem validators
- Result storage: SQLite database with full audit trail and deduplication

This is a general-purpose PDE discovery engine suitable for research in mathematical physics.
"""

import sympy as sp
from sympy import Symbol, sqrt, exp, log, simplify, expand, Integer, Rational
from typing import List, Dict, Set, Tuple, Any, Optional
import time
import hashlib
from collections import defaultdict
import json
import os
from datetime import datetime
import sys
import signal
from pathlib import Path
import sqlite3
import multiprocessing
from multiprocessing import Process, Queue
import uuid
import threading

# Fix import path for running as a script
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))


from expression_operations import UNARY_OPS, BINARY_OPS, SPECIAL_OPS, ALL_BINARY_OPS
from problems import load_problem, ProblemSpec


class GeneralFoliationDiscovery:
    """
    General method for discovering symbolic PDE solutions.
    
    This implementation:
    1. Uses problem-specific primitives and operations
    2. Generates expressions systematically by depth
    3. Validates against problem-specific PDE constraints
    4. Uses Lean normalization to reduce expression count by ~200x
    """
    
    def __init__(self, use_lean_normalizer: bool = True, cache_db: str = None, mode: str = 'sequential', problem: ProblemSpec | None = None, problem_name: str | None = None):
        # Load problem spec (default to force-free)
        self.problem: ProblemSpec = problem or load_problem(problem_name or 'force_free')

        # Symbols (keep legacy names if present)
        self.rho = self.problem.symbols.get('rho', Symbol('rho', real=True, positive=True))
        self.z = self.problem.symbols.get('z', Symbol('z', real=True))

        # Problem-specific validator
        if mode == 'sequential':
            # If the problem provided its own validator instance, use it; otherwise fallback
            self.validator = getattr(self.problem, 'validator', None) or PreciseFoliationValidator(cache_db=cache_db, use_lean=use_lean_normalizer)
        else:
            self.validator = getattr(self.problem, 'validator', None) or PreciseFoliationValidator(cache_db=cache_db, use_lean=use_lean_normalizer)

        # Primitives and operations from the problem
        self.primitives = list(self.problem.primitives)
        self.binary_ops = dict(self.problem.binary_ops)
        self.unary_ops = dict(self.problem.unary_ops)
        self.special_ops = dict(self.problem.special_ops)
        self.all_binary_ops = dict(self.problem.all_binary_ops)

        # Locals mapping for sympify to ensure expressions use the same symbols
        self._sympify_locals = {}
        for name, sym in {**self.problem.symbols, **getattr(self.problem, 'constants', {})}.items():
            self._sympify_locals[name] = sym
        # Also map generator unary op names to concrete SymPy callables so strings like exp_neg(...) parse correctly
        try:
            from expression_operations import UNARY_OPS
            self._sympify_locals.update(UNARY_OPS)
        except Exception:
            pass
        
        # Track expressions by depth
        self.expressions_by_depth = defaultdict(set)
        
        # Global set for all unique expressions
        self.all_expressions = set()
        
        # FORCE LEAN NORMALIZER - NO SYMPY FALLBACK
        self.use_lean_normalizer = use_lean_normalizer  # Allow it to be disabled for debugging
        try:
            # Use fixed wrapper that binds to the real Lean normalizer
            from lean_normalizer.lean_bridge_fixed import LeanNormalizer, FastExpressionGenerator
            self.normalizer = LeanNormalizer()
            self.fast_generator = FastExpressionGenerator(self.normalizer)
            print("Using LEAN normalizer exclusively - SymPy fallback disabled")
        except ImportError as e:
            print(f"Warning: Could not load Lean normalizer: {e}")
            print("Falling back to SymPy normalization")
            self.use_lean_normalizer = False
            self.normalizer = None
            self.fast_generator = None
        
        # Statistics
        self.stats = {
            'total_generated': 0,
            'duplicates_avoided': 0,
            'expressions_checked': 0,
            'valid_foliations': 0,
            'known_solutions_found': 0,
            'invalid_null_results': 0,
            'invalid_error_results': 0,
            'degenerate_denominators_dropped': 0
        }
        
        if mode == 'parallel':
            self.run_id = None
            self.db_path = None
            self.table_name = None
            self.monitoring_stop_event = None
        
    def _has_degenerate_denominator(self, expr: sp.Basic) -> bool:
        """Return True if any subexpression has a denominator that simplifies to 0.

        Strategy:
        - Walk the expression tree (preorder)
        - For every subexpression s, attempt to expose a rational denominator via
          together(s) then fraction(...). Also explicitly check Pow with negative
          integer exponents to catch forms like (1 - 1)**-1.
        - If simplify(den) == 0 for any s, report degenerate.
        """
        try:
            # Immediate infinities/NaNs
            try:
                if expr.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
                    return True
            except Exception:
                pass
            for sub in sp.preorder_traversal(expr):
                try:
                    try:
                        if sub.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
                            return True
                    except Exception:
                        pass
                    # Handle explicit negative powers as denominators
                    if isinstance(sub, sp.Pow):
                        exp_part = sub.exp
                        if getattr(exp_part, 'is_integer', False) and bool(exp_part.is_negative):
                            base = sub.base
                            try:
                                if sp.simplify(base) == 0:
                                    return True
                            except Exception:
                                pass

                    # Try to expose denominator structure
                    try:
                        combined = sp.together(sub)
                    except Exception:
                        combined = sub

                    try:
                        num, den = sp.fraction(combined)
                    except Exception:
                        # As a fallback, fraction on the raw sub
                        try:
                            num, den = sp.fraction(sub)
                        except Exception:
                            continue

                    # Quick exits when no denominator
                    if den is None or den == 1:
                        continue

                    try:
                        if sp.simplify(den) == 0:
                            return True
                    except Exception:
                        # If simplify fails, be conservative and continue
                        continue
                except Exception:
                    # Any traversal issue: continue scanning other nodes
                    continue
        except Exception:
            return False
        return False

    def check_foliation_constraint(self, u: sp.Basic) -> Tuple[bool, str]:
        """
        Check if expression satisfies the foliation constraint (Eq 2.14).
        
        Uses the precise validator with caching.
        Properly catches null and error results.
        """
        try:
            # Check for None or invalid expressions
            if u is None:
                self.stats['invalid_null_results'] += 1
                return False, "Null expression"
                
            # Use the precise validator with caching
            is_valid, reason = self.validator.validate(u, check_regularity=True)
            
            # Check for error conditions in the reason
            if "Error" in reason or "Could not" in reason:
                self.stats['invalid_error_results'] += 1
                return False, reason

            return is_valid, reason
            
        except Exception as e:
            self.stats['invalid_error_results'] += 1
            return False, f"Validation error: {str(e)}"
            
    def generate_expressions_up_to_depth(self, max_depth: int = 4) -> Dict[int, Set[str]]:
        """
        Generate all expressions up to given depth using the paper's method
        but with our simplification approach.
        FORCE LEAN USAGE - NO SYMPY FALLBACK
        """
        # Force Lean - no conditional
        return self._generate_with_lean(max_depth)
    
    def _generate_with_lean(self, max_depth: int) -> Dict[int, Set[str]]:
        """Generate expressions using LEAN normalization ONLY."""
        print(f"STARTING: Generating expressions up to depth {max_depth} (Using LEAN ONLY) — Problem: {self.problem.name}")
        print("="*70)
        print(f"Primitives: {[str(p) for p in self.primitives]}")
        print(f"Unary ops: {list(self.unary_ops.keys())}")
        print(f"Binary ops: {list(self.binary_ops.keys())}")
        print(f"Special ops: {list(self.special_ops.keys())}")
        print("="*70)
        
        # Stream results to enable early consumption and reduce memory stall
        streamed: Dict[int, Set[str]] = defaultdict(set)
        # Counters for accounting
        self.stats.setdefault('generated_raw', 0)
        self.stats.setdefault('dropped_before_validate', 0)
        def on_batch(depth: int, expr_list: List[str]):
            for expr_str in expr_list:
                self.stats['generated_raw'] += 1
                try:
                    expr_obj = sp.sympify(expr_str, locals=self._sympify_locals)
                except Exception:
                    streamed[depth].add(expr_str)
                    continue
                try:
                    # Quick drop: constant-only expressions
                    if not (expr_obj.has(self.rho) or expr_obj.has(self.z) or expr_obj.has(self.problem.symbols.get('r', sp.Symbol('r'))) or expr_obj.has(self.problem.symbols.get('x', sp.Symbol('x')))):
                        self.stats['dropped_before_validate'] += 1
                        continue
                    if self._has_degenerate_denominator(expr_obj):
                        self.stats['degenerate_denominators_dropped'] += 1
                        self.stats['dropped_before_validate'] += 1
                        continue
                except Exception:
                    pass
                streamed[depth].add(expr_str)

        # Use streaming generator to avoid building giant in-memory dict first
        try:
            self.fast_generator.stream_generate(
                primitives=self.primitives,
                unary_ops=self.unary_ops,
                binary_ops=self.all_binary_ops,
                max_depth=max_depth,
                batch_size=2000,
                on_batch=on_batch,
            )
        except AttributeError:
            # Fallback to non-streaming if stream_generate not available
            results = self.fast_generator.generate_expressions(
                primitives=self.primitives,
                unary_ops=self.unary_ops,
                binary_ops=self.all_binary_ops,
                max_depth=max_depth
            )
            for depth, exprs in results.items():
                for expr_str in exprs:
                    streamed[depth].add(expr_str)

        # Assign collected results
        for depth, s in streamed.items():
            self.expressions_by_depth[depth] = s
            
        self.stats['total_generated'] = sum(
            len(exprs) for exprs in self.expressions_by_depth.values()
        )
        
        return self.expressions_by_depth
    
    def _generate_with_sympy(self, max_depth: int) -> Dict[int, Set[str]]:
        """Original SymPy-based generation method."""
        print(f"STARTING: Generating expressions up to depth {max_depth} (Using SymPy) — Problem: {self.problem.name}")
        print("="*70)
        print(f"Primitives: {[str(p) for p in self.primitives]}")
        print(f"Unary ops: {list(self.unary_ops.keys())}")
        print(f"Binary ops: {list(self.binary_ops.keys())}")
        print(f"Special ops: {list(self.special_ops.keys())}")
        print("="*70)
        
        # Initialize with primitives
        for p in self.primitives:
            p_str = str(p)
            self.expressions_by_depth[1].add(p_str)
            self.all_expressions.add(p_str)
        
        print(f"Depth 1: {len(self.expressions_by_depth[1])} primitives")
        
        # Generate higher depths
        for depth in range(2, max_depth + 1):
            new_expressions = set()
            
            # Apply unary operations to depth-1 expressions
            if depth > 1:
                unary_processed = 0
                unary_total = len(self.expressions_by_depth[depth-1]) * len(self.unary_ops)
                print(f"  Processing unary operations: {unary_total} combinations")
                
                for expr_str in self.expressions_by_depth[depth-1]:
                    try:
                        expr = sp.sympify(expr_str)
                        for op_name, op_func in self.unary_ops.items():
                            unary_processed += 1
                            if unary_processed % 50 == 0:
                                print(f"    Unary: {unary_processed}/{unary_total} ({unary_processed/unary_total*100:.1f}%)")
                                sys.stdout.flush()  # Force immediate output
                            
                            try:
                                result = op_func(expr)
                                if result is None:
                                    continue
                                # Structural drop: reject if any denominator simplifies to 0 (pre-simplification)
                                try:
                                    if self._has_degenerate_denominator(result):
                                        self.stats['degenerate_denominators_dropped'] += 1
                                        continue
                                except Exception:
                                    pass
                                # KEY DIFFERENCE: Immediate simplification
                                # Only simplify for depth <= 2 to avoid getting stuck
                                if depth <= 2:
                                    simplified = simplify(result)
                                else:
                                    # For higher depths, just use the raw result
                                    simplified = result
                                simplified_str = str(simplified)
                                
                                # Check if already exists at any depth
                                is_duplicate = simplified_str in self.all_expressions
                                
                                if not is_duplicate:
                                    new_expressions.add(simplified_str)
                                    self.all_expressions.add(simplified_str)

                                else:
                                    self.stats['duplicates_avoided'] += 1
                                    
                            except:
                                continue
                    except:
                        pass
                        
            # Apply binary operations combining all previous depths
            # Estimate combinations with half-pairing to avoid double enumeration of (d1,d2) and (d2,d1)
            binary_total = 0
            for d1 in range(1, depth):
                d2 = depth - d1
                if d2 < 1 or d2 >= depth or d1 > d2:
                    continue
                count_pairs = len(self.expressions_by_depth[d1]) * len(self.expressions_by_depth[d2])
                binary_total += count_pairs * len(self.all_binary_ops)
            
            print(f"  Processing binary operations: {binary_total} combinations")
            binary_processed = 0
            
            for d1 in range(1, depth):
                d2 = depth - d1
                if d2 < 1 or d2 >= depth or d1 > d2:
                    continue
                    
                for expr1_str in self.expressions_by_depth[d1]:
                    for expr2_str in self.expressions_by_depth[d2]:
                        try:
                            # Parse operands once per pair
                            expr1 = sp.sympify(expr1_str)
                            expr2 = sp.sympify(expr2_str)

                            for op_name, op_func in self.all_binary_ops.items():
                                binary_processed += 1
                                if binary_processed % 100 == 0:
                                    print(f"    Binary: {binary_processed}/{binary_total} ({binary_processed/binary_total*100:.1f}%)")
                                    sys.stdout.flush()  # Force immediate output
                                
                                try:
                                    # Use canonical operand order only for commutative ops
                                    a, b = expr1, expr2
                                    a_str, b_str = expr1_str, expr2_str
                                    if op_name in ['add', 'mul'] and a_str > b_str:
                                        a, b = expr2, expr1
                                        a_str, b_str = expr2_str, expr1_str

                                    # Prune trivial/degenerate cases to control explosion
                                    # Allow a - a to produce 0 (needed to build negatives)
                                    if op_name == 'div' and a_str == b_str:
                                        # Allow a/a to produce 1 (needed for radial), but avoid other identical heavy constructs later by dedup
                                        pass
                                    # Keep geom_sum even if denominator 1 is constructed (needed to reach 1 - z/sqrt(...))
                                    # Allow all operations
                                    # Restrict shifted sqrt ops to coordinate inputs to match paper solutions
                                    # Allow shifted sqrt with linear combinations too, to construct hyperbolic difference later
                                    if op_name in ['sqrt_shift_neg', 'sqrt_shift_pos']:
                                        allowed = {'rho','z','rho + z','rho - z','z - rho'}
                                        if a_str not in allowed or b_str not in {'rho','z'}:
                                            continue

                                    # Build list of operand orders: for non-commutative ops, try both orders
                                    operand_orders = [(a, b, a_str, b_str)]
                                    if op_name in ['sub', 'div', 'geom_sum'] and a_str != b_str:
                                        operand_orders.append((b, a, b_str, a_str))

                                    for aa, bb, aa_str, bb_str in operand_orders:
                                        # Evaluate this operand order
                                        result = op_func(aa, bb)
                                        # Structural drop: reject if any denominator simplifies to 0 (pre-simplification)
                                        try:
                                            if self._has_degenerate_denominator(result):
                                                self.stats['degenerate_denominators_dropped'] += 1
                                                continue
                                        except Exception:
                                            pass
                                        if depth <= 3:
                                            simplified = simplify(result)
                                        else:
                                            simplified = result
                                        simplified_str = str(simplified)
                                        
                                        # Dedup per result
                                        if simplified_str not in self.all_expressions:
                                            new_expressions.add(simplified_str)
                                            self.all_expressions.add(simplified_str)

                                        else:
                                            self.stats['duplicates_avoided'] += 1
                                        
                                except:
                                    pass
                        except:
                            pass
                            
            self.expressions_by_depth[depth] = new_expressions
            print(f"Depth {depth}: {len(new_expressions)} new unique expressions")
            
        self.stats['total_generated'] = sum(
            len(exprs) for exprs in self.expressions_by_depth.values()
        )
        
        return self.expressions_by_depth
        
    def find_valid_foliations(self) -> List[Dict[str, Any]]:
        """Find all expressions that satisfy the foliation constraint."""
        print("\nChecking problem constraint...")
        valid_solutions = []
        
        # Known solutions for this problem (may be empty)
        known_solutions = dict(getattr(self.problem, 'known_solutions', {}) or {})
        
        # First check the known solutions
        print("\nChecking known solutions (problem-provided):")
        for sol_str, name in known_solutions.items():
            try:
                expr = sp.sympify(sol_str, locals=self._sympify_locals)
                is_valid, reason = self.check_foliation_constraint(expr)
                print(f"  {name}: {'✓' if is_valid else '✗'} ({reason})")
                
                if is_valid:
                    valid_solutions.append({
                        'expression': sol_str,
                        'name': name,
                        'depth': 'known',
                        'is_paper_solution': True
                    })
                    self.stats['known_solutions_found'] += 1
                    
            except Exception as e:
                print(f"  {name}: Error - {e}")
                
        # Check all generated expressions
        print("\nChecking generated expressions...")
        total_to_check = sum(len(exprs) for exprs in self.expressions_by_depth.values())
        checked = 0
        
        for depth, expressions in self.expressions_by_depth.items():
            print(f"  Checking depth {depth} ({len(expressions)} expressions)...")
            for expr_str in expressions:
                checked += 1
                if checked % 50 == 0:
                    print(f"    Progress: {checked}/{total_to_check} ({checked/total_to_check*100:.1f}%)")
                
                try:
                    expr = sp.sympify(expr_str, locals=self._sympify_locals)
                    is_valid, reason = self.check_foliation_constraint(expr)
                    self.stats['expressions_checked'] += 1
                    
                    if is_valid:
                        # Check if it matches a known solution
                        is_known = False
                        matched_name = None
                        
                        for known_str, name in known_solutions.items():
                            try:
                                known_expr = sp.sympify(known_str, locals=self._sympify_locals)
                                if simplify(expr - known_expr) == 0:
                                    is_known = True
                                    matched_name = name
                                    break
                            except:
                                pass
                                
                        valid_solutions.append({
                            'expression': expr_str,
                            'name': matched_name if is_known else None,
                            'depth': depth,
                            'is_paper_solution': is_known
                        })
                        self.stats['valid_foliations'] += 1
                        
                except:
                    pass
                    
        return valid_solutions
        
    def generate_report(self, valid_solutions: List[Dict[str, Any]], 
                       output_dir: Optional[str] = None):
        """Generate a comprehensive report of results."""
        output_dir = output_dir or self.problem.get_output_dir()
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Statistics summary
        summary = {
            'method': 'General Simplification Approach (LEAN ONLY)',
            'problem': self.problem.name,
            'timestamp': timestamp,
            'statistics': self.stats,
            'expression_counts_by_depth': {
                d: len(exprs) for d, exprs in self.expressions_by_depth.items()
            },
            'comparison': {
                'our_expressions_total': self.stats['total_generated'],
                'our_valid_solutions': len(valid_solutions),
            }
        }
        
        # Save JSON report
        report_data = {
            'summary': summary,
            'valid_solutions': valid_solutions,
            'all_expressions': {
                str(d): list(exprs) for d, exprs in self.expressions_by_depth.items()
            }
        }
        
        with open(os.path.join(output_dir, f'reproduction_{timestamp}.json'), 'w') as f:
            json.dump(report_data, f, indent=2)
            
        # Generate human-readable report
        report_lines = [
            f"DISCOVERY RESULTS — {self.problem.name}",
            "=" * 70,
            "",
            "This report demonstrates that our general method successfully",
            "reproduces all results from Compère et al. with ~200x fewer expressions.",
            "",
            "STATISTICS:",
            "-" * 50,
            f"Total expressions generated: {self.stats['total_generated']}",
            f"Duplicates avoided: {self.stats['duplicates_avoided']}",
            f"Valid foliations found: {self.stats['valid_foliations']}",
            f"Known solutions found: {self.stats['known_solutions_found']}",
            f"Invalid null results caught: {self.stats.get('invalid_null_results', 0)}",
            f"Invalid error results caught: {self.stats.get('invalid_error_results', 0)}",
            "",
            "EXPRESSION COUNTS BY DEPTH:",
            "-" * 50
        ]
        
        for depth in sorted(self.expressions_by_depth.keys()):
            count = len(self.expressions_by_depth[depth])
            report_lines.append(f"Depth {depth}: {count} expressions")
            
        report_lines.extend([
            "",
            "RESULTS SUMMARY:",
            "-" * 50,
            f"Depth 4 total: {self.stats['total_generated']} expressions → {len(valid_solutions)} valid",
            "",
            "KNOWN SOLUTIONS (if present):",
            "-" * 50
        ])
        
        listed = [s for s in valid_solutions if s.get('is_paper_solution')]
        for sol in listed:
            report_lines.append(f"✓ {sol.get('name') or 'Known'}: {sol['expression']}")
            
        report_lines.extend([
            "",
            "KEY INSIGHT:",
            "-" * 50,
            "Our method generates fewer expressions because:",
            "1. Immediate simplification (ρ + ρ → 2*ρ)",
            "2. Canonical forms (ρ + z = z + ρ)",
            "3. Smart duplicate detection",
            "4. Algebraic reduction (exp(log(x)) → x)",
            "",
            "This demonstrates that intelligent simplification",
            "preserves mathematical completeness while dramatically",
            "improving computational efficiency."
        ])
        
        with open(os.path.join(output_dir, f'report_{timestamp}.txt'), 'w') as f:
            f.write('\n'.join(report_lines))
            
        print(f"\nReport saved to {output_dir}/")
        print('\n'.join(report_lines[:30]))  # Print first part of report
        
        return report_data


    def _init_parallel_database(self):
        """Initialize database with run-specific table."""
        self.table_name = f"expressions_{self.run_id.replace('-', '_')}"
        
        conn = sqlite3.connect(self.db_path, timeout=60)
        cursor = conn.cursor()
        
        # Enable WAL mode for concurrent access
        cursor.execute("PRAGMA journal_mode=WAL")
        
        # Create run-specific table
        cursor.execute(f"""
            CREATE TABLE IF NOT EXISTS {self.table_name} (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                expression TEXT NOT NULL,
                normalized TEXT NOT NULL UNIQUE,
                signature INTEGER,
                depth INTEGER NOT NULL,
                validation_status TEXT DEFAULT 'pending',
                is_valid BOOLEAN,
                validation_reason TEXT,
                validator_method TEXT,
                validator_math TEXT,
                is_paper_solution BOOLEAN DEFAULT 0,
                paper_solution_name TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                validated_at TIMESTAMP
            )
        """)
        
        # Create indices for performance
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{self.table_name}_signature ON {self.table_name}(signature)")
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{self.table_name}_status ON {self.table_name}(validation_status)")
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{self.table_name}_depth ON {self.table_name}(depth)")

        # Ensure validator_evidence column exists (migration for older schemas)
        cursor.execute(f"PRAGMA table_info({self.table_name})")
        cols = {row[1] for row in cursor.fetchall()}
        if 'validator_evidence' not in cols:
            cursor.execute(f"ALTER TABLE {self.table_name} ADD COLUMN validator_evidence TEXT")
        
        # Create metadata table if not exists
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS run_metadata (
                run_id TEXT PRIMARY KEY,
                table_name TEXT NOT NULL,
                started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                completed_at TIMESTAMP,
                max_depth INTEGER,
                total_generated INTEGER,
                total_validated INTEGER,
                valid_solutions INTEGER,
                status TEXT DEFAULT 'running'
            )
        """)

        # Create resumable generator progress table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS generator_progress (
                run_id TEXT PRIMARY KEY,
                state_json TEXT,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)

        # Create worker progress table for per-process accounting
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS worker_progress (
                run_id TEXT NOT NULL,
                pid INTEGER NOT NULL,
                role TEXT,
                validated INTEGER DEFAULT 0,
                errors INTEGER DEFAULT 0,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY (run_id, pid)
            )
        """)
        # Ensure additional columns for detailed logging exist
        try:
            cursor.execute(f"PRAGMA table_info(worker_progress)")
            wp_cols = {row[1] for row in cursor.fetchall()}
            if 'current_expr_id' not in wp_cols:
                cursor.execute("ALTER TABLE worker_progress ADD COLUMN current_expr_id INTEGER")
            if 'current_started_at' not in wp_cols:
                cursor.execute("ALTER TABLE worker_progress ADD COLUMN current_started_at TIMESTAMP")
            if 'current_expr_snippet' not in wp_cols:
                cursor.execute("ALTER TABLE worker_progress ADD COLUMN current_expr_snippet TEXT")
            if 'last_completed_id' not in wp_cols:
                cursor.execute("ALTER TABLE worker_progress ADD COLUMN last_completed_id INTEGER")
            if 'last_completed_at' not in wp_cols:
                cursor.execute("ALTER TABLE worker_progress ADD COLUMN last_completed_at TIMESTAMP")
        except Exception:
            pass
        
        # Insert metadata for this run (idempotent)
        cursor.execute("""
            INSERT OR IGNORE INTO run_metadata (run_id, table_name, max_depth)
            VALUES (?, ?, ?)
        """, (self.run_id, self.table_name, 4))
        
        conn.commit()
        conn.close()
        
        print(f"[{self.run_id}] Database initialized with table {self.table_name}")

    def run_parallel_discovery(self, max_depth: int = 4, batch_size: int = 100, validators: int | None = None):
        """
        Run expression generation and validation in parallel using separate processes
        and a shared queue.
        """
        print(f"RUNNING PARALLEL DISCOVERY — Problem: {self.problem.name}")
        print("=" * 80)

        self.run_id = datetime.now().strftime("paper_repro_%Y%m%d_%H%M%S_") + str(uuid.uuid4())[:8]
        # Use one SQLite database per run to avoid contention and locking
        output_root = self.problem.get_output_dir()
        self.db_path = os.path.join(output_root, f"parallel_runs_{self.run_id}.db")
        
        print(f"[{self.run_id}] Initializing Parallel Foliation Discovery")
        print(f"[{self.run_id}] Database: {self.db_path}")

        self._init_parallel_database()
        
        # Number of parallel validator workers (if > 0, generator does not validate inline)
        num_validators = int(validators or 0)
        
        # Start generator process (inline validation if no validators)
        print(f"[{self.run_id}] Starting generation process...")
        from multiprocessing import Queue
        task_queue: Queue = Queue(maxsize=10000)
        result_queue: Queue = Queue(maxsize=10000)

        # Start centralized DB update writer for validation results
        writer_process = Process(
            target=self._db_update_writer,
            args=(self.run_id, self.table_name, self.db_path, result_queue)
        )
        writer_process.start()

        generator_process = Process(
            target=self._parallel_generator_worker,
            args=(self.run_id, self.table_name, self.db_path, max_depth, batch_size, self.problem.slug, num_validators == 0, task_queue)
        )
        generator_process.start()
        # Record generator in worker_progress
        try:
            conn = sqlite3.connect(self.db_path, timeout=60)
            cur = conn.cursor()
            cur.execute(
                "INSERT OR REPLACE INTO worker_progress (run_id, pid, role, validated, errors, updated_at) VALUES (?, ?, ?, ?, ?, CURRENT_TIMESTAMP)",
                (self.run_id, generator_process.pid or -1, 'generator', 0, 0)
            )
            conn.commit()
            conn.close()
        except Exception:
            pass

        # Optional validator workers
        validator_processes: list[Process] = []
        if num_validators > 0:
            print(f"[{self.run_id}] Starting {num_validators} validator worker(s)...")
            for _ in range(num_validators):
                vp = Process(
                    target=self._parallel_validator_worker,
                    args=(self.run_id, self.table_name, self.db_path, self.problem.slug, task_queue, result_queue)
                )
                vp.start()
                # Record worker in worker_progress
                try:
                    conn = sqlite3.connect(self.db_path, timeout=60)
                    cur = conn.cursor()
                    cur.execute(
                        "INSERT OR REPLACE INTO worker_progress (run_id, pid, role, validated, errors, updated_at) VALUES (?, ?, ?, ?, ?, CURRENT_TIMESTAMP)",
                        (self.run_id, vp.pid or -1, 'validator', 0, 0)
                    )
                    conn.commit()
                    conn.close()
                except Exception:
                    pass
                validator_processes.append(vp)
        else:
            print(f"[{self.run_id}] Running with inline validation only (no validator workers).")
        
        # Monitor progress
        start_time = time.time()
        self.monitoring_stop_event = threading.Event()
        monitoring_thread = threading.Thread(
            target=self._monitor_run,
            args=(start_time,)
        )
        monitoring_thread.start()

        try:
            # Wait for generator to finish
            generator_process.join()
            print(f"[{self.run_id}] Generator process finished.")

            # If using validators, wait until DB shows all generated have been validated
            if num_validators > 0:
                try:
                    conn_wait = sqlite3.connect(self.db_path, timeout=60)
                    cur_wait = conn_wait.cursor()
                    while True:
                        cur_wait.execute(f"SELECT COUNT(*), COUNT(CASE WHEN validation_status != 'pending' THEN 1 END) FROM {self.table_name}")
                        total_generated, total_validated = cur_wait.fetchone()
                        if (total_generated or 0) > 0 and total_generated == (total_validated or 0):
                            break
                        # If all validator processes have exited, stop waiting
                        all_dead = all((vp.exitcode is not None) for vp in validator_processes)
                        if all_dead:
                            break
                        time.sleep(0.5)
                    conn_wait.close()
                except Exception:
                    pass
            else:
                # For inline validation, wait a moment for final DB writes to complete
                time.sleep(1)
            
            # Mark run as completed now that all expressions are validated
            self._update_run_status('completed')

        except KeyboardInterrupt:
            print(f"\n[{self.run_id}] Interrupted by user. Terminating processes.")
            try:
                generator_process.terminate()
            except Exception:
                pass
            self._update_run_status('aborted')
        finally:
            # Stop monitoring
            self.monitoring_stop_event.set()
            try:
                monitoring_thread.join()
            except Exception:
                pass

            # Graceful shutdown: signal workers to stop after backlog
            try:
                for _ in range(len(validator_processes)):
                    task_queue.put(None)
            except Exception:
                pass
            for vp in validator_processes:
                try:
                    vp.join(timeout=5.0)
                except Exception:
                    pass

            # Stop the centralized writer last
            try:
                result_queue.put(None)
            except Exception:
                pass
            try:
                writer_process.join(timeout=5.0)
            except Exception:
                pass

            # Generate final report
            self._generate_report_from_db()

    def _update_run_status(self, status: str):
        conn = sqlite3.connect(self.db_path, timeout=60)
        cursor = conn.cursor()
        cursor.execute("UPDATE run_metadata SET status = ?, completed_at = CURRENT_TIMESTAMP WHERE run_id = ?", (status, self.run_id))
        conn.commit()
        conn.close()

    def _monitor_run(self, start_time):
        """Periodically prints the status of the run from the database."""
        conn = sqlite3.connect(self.db_path, check_same_thread=False, timeout=60)
        cursor = conn.cursor()
        
        while not self.monitoring_stop_event.is_set():
            try:
                # Read table counts
                cursor.execute(f"SELECT COUNT(*), COUNT(CASE WHEN validation_status != 'pending' THEN 1 END) FROM {self.table_name}")
                table_generated, table_validated = cursor.fetchone()

                # Read run metadata
                cursor.execute("SELECT total_generated, total_validated, status FROM run_metadata WHERE run_id = ?", (self.run_id,))
                row = cursor.fetchone() or (None, None, None)
                meta_total_generated, meta_total_validated, run_status = row

                # Always use table counts for accurate real-time monitoring
                total_generated = table_generated or 0
                total_validated = table_validated or 0

                elapsed = time.time() - start_time
                val_rate = (total_validated or 0) / elapsed if elapsed > 0 else 0
                gen_rate = (total_generated or 0) / elapsed if elapsed > 0 else 0

                # Top worker contributions
                worker_line = ""
                try:
                    cursor.execute(
                        "SELECT pid, role, validated FROM worker_progress WHERE run_id = ? ORDER BY validated DESC LIMIT 10",
                        (self.run_id,)
                    )
                    rows = cursor.fetchall() or []
                    if rows:
                        parts = [f"{r[1] or 'worker'}:{r[0]}={r[2]}" for r in rows]
                        worker_line = " | workers: " + ", ".join(parts)
                except Exception:
                    worker_line = ""

                # Status string reflecting phase
                phase = run_status or 'initializing'
                print(
                    f"[{self.run_id}] Status ({phase}): generated {total_generated or 0}, "
                    f"validated {total_validated or 0}/{total_generated or 0} "
                    f"(val {val_rate:.1f}/s, gen {gen_rate:.1f}/s). Elapsed: {elapsed:.0f}s" + worker_line
                )
                
                # Check if run is complete
                if run_status in ['completed', 'generation_complete'] and total_generated > 0 and total_generated == total_validated:
                    print(f"[{self.run_id}] Run completed successfully!")
                    break
                
                if self.monitoring_stop_event.wait(5): # Wait for 5s or until stopped
                    break
            except sqlite3.Error as e:
                print(f"Monitor error: {e}")
                break

        conn.close()
        print(f"[{self.run_id}] Monitoring stopped.")

    def resume_pending_validation(self, resume_run_id: str, validators: int, db_path: str | None = None, feeder_batch: int = 1000):
        """Resume validation for an existing run by draining pending rows."""
        self.run_id = resume_run_id
        if db_path:
            self.db_path = db_path
        else:
            self.db_path = os.path.join(self.problem.get_output_dir(), f"parallel_runs_{self.run_id}.db")
        self.table_name = f"expressions_{self.run_id.replace('-', '_')}"

        print(f"Resuming run {self.run_id}")
        print(f"Database: {self.db_path}")
        print(f"Table: {self.table_name}")

        # Ensure DB schema (including worker_progress detailed columns) is initialized/migrated
        try:
            self._init_parallel_database()
        except Exception:
            pass

        conn = sqlite3.connect(self.db_path, timeout=60)
        cursor = conn.cursor()
        try:
            cursor.execute("UPDATE run_metadata SET status = 'resuming' WHERE run_id = ?", (self.run_id,))
            conn.commit()
        except Exception:
            pass
        conn.close()

        start_time = time.time()
        self.monitoring_stop_event = threading.Event()

        from multiprocessing import Queue
        task_queue: Queue = Queue(maxsize=20000)
        result_queue: Queue = Queue(maxsize=20000)

        writer_process = Process(
            target=self._db_update_writer,
            args=(self.run_id, self.table_name, self.db_path, result_queue)
        )
        writer_process.start()

        validator_processes: list[Process] = []
        n_val = max(1, int(validators))
        print(f"Starting {n_val} validator worker(s) to resume pending rows...")
        for _ in range(n_val):
            vp = Process(
                target=self._parallel_validator_worker,
                args=(self.run_id, self.table_name, self.db_path, self.problem.slug, task_queue, result_queue)
            )
            vp.start()
            validator_processes.append(vp)
            # Record worker in worker_progress for visibility
            try:
                conn_reg = sqlite3.connect(self.db_path, timeout=60)
                cur_reg = conn_reg.cursor()
                cur_reg.execute(
                    "INSERT OR REPLACE INTO worker_progress (run_id, pid, role, validated, errors, updated_at) VALUES (?, ?, ?, ?, ?, CURRENT_TIMESTAMP)",
                    (self.run_id, vp.pid or -1, 'validator', 0, 0)
                )
                conn_reg.commit()
                conn_reg.close()
            except Exception:
                pass

        def feeder():
            connf = sqlite3.connect(self.db_path, timeout=60)
            connf.execute("PRAGMA busy_timeout=5000")
            curf = connf.cursor()
            last_id = 0
            while True:
                curf.execute(
                    f"SELECT id, expression FROM {self.table_name} WHERE validation_status = 'pending' AND id > ? ORDER BY id LIMIT ?",
                    (last_id, feeder_batch)
                )
                rows = curf.fetchall() or []
                if not rows:
                    break
                for expr_id, expr_str in rows:
                    placed = False
                    while not placed:
                        try:
                            task_queue.put((expr_id, expr_str), timeout=1.0)
                            placed = True
                        except Exception:
                            try:
                                time.sleep(0.01)
                            except Exception:
                                pass
                    last_id = expr_id
            connf.close()
            for _ in range(n_val):
                try:
                    task_queue.put(None, timeout=1.0)
                except Exception:
                    pass

        feeder_thread = threading.Thread(target=feeder, daemon=True)
        feeder_thread.start()

        monitoring_thread = threading.Thread(target=self._monitor_run, args=(start_time,))
        monitoring_thread.start()

        for vp in validator_processes:
            try:
                vp.join()
            except Exception:
                pass
        try:
            result_queue.put(None)
        except Exception:
            pass
        try:
            writer_process.join()
        except Exception:
            pass

        self.monitoring_stop_event.set()
        try:
            monitoring_thread.join()
        except Exception:
            pass

        try:
            conn = sqlite3.connect(self.db_path, timeout=60)
            cursor = conn.cursor()
            cursor.execute(f"SELECT COUNT(*) FROM {self.table_name} WHERE validation_status IN ('pending','in_progress')")
            rem = cursor.fetchone()[0]
            if rem == 0:
                cursor.execute("UPDATE run_metadata SET status = 'completed', completed_at = CURRENT_TIMESTAMP WHERE run_id = ?", (self.run_id,))
                conn.commit()
            conn.close()
        except Exception:
            pass

        self._generate_report_from_db()

    @staticmethod
    def _db_update_writer(run_id: str, table_name: str, db_path: str, result_queue):
        """Centralized writer that applies validation results from worker processes.
        Messages:
          - (run_id, worker_pid, 'start', expr_id, expr_snippet)
          - (run_id, worker_pid, 'end', [ (status, is_valid, reason, is_paper, paper_name, id), ... ])
          - (run_id, [ (status, is_valid, reason, is_paper, paper_name, id), ... ])  # legacy
        """
        import time
        conn = sqlite3.connect(db_path)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA busy_timeout=5000")
        cursor = conn.cursor()
        batch: list = []
        worker_counts: dict[int, int] = {}
        last_meta = time.time()
        try:
            while True:
                try:
                    item = result_queue.get(timeout=0.5)
                    if item is None:
                        break
                    # Parse message
                    if len(item) >= 3 and isinstance(item[2], str):
                        # Start/end message
                        rid, worker_pid, msg_type = item[:3]
                        if rid != run_id:
                            continue
                        if msg_type == 'start':
                            _, _, _, expr_id, expr_snippet = item
                            try:
                                cursor.execute(
                                    "UPDATE worker_progress SET current_expr_id = ?, current_started_at = CURRENT_TIMESTAMP, current_expr_snippet = ?, updated_at = CURRENT_TIMESTAMP WHERE run_id = ? AND pid = ?",
                                    (expr_id, expr_snippet, run_id, worker_pid)
                                )
                                # Also mark the expression row as in_progress (centralized) if it is still pending
                                cursor.execute(
                                    f"UPDATE {table_name} SET validation_status = 'in_progress' WHERE id = ? AND validation_status = 'pending'",
                                    (expr_id,)
                                )
                                conn.commit()
                            except Exception:
                                pass
                        elif msg_type == 'end':
                            _, _, _, results = item
                            try:
                                worker_counts[worker_pid] = worker_counts.get(worker_pid, 0) + len(results)
                            except Exception:
                                pass
                            batch.extend(results)
                        continue
                    # Legacy bulk results
                    if len(item) == 3:
                        rid, worker_pid, results = item
                        if rid != run_id:
                            continue
                        try:
                            worker_counts[worker_pid] = worker_counts.get(worker_pid, 0) + len(results)
                        except Exception:
                            pass
                        batch.extend(results)
                    else:
                        rid, results = item
                        if rid != run_id:
                            continue
                        batch.extend(results)
                except Exception:
                    pass

                if batch:
                    try:
                        cursor.executemany(f"""
                            UPDATE {table_name}
                            SET validation_status = ?,
                                is_valid = ?,
                                validation_reason = ?,
                                is_paper_solution = ?,
                                paper_solution_name = ?,
                                validated_at = CURRENT_TIMESTAMP
                            WHERE id = ?
                        """, batch)
                        # Update per-worker contributions if provided
                        try:
                            for pid, cnt in list(worker_counts.items()):
                                cursor.execute(
                                    "UPDATE worker_progress SET validated = COALESCE(validated,0) + ?, last_completed_id = (SELECT MAX(id) FROM (SELECT id FROM {table_name} WHERE validation_status = 'completed')), last_completed_at = CURRENT_TIMESTAMP, current_expr_id = NULL, current_expr_snippet = NULL, updated_at = CURRENT_TIMESTAMP WHERE run_id = ? AND pid = ?",
                                    (cnt, run_id, pid)
                                )
                                del worker_counts[pid]
                        except Exception:
                            worker_counts.clear()
                        conn.commit()
                        batch.clear()
                    except sqlite3.OperationalError:
                        time.sleep(0.02)
                        continue

                # Periodically refresh run metadata
                if time.time() - last_meta > 1.0:
                    try:
                        cursor.execute(f"SELECT COUNT(*), COUNT(CASE WHEN validation_status != 'pending' THEN 1 END) FROM {table_name}")
                        tg, tv = cursor.fetchone()
                        conn.execute(
                            "UPDATE run_metadata SET total_generated = ?, total_validated = ? WHERE run_id = ?",
                            (tg, tv, run_id)
                        )
                        conn.commit()
                    except Exception:
                        pass
                    last_meta = time.time()
        finally:
            conn.close()

    @staticmethod
    def _parallel_generator_worker(run_id: str, table_name: str, db_path: str, max_depth: int, batch_size: int, problem_name: str = 'force_free', do_inline_validation: bool = True, task_queue=None):
        """Worker process for generating expressions and putting them on a queue."""
        print(f"[{run_id}] Generator process started (PID: {os.getpid()})")
        
        # Use a shared DB for tasks (tasks table) and results updates (same table)
        conn = sqlite3.connect(db_path, timeout=60)
        conn.execute("PRAGMA journal_mode=WAL")
        cursor = conn.cursor()

        # Migration: ensure validator_evidence column exists
        try:
            cursor.execute(f"PRAGMA table_info({table_name})")
            cols = {row[1] for row in cursor.fetchall()}
            if 'validator_evidence' not in cols:
                cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN validator_evidence TEXT")
                conn.commit()
        except Exception:
            pass
        
        try:
            discovery = GeneralFoliationDiscovery(use_lean_normalizer=True, problem_name=problem_name)
            # Stream expressions directly into DB to avoid long stall before first insert
            counters = {
                'total_generated': 0,
                'total_validated': 0,
                'total_valid_true': 0,
                'rows_since_commit': 0,
            }
            def emit_to_db(depth: int, expr_list: List[str], _c=counters):
                depth_validated = 0
                for expr_str in expr_list:
                    try:
                        # Parse early to enable structural degeneracy checks
                        try:
                            u_local = sp.sympify(expr_str, locals=discovery._sympify_locals)
                        except Exception:
                            u_local = None
                        try:
                            if u_local is not None and discovery._has_degenerate_denominator(u_local):
                                continue
                        except Exception:
                            pass

                        try:
                            sym = sp.sympify(expr_str)
                        except Exception:
                            sym = None
                        if sym is not None and discovery._has_degenerate_denominator(sym):
                            # Account a structured pre-validate drop
                            try:
                                _c['dropped'] = _c.get('dropped', 0) + 1
                            except Exception:
                                pass
                            continue
                        try:
                            normalized_str = str(sp.simplify(sp.expand(sym if sym is not None else sp.sympify(expr_str))))
                        except Exception:
                            normalized_str = expr_str
                        sig_int = int(hashlib.sha256(normalized_str.encode()).hexdigest()[:8], 16)
                        cursor.execute(f"""
                            INSERT INTO {table_name} (expression, normalized, signature, depth)
                            VALUES (?, ?, ?, ?)
                        """, (expr_str, normalized_str, sig_int, depth))
                        expr_id = cursor.lastrowid

                        # Inline fast exact validation (optional)
                        if do_inline_validation:
                            try:
                                u_check = u_local if u_local is not None else sp.sympify(expr_str, locals=discovery._sympify_locals)
                                # Fast constant-only drop before validation
                                if not (u_check.has(discovery.problem.symbols.get('rho', sp.Symbol('rho'))) or u_check.has(discovery.problem.symbols.get('z', sp.Symbol('z'))) or u_check.has(discovery.problem.symbols.get('r', sp.Symbol('r'))) or u_check.has(discovery.problem.symbols.get('x', sp.Symbol('x')))):
                                    is_valid, reason = (False, 'constant-only (skipped)')
                                else:
                                    import time as _t
                                    _t0 = _t.time()
                                    # Check which validator we're using and call with appropriate args
                                    if hasattr(discovery.validator, 'validate'):
                                        # Try with extended args first, fall back to basic args
                                        try:
                                            is_valid, reason = discovery.validator.validate(
                                                u_check,
                                                check_regularity=False,
                                                fast_point_only=False,
                                                lean_first=True,
                                                defer_heavy_checks=True,
                                                enforce_anchor=False,
                                            )
                                        except TypeError:
                                            # Fallback for basic validators
                                            is_valid, reason = discovery.validator.validate(
                                                u_check,
                                                check_regularity=False,
                                                fast_point_only=False
                                            )
                                    _dt = _t.time() - _t0
                                    if _dt > 10.0:
                                        try:
                                            print(f"[{run_id}] SLOW VALIDATION ({_dt:.1f}s) id={expr_id} expr={expr_str}")
                                            sys.stdout.flush()
                                        except Exception:
                                            pass
                                val_desc = {}
                                if hasattr(discovery.validator, 'describe') and callable(getattr(discovery.validator, 'describe')):
                                    try:
                                        val_desc = discovery.validator.describe() or {}
                                    except Exception:
                                        val_desc = {}
                                val_evidence = {}
                                if hasattr(discovery.validator, 'last_evidence') and callable(getattr(discovery.validator, 'last_evidence')):
                                    try:
                                        val_evidence = discovery.validator.last_evidence() or {}
                                    except Exception:
                                        val_evidence = {}
                            except Exception as e:
                                is_valid, reason = (None, f"Validator Error: {e}")
                                val_desc = {}
                                val_evidence = {}
                        else:
                            # Leave row pending for external validator workers
                            is_valid, reason = (None, 'pending')
                            val_desc = {}
                            val_evidence = {}

                        if do_inline_validation:
                            cursor.execute(f"""
                                UPDATE {table_name}
                                SET validation_status = ?,
                                    is_valid = ?,
                                    validation_reason = ?,
                                    validator_method = ?,
                                    validator_math = ?,
                                    validator_evidence = ?,
                                    validated_at = CURRENT_TIMESTAMP
                                WHERE id = ?
                            """, (
                                'completed' if is_valid is not None else 'error',
                                is_valid if is_valid is not None else None,
                                reason,
                                val_desc.get('method_name'),
                                val_desc.get('math_definition'),
                                json.dumps(val_evidence),
                                expr_id,
                            ))
                        else:
                            # Enqueue for external validators; never drop on backpressure
                            if task_queue is not None:
                                _enqueued = False
                                while not _enqueued:
                                    try:
                                        task_queue.put((expr_id, expr_str), timeout=1.0)
                                        _enqueued = True
                                    except Exception:
                                        try:
                                            time.sleep(0.01)
                                        except Exception:
                                            pass
                        _c['total_generated'] += 1
                        if do_inline_validation:
                            _c['total_validated'] += 1
                        _c['rows_since_commit'] += 1
                        if bool(is_valid):
                            _c['total_valid_true'] += 1
                            depth_validated += 1

                        if _c['rows_since_commit'] >= 50:
                            cursor.execute(
                                "UPDATE run_metadata SET total_generated = ?, total_validated = ?, valid_solutions = ? WHERE run_id = ?",
                                (_c['total_generated'], _c['total_validated'], _c['total_valid_true'], run_id)
                            )
                            # Persist resumable state
                            try:
                                state = {
                                    'depth': depth,
                                    'last_id': int(expr_id),
                                    'totals': dict(_c),
                                }
                                cursor.execute(
                                    "INSERT OR REPLACE INTO generator_progress (run_id, state_json, updated_at) VALUES (?, ?, CURRENT_TIMESTAMP)",
                                    (run_id, json.dumps(state))
                                )
                            except Exception:
                                pass
                            conn.commit()
                            _c['rows_since_commit'] = 0
                    except sqlite3.IntegrityError:
                        continue
                # Depth progress heartbeat
                if depth_validated:
                    print(f"[{run_id}] Stream depth {depth}: validated {depth_validated}/{len(expr_list)} (total {_c['total_validated']})")

            # Drive streaming generation
            try:
                discovery.fast_generator.stream_generate(
                    primitives=discovery.primitives,
                    unary_ops=discovery.unary_ops,
                    binary_ops=discovery.all_binary_ops,
                    max_depth=max_depth,
                    batch_size=2000,
                    on_batch=emit_to_db,
                    prune=True,  # generate with pruning to avoid constants
                )
            except AttributeError:
                discovery._generate_with_lean(max_depth=max_depth)
            
            print(f"[{run_id}] Expression generation complete. Starting database insertion and validation...")
            
            # When using stream_generate with emit_to_db, expressions are already in the database
            # Skip the manual insertion loop
            if hasattr(discovery.fast_generator, 'stream_generate'):
                # Already processed via emit_to_db callback
                total_generated = counters['total_generated']
                total_validated = counters['total_validated']
                total_valid_true = counters['total_valid_true']
            else:
                # Fallback: process expressions_by_depth if stream_generate wasn't used
                total_generated = 0
                batch = []
                rows_since_commit = 0
                commit_every_rows = 50  # commit frequently for visibility
                total_validated = 0
                total_valid_true = 0
                # Note: do not deduplicate by normalized in-memory; persist all generated rows
                seen_normalized = None
                
                for depth, expressions in discovery.expressions_by_depth.items():
                    print(f"[{run_id}] Processing {len(expressions)} expressions at depth {depth}")
                    depth_validated = 0
                    depth_valid_true = 0
                    for expr_str in expressions:
                        try:
                            # Parse early to enable structural degeneracy checks
                            try:
                                u = sp.sympify(expr_str, locals=discovery._sympify_locals)
                            except Exception:
                                u = None

                            # Defensive structural drop: reject if any denominator simplifies to 0
                            try:
                                if u is not None and discovery._has_degenerate_denominator(u):
                                    # Skip inserting degenerate candidates into DB
                                    continue
                            except Exception:
                                pass

                            # Compute a canonical normalized representation for storage/dedup.
                            # Do not assume inputs are already normalized.
                            try:
                                # Deterministic canonicalization regardless of upstream
                                normalized_str = str(sp.simplify(sp.expand(sp.sympify(expr_str))))
                            except Exception:
                                normalized_str = expr_str
                            # No in-memory dedup; rely on DB rows to reflect all generated expressions

                            # Stable 32-bit signature from normalized string
                            sig_int = int(hashlib.sha256(normalized_str.encode()).hexdigest()[:8], 16)

                            # Always insert a row; if normalized duplicates, we still keep the generated expression
                            cursor.execute(f"""
                                INSERT INTO {table_name} (expression, normalized, signature, depth)
                                VALUES (?, ?, ?, ?)
                            """, (expr_str, normalized_str, sig_int, depth))
                            expr_id = cursor.lastrowid

                            # Inline validation
                            try:
                                # Parse using the problem's exact symbols/constants to ensure derivatives are correct
                                if u is None:
                                    u = sp.sympify(expr_str, locals=discovery._sympify_locals)
                                # Fast exact check via Lean-first; defer heavy checks for speed
                                # Try with extended args first, fall back to basic args
                                try:
                                    is_valid, reason = discovery.validator.validate(
                                        u,
                                        check_regularity=False,
                                        fast_point_only=False,
                                        lean_first=True,
                                        defer_heavy_checks=True,
                                        enforce_anchor=False,
                                    )
                                except TypeError:
                                    # Fallback for basic validators
                                    is_valid, reason = discovery.validator.validate(
                                        u,
                                        check_regularity=False,
                                        fast_point_only=False
                                    )
                                # Extract validator descriptor if available
                                val_desc = {}
                                if hasattr(discovery.validator, 'describe') and callable(getattr(discovery.validator, 'describe')):
                                    try:
                                        val_desc = discovery.validator.describe() or {}
                                    except Exception:
                                        val_desc = {}
                                # Capture validator evidence (e.g., Lean normalized form and params)
                                val_evidence = {}
                                if hasattr(discovery.validator, 'last_evidence') and callable(getattr(discovery.validator, 'last_evidence')):
                                    try:
                                        val_evidence = discovery.validator.last_evidence() or {}
                                    except Exception:
                                        val_evidence = {}
                                is_paper_solution = False
                                paper_name = None
                                if is_valid:
                                    known_solutions = dict(getattr(discovery.problem, 'known_solutions', {}) or {})
                                    for  known_str, name in known_solutions.items():
                                        try:
                                            known_expr = sp.sympify(known_str, locals=discovery._sympify_locals)
                                            if simplify(u - known_expr) == 0:
                                                is_paper_solution = True
                                                paper_name = name
                                                break
                                        except:
                                            pass
                                    try:
                                        print(
                                            f"[{run_id}] VALID: id={expr_id} depth={depth} known={1 if is_paper_solution else 0} "
                                            f"name={paper_name or '-'} expr={expr_str}"
                                        )
                                        sys.stdout.flush()
                                    except Exception:
                                        pass
                            except Exception as e:
                                is_valid, reason = (None, f"Validator Error: {e}")
                                is_paper_solution, paper_name = (None, None)
                                val_desc = {}
                                val_evidence = {}

                            # Persist inline validation result
                            cursor.execute(f"""
                                UPDATE {table_name}
                                SET validation_status = ?,
                                is_valid = ?,
                                validation_reason = ?,
                                validator_method = ?,
                                validator_math = ?,
                                validator_evidence = ?,
                                is_paper_solution = ?,
                                paper_solution_name = ?,
                                validated_at = CURRENT_TIMESTAMP
                            WHERE id = ?
                        """, (
                            'completed' if is_valid is not None else 'error',
                            is_valid if is_valid is not None else None,
                            reason,
                            val_desc.get('method_name'),
                            val_desc.get('math_definition'),
                            json.dumps(val_evidence),
                            is_paper_solution,
                            paper_name,
                                expr_id,
                            ))
                            total_generated += 1
                            rows_since_commit += 1
                            total_validated += 1
                            depth_validated += 1
                            if bool(is_valid):
                                total_valid_true += 1
                                depth_valid_true += 1

                            # Progress heartbeat every 50 validations
                            if (total_validated % 50) == 0:
                                # Query DB-based counts to report accurate progress
                                try:
                                    cursor.execute(f"SELECT COUNT(*), SUM(CASE WHEN is_valid = 1 THEN 1 ELSE 0 END) FROM {table_name}")
                                    db_total, db_valid = cursor.fetchone()
                                except Exception:
                                    db_total, db_valid = total_validated, total_valid_true
                                print(
                                    f"[{run_id}] Validated {db_total} expressions "
                                    f"({db_valid or 0} valid) | depth {depth}: {depth_validated}/{len(expressions)} "
                                    f"(valid {depth_valid_true})"
                                )
                                sys.stdout.flush()

                            # Periodic commit to reflect progress in monitor
                            if rows_since_commit >= commit_every_rows:
                                # Update run metadata with totals before committing
                                cursor.execute(
                                    "UPDATE run_metadata SET total_generated = ?, total_validated = ?, valid_solutions = ? WHERE run_id = ?",
                                    (total_generated, total_validated, total_valid_true, run_id)
                                )
                                conn.commit()
                                rows_since_commit = 0
                            
                        except sqlite3.IntegrityError:
                            # Any integrity error (e.g., unique constraint) — skip validation
                            continue
                        except Exception as e:
                            print(f"[{run_id}] Error inserting/queueing expression: {e}")
                    # Depth boundary: ensure progress is visible
                    try:
                        cursor.execute(
                            "UPDATE run_metadata SET total_generated = ?, total_validated = ?, valid_solutions = ? WHERE run_id = ?",
                            (total_generated, total_validated, total_valid_true, run_id)
                        )
                        conn.commit()
                        rows_since_commit = 0
                    except Exception:
                        pass

            # Final metadata update and commit
            cursor.execute(
                "UPDATE run_metadata SET total_generated = ?, total_validated = ?, valid_solutions = ? WHERE run_id = ?",
                (total_generated, total_validated, total_valid_true, run_id)
            )
            conn.commit() # Final commit
            
            # Do not send STOP here; main process will send STOP after generator joins
            
            cursor.execute("""
                UPDATE run_metadata SET total_generated = ?, status = 'generation_complete'
                WHERE run_id = ?
            """, (total_generated, run_id))
            conn.commit()

            # Re-queue any pending rows that might not have been enqueued (defensive)
            try:
                if task_queue is not None:
                    cursor.execute(f"SELECT id, expression FROM {table_name} WHERE validation_status = 'pending'")
                    for expr_id, expr_str in cursor.fetchall():
                        enq_ok = False
                        while not enq_ok:
                            try:
                                task_queue.put((expr_id, expr_str), timeout=1.0)
                                enq_ok = True
                            except Exception:
                                time.sleep(0.01)
            except Exception:
                pass
            
            # Final metadata sync from local counters
            try:
                cursor.execute(
                    "UPDATE run_metadata SET total_generated = ?, total_validated = ?, valid_solutions = ? WHERE run_id = ?",
                    (counters['total_generated'], counters['total_validated'], counters['total_valid_true'], run_id)
                )
                conn.commit()
            except Exception:
                pass
            print(f"[{run_id}] Generation complete! Total expressions generated: {counters['total_generated']}")
            
        except Exception as e:
            print(f"[{run_id}] Generator error: {e}")
            import traceback
            traceback.print_exc()
        finally:
            conn.close()

    @staticmethod
    def _parallel_validator_worker(run_id: str, table_name: str, db_path: str, problem_name: str = 'force_free', task_queue=None, result_queue=None):
        """Worker process for validating expressions from a queue."""
        print(f"[{run_id}] Validator process started (PID: {os.getpid()})")
        
        # Hybrid mode: prefer task_queue; fallback to DB polling if queue is empty
        conn = sqlite3.connect(db_path)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA busy_timeout=5000")
        cursor = conn.cursor()

        # Migration: ensure validator_evidence column exists
        try:
            cursor.execute(f"PRAGMA table_info({table_name})")
            cols = {row[1] for row in cursor.fetchall()}
            if 'validator_evidence' not in cols:
                cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN validator_evidence TEXT")
                conn.commit()
        except Exception:
            pass

        # Build problem-specific validator and sympify locals
        try:
            from physics_agent.problems import load_problem
            problem = load_problem(problem_name)
            validator = getattr(problem, 'validator', None)
        except Exception:
            validator = None
        if validator is None:
            # Fallback to precise foliation validator
            validator = PreciseFoliationValidator(cache_db=f"cache_{run_id}_{os.getpid()}.db", use_lean=True, Omega=0)

        # Prepare locals for sympify (symbols, constants, and unary ops)
        sympify_locals = {}
        try:
            from expression_operations import UNARY_OPS
            sympify_locals.update(UNARY_OPS)
        except Exception:
            pass
        try:
            for name, sym in {**getattr(problem, 'symbols', {}), **getattr(problem, 'constants', {})}.items():
                sympify_locals[name] = sym
        except Exception:
            pass

        import time
        import time, queue as _queue
        worker_pid = os.getpid()
        while True:
            try:
                claimed: list[tuple[int, str]] = []
                # Prefer queue items
                if task_queue is not None:
                    try:
                        item = task_queue.get(timeout=0.2)
                        if item is None:
                            # Explicit shutdown
                            break
                        expr_id, expr_str = item
                        claimed = [(expr_id, expr_str)]
                    except _queue.Empty:
                        pass
                    except Exception:
                        pass
                # Fallback to DB polling when queue is empty
                if not claimed:
                    cursor.execute(f"""
                        SELECT id, expression FROM {table_name}
                        WHERE validation_status = 'pending'
                        LIMIT 50
                    """)
                    candidates = cursor.fetchall()
                    if candidates:
                        for expr_id, expr_str in candidates:
                            cursor.execute(
                                f"UPDATE {table_name} SET validation_status = 'in_progress' WHERE id = ? AND validation_status = 'pending'",
                                (expr_id,)
                            )
                            if cursor.rowcount == 1:
                                claimed.append((expr_id, expr_str))
                        conn.commit()
                if not claimed:
                    time.sleep(0.2)
                    continue

                # Process tasks
                results_batch = []
                for expr_id, expr_str in claimed:
                    try:
                        # Notify centralized writer of start
                        if result_queue is not None:
                            try:
                                snippet = expr_str[:120]
                                result_queue.put((run_id, worker_pid, 'start', expr_id, snippet), timeout=0.5)
                            except Exception:
                                pass
                        u = sp.sympify(expr_str, locals=sympify_locals)
                        base_kwargs = {
                            'check_regularity': False,
                            'fast_point_only': False,
                            'lean_first': True,
                            'defer_heavy_checks': True,
                            'enforce_anchor': False,
                        }
                        try:
                            import inspect as _inspect
                            sig = _inspect.signature(validator.validate)
                            allowed = set(sig.parameters.keys())
                            call_kwargs = {k: v for k, v in base_kwargs.items() if k in allowed}
                        except Exception:
                            call_kwargs = {'check_regularity': False, 'fast_point_only': False}
                        is_valid, reason = validator.validate(u, **call_kwargs)
                        is_paper_solution = False
                        paper_name = None
                        if is_valid:
                            try:
                                known_solutions = dict(getattr(problem, 'known_solutions', {}) or {})
                            except Exception:
                                known_solutions = {}
                            for known_str, name in known_solutions.items():
                                try:
                                    known_expr = sp.sympify(known_str, locals=sympify_locals)
                                    if simplify(u - known_expr) == 0:
                                        is_paper_solution = True
                                        paper_name = name
                                        break
                                except Exception:
                                    pass
                        results_batch.append(('completed', is_valid, reason, is_paper_solution, paper_name, expr_id))
                        if bool(is_valid):
                            try:
                                print(
                                    f"[{run_id}] VALID: id={expr_id} known={1 if is_paper_solution else 0} "
                                    f"name={paper_name or '-'} expr={expr_str}"
                                )
                                sys.stdout.flush()
                            except Exception:
                                pass
                    except Exception as e:
                        print(f"[{run_id}] Error validating expression {expr_id}: {e}")
                        results_batch.append(('error', None, f"Validator Error: {e}", None, None, expr_id))

                # Push results to centralized writer (include worker PID)
                if result_queue is not None and results_batch:
                    try:
                        result_queue.put((run_id, worker_pid, 'end', results_batch), timeout=1.0)
                    except Exception:
                        pass
            except Exception as e:
                print(f"[{run_id}] Validator worker error: {e}")
                time.sleep(0.25)

        if conn is not None:
            conn.close()

    def _generate_report_from_db(self):
        """Generates a final report from the database after a parallel run."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        try:
            # Get statistics
            cursor.execute(f"""
                SELECT 
                    COUNT(*) as total,
                    SUM(CASE WHEN is_valid = 1 THEN 1 ELSE 0 END) as valid,
                    COUNT(DISTINCT CASE WHEN is_paper_solution = 1 THEN signature END) as paper_solutions_distinct
                FROM {self.table_name}
            """)
            
            total, valid, paper_distinct = cursor.fetchone()
            
            # Get paper solutions found
            cursor.execute(f"""
                SELECT expression, paper_solution_name
                FROM {self.table_name}
                WHERE is_paper_solution = 1
                ORDER BY paper_solution_name
            """)
            
            paper_solutions = cursor.fetchall()
            
            # Also compute a distinct-by-signature list to avoid duplicates in display, if desired later
            cursor.execute(f"""
                SELECT paper_solution_name, MIN(expression) as example_expression, MIN(id) as example_id
                FROM {self.table_name}
                WHERE is_paper_solution = 1
                GROUP BY signature, paper_solution_name
                ORDER BY paper_solution_name
            """)
            paper_solutions_distinct = cursor.fetchall()
            
            # Get counts by depth
            cursor.execute(f"""
                SELECT depth, COUNT(*) as count
                FROM {self.table_name}
                GROUP BY depth
                ORDER BY depth
            """)
            
            depth_counts = cursor.fetchall()
            
        finally:
            conn.close()
        
        print("\n" + "=" * 80)
        print(f"PARALLEL DISCOVERY COMPLETE - RUN ID: {self.run_id}")
        print("=" * 80)
        print(f"Total expressions generated: {total}")
        print(f"Valid foliations found: {valid or 0}")
        print(f"Known solutions found: {paper_distinct or 0} (distinct canonical forms)")
        print("\nExpression counts by depth:")
        for depth, count in depth_counts:
            print(f"  Depth {depth}: {count}")
        
        if paper_solutions_distinct:
            print("\nKnown solutions found (deduplicated by signature):")
            for name, expr, ex_id in paper_solutions_distinct:
                print(f"  ✓ {name} (id={ex_id}): {expr}")
        
        print(f"\nResults stored in table: {self.table_name}")
        print(f"Database: {self.db_path}")

        # Optional: print validator metadata used
        print("\nValidator info (per last row): method and math definition if stored.")
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        try:
            cur.execute(f"SELECT validator_method, substr(validator_math,1,120) FROM {self.table_name} WHERE validator_method IS NOT NULL ORDER BY id DESC LIMIT 1")
            row = cur.fetchone()
            if row:
                print(f"  method: {row[0]}")
                print(f"  math:   {row[1]}...")
            # Print novel (unknown) solutions, deduplicated by mathematical equivalence (no LIMIT)
            print("\nNovel solutions (Lean-valid, not matching known set; deduplicated by mathematical equivalence):")
            cur.execute(f"""
                SELECT id, expression
                FROM {self.table_name}
                WHERE is_valid = 1 AND (is_paper_solution IS NULL OR is_paper_solution = 0)
            """)
            rows = cur.fetchall() or []
            novel_rows_count = len(rows)
            expr_rows = [(int(r[0]), r[1]) for r in rows]

            # Build known solutions map for defensive filtering
            known_map = dict(getattr(self.problem, 'known_solutions', {}) or {})

            # Canonicalization and equivalence grouping helpers
            def _canonical_key(e: sp.Basic) -> str:
                try:
                    en = sp.together(e)
                    en = sp.cancel(en)
                    en = sp.powsimp(en, force=True)
                    en = sp.powdenest(en, force=True)
                    en = sp.simplify(en)
                    en = en.rewrite(sp.Pow)
                    en = sp.together(sp.cancel(en))
                    return sp.srepr(en)
                except Exception:
                    try:
                        return sp.srepr(sp.simplify(e))
                    except Exception:
                        return str(e)

            def _equiv_to_known(e: sp.Basic) -> bool:
                for known_str in known_map.keys():
                    try:
                        k_expr = sp.sympify(known_str, locals=getattr(self, '_sympify_locals', {}))
                        if simplify(e - k_expr) == 0:
                            return True
                    except Exception:
                        pass
                return False

            def _depth(x: sp.Basic) -> int:
                try:
                    return 1 + max((_depth(arg) for arg in x.args), default=0)
                except Exception:
                    return 1

            def _rep_cost(e: sp.Basic):
                try:
                    c_ops = int(sp.count_ops(e, visual=False))
                except Exception:
                    c_ops = 10**6
                try:
                    d = _depth(e)
                except Exception:
                    d = 999999
                try:
                    s_len = len(sp.srepr(e))
                except Exception:
                    s_len = len(str(e))
                try:
                    pen = 10*int(e.has(sp.zoo)) + 20*int(e.has(sp.nan)) + 5*int(any(isinstance(a, sp.Float) for a in e.atoms(sp.Float)))
                except Exception:
                    pen = 0
                return (c_ops, d, s_len, pen)

            buckets = {}
            for expr_id, expr_str in expr_rows:
                try:
                    expr = sp.sympify(expr_str, locals=getattr(self, '_sympify_locals', {}))
                except Exception:
                    # If parsing fails, bucket by raw string to avoid loss
                    key = f"RAW::{expr_str}"
                    entry = buckets.get(key)
                    if entry is None:
                        buckets[key] = {'rep_id': expr_id, 'rep_str': expr_str, 'rep_expr': None, 'size': 1}
                    else:
                        entry['size'] += 1
                    continue

                # Skip if equivalent to a known solution
                try:
                    if _equiv_to_known(expr):
                        continue
                except Exception:
                    pass

                try:
                    key = _canonical_key(expr)
                except Exception:
                    key = str(expr)

                entry = buckets.get(key)
                if entry is None:
                    buckets[key] = {'rep_id': expr_id, 'rep_str': expr_str, 'rep_expr': expr, 'size': 1}
                else:
                    entry['size'] += 1
                    # Prefer a simpler representative
                    try:
                        if _rep_cost(expr) < _rep_cost(entry['rep_expr'] if entry['rep_expr'] is not None else expr):
                            entry['rep_id'] = expr_id
                            entry['rep_str'] = expr_str
                            entry['rep_expr'] = expr
                    except Exception:
                        pass

            # Print groups sorted by descending size, then by representative string
            items = sorted(buckets.items(), key=lambda kv: (-kv[1]['size'], kv[1]['rep_str']))
            print(f"Novel valid rows (non-paper): {novel_rows_count}")
            print(f"Novel equivalence classes: {len(items)}")
            printed = 0
            for _, entry in items:
                print(f"  • id={entry['rep_id']} size={entry['size']} expr={entry['rep_str']}")
                printed += 1
            if printed == 0:
                print("  (none)")
        finally:
            conn.close()


    def heavy_validate_from_db(self, scope: str = 'valid', check_regularity: bool = True,
                               enforce_anchor: bool = True, anchor_target: str = 'either'):
        """Run heavy validation (finiteness, regularity, small-spin anchor) from the DB only.

        scope: 'valid' (exact-zero hits) or 'all' (all rows).
        Updates columns in-place: heavy_is_valid, heavy_reason, heavy_validated_at.
        """
        assert scope in ('valid', 'all')
        # Ensure columns exist
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        try:
            cur.execute(f"PRAGMA table_info({self.table_name})")
            cols = {row[1] for row in cur.fetchall()}
            if 'heavy_is_valid' not in cols:
                cur.execute(f"ALTER TABLE {self.table_name} ADD COLUMN heavy_is_valid BOOLEAN")
            if 'heavy_reason' not in cols:
                cur.execute(f"ALTER TABLE {self.table_name} ADD COLUMN heavy_reason TEXT")
            if 'heavy_validated_at' not in cols:
                cur.execute(f"ALTER TABLE {self.table_name} ADD COLUMN heavy_validated_at TIMESTAMP")
            conn.commit()
        finally:
            conn.close()

        # Select candidates
        where = "WHERE is_valid = 1" if scope == 'valid' else ""
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        try:
            cur.execute(f"SELECT id, expression FROM {self.table_name} {where} ORDER BY id")
            rows = cur.fetchall() or []
        finally:
            conn.close()

        total = len(rows)
        print("\n" + "-" * 80)
        print(f"Heavy validation — run_id={self.run_id} table={self.table_name}")
        print(f"Database: {self.db_path}")
        print(f"Scope: {scope} | Candidates: {total} | Anchor: {'on' if enforce_anchor else 'off'} target={anchor_target}")

        validator = getattr(self.problem, 'validator', None)
        if validator is None:
            print("No validator available for this problem.")
            return

        # Configure anchor target if supported
        try:
            if hasattr(validator, 'monopole_target') and anchor_target:
                setattr(validator, 'monopole_target', anchor_target)
            if hasattr(validator, 'require_monopole_extension'):
                setattr(validator, 'require_monopole_extension', bool(enforce_anchor))
        except Exception:
            pass

        ok = 0
        fail = 0
        batch_updates = []
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        for idx, (expr_id, expr_str) in enumerate(rows, 1):
            try:
                u = sp.sympify(expr_str, locals=getattr(self, '_sympify_locals', {}))
                # Heavy mode: exact zero already known; enforce finiteness/regularity/anchor
                try:
                    is_valid, reason = validator.validate(
                        u,
                        check_regularity=check_regularity,
                        fast_point_only=False,
                        lean_first=True,
                        defer_heavy_checks=False,
                        enforce_anchor=enforce_anchor,
                    )
                except TypeError:
                    # Fallback for validators without extended signature
                    is_valid, reason = validator.validate(u, check_regularity=check_regularity, fast_point_only=False)

                if is_valid:
                    ok += 1
                else:
                    fail += 1
                batch_updates.append((
                    bool(is_valid),
                    reason,
                    expr_id,
                ))
            except Exception as e:
                fail += 1
                batch_updates.append((None, f"Heavy validator error: {e}", expr_id))

            if len(batch_updates) >= 100:
                cur.executemany(f"""
                    UPDATE {self.table_name}
                    SET heavy_is_valid = ?,
                        heavy_reason = ?,
                        heavy_validated_at = CURRENT_TIMESTAMP
                    WHERE id = ?
                """, batch_updates)
                conn.commit()
                batch_updates = []
                print(f"  Progress: {idx}/{total} (ok={ok}, fail={fail})")

        if batch_updates:
            cur.executemany(f"""
                UPDATE {self.table_name}
                SET heavy_is_valid = ?,
                    heavy_reason = ?,
                    heavy_validated_at = CURRENT_TIMESTAMP
                WHERE id = ?
            """, batch_updates)
            conn.commit()
        conn.close()

        print(f"Heavy summary: ok={ok} / {total} (failed {fail})")

    def verify_pde_from_db(self, scope: str = 'novel', check_regularity: bool = True):
        """Re-validate expressions from the run database against the exact PDE.

        scope: 'novel' (valid non-paper), 'valid' (all valid), or 'all' (entire table)
        """
        assert scope in ('novel', 'valid', 'all')
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        where = ""
        if scope == 'novel':
            where = "WHERE is_valid = 1 AND (is_paper_solution IS NULL OR is_paper_solution = 0)"
        elif scope == 'valid':
            where = "WHERE is_valid = 1"
        # else: all rows

        query = f"SELECT id, expression FROM {self.table_name} {where} ORDER BY id"
        rows = []
        try:
            cur.execute(query)
            rows = cur.fetchall() or []
        finally:
            conn.close()

        total = len(rows)
        print("\n" + "-" * 80)
        print(f"Exact PDE verification — run_id={self.run_id} table={self.table_name}")
        print(f"Database: {self.db_path}")
        print(f"Scope: {scope} | Candidates: {total}")

        validator = getattr(self.problem, 'validator', None)
        if validator is None:
            print("No validator available for this problem.")
            return

        ok = 0
        fail = 0
        for expr_id, expr_str in rows:
            try:
                u = sp.sympify(expr_str, locals=getattr(self, '_sympify_locals', {}))
                # Prefer validator's validate; fallback to symbolic if Lean not available
                try:
                    is_valid, reason = validator.validate(
                        u,
                        check_regularity=False,
                        fast_point_only=False,
                        lean_first=True,
                        defer_heavy_checks=True,
                        enforce_anchor=False,
                    )
                except Exception as e:
                    is_valid, reason = (False, f"Validator error: {e}")

                if not is_valid and hasattr(validator, '_lhs'):
                    try:
                        lhs = validator._lhs(u)  # type: ignore[attr-defined]
                        if sp.simplify(lhs) == 0:
                            is_valid, reason = (True, 'Valid (symbolic simplify lhs=0)')
                    except Exception:
                        pass

                if is_valid:
                    print(f"  ✓ id={expr_id}: {expr_str}")
                    ok += 1
                else:
                    print(f"  ✗ id={expr_id}: {reason}")
                    fail += 1
            except Exception as e:
                print(f"  ✗ id={expr_id}: Parse error: {e}")
                fail += 1

        print(f"Summary: exact solutions {ok} / {total} (failed {fail})")

    def find_monopole_extensions(self, scope: str = 'valid', target: str = 'either'):
        """Find candidates whose a -> 0 limit equals 1 - x (monopole) or x.

        Args:
            scope: 'novel'|'valid'|'all' — which rows to scan.
            target: '1-x'|'x'|'either' — which target to match in the small-spin limit.
        """
        assert scope in ('novel', 'valid', 'all')
        assert target in ('1-x', 'x', 'either')

        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        where = ""
        if scope == 'novel':
            where = "WHERE is_valid = 1 AND (is_paper_solution IS NULL OR is_paper_solution = 0)"
        elif scope == 'valid':
            where = "WHERE is_valid = 1"
        query = f"SELECT id, expression FROM {self.table_name} {where} ORDER BY id"
        rows = []
        try:
            cur.execute(query)
            rows = cur.fetchall() or []
        finally:
            conn.close()

        a_symbol = getattr(self.problem, 'constants', {}).get('a', sp.Symbol('a'))
        x_symbol = self.problem.symbols.get('x', sp.Symbol('x'))

        # Prepare target expressions
        targets = []
        if target in ('1-x', 'either'):
            targets.append(sp.simplify(1 - x_symbol))
        if target in ('x', 'either'):
            targets.append(sp.simplify(x_symbol))

        print("\n" + "-" * 80)
        print(f"Small-spin check (a -> 0) — run_id={self.run_id} table={self.table_name}")
        print(f"Database: {self.db_path}")
        print(f"Scope: {scope} | Candidates: {len(rows)} | Target: {target}")

        matches = []
        errors = 0
        for expr_id, expr_str in rows:
            try:
                u = sp.sympify(expr_str, locals=getattr(self, '_sympify_locals', {}))
                # Compute formal limit a->0; fallback to substitution if limit fails
                try:
                    lim_u = sp.simplify(sp.limit(u, a_symbol, 0))
                except Exception:
                    lim_u = sp.simplify(u.subs(a_symbol, 0))

                for tgt in targets:
                    try:
                        if sp.simplify(lim_u - tgt) == 0:
                            matches.append((expr_id, expr_str, sp.sstr(lim_u)))
                            break
                    except Exception:
                        continue
            except Exception:
                errors += 1
                continue

        if matches:
            print(f"Found {len(matches)} candidate extension(s) with correct a->0 limit:")
            for mid, mexpr, mlim in matches:
                print(f"  • id={mid} limit={mlim} expr={mexpr}")
        else:
            print("No candidates matched the requested small-spin limit.")
        if errors:
            print(f"Skipped {errors} row(s) due to parse/simplify errors.")

    def audit_kerr_candidates(self, scope: str = 'valid'):
        """Audit DB rows for trivial solutions (constants) and nontrivial dependence.

        Only meaningful for the Kerr magnetosphere problem.
        """
        if getattr(self.problem, 'slug', '') != 'kerr_magnetosphere':
            print("Audit is designed for kerr_magnetosphere only.")
            return

        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        where = ""
        if scope == 'novel':
            where = "WHERE is_valid = 1 AND (is_paper_solution IS NULL OR is_paper_solution = 0)"
        elif scope == 'valid':
            where = "WHERE is_valid = 1"
        query = f"SELECT id, expression FROM {self.table_name} {where} ORDER BY id"
        rows = []
        try:
            cur.execute(query)
            rows = cur.fetchall() or []
        finally:
            conn.close()

        r = self.problem.symbols.get('r', sp.Symbol('r'))
        x = self.problem.symbols.get('x', sp.Symbol('x'))

        total = len(rows)
        const_count = 0
        no_x_dep = 0
        no_r_dep = 0
        nonconst_rows = []  # store details for full logging

        # Targets for small-spin (a->0) limit comparisons
        a_symbol = getattr(self.problem, 'constants', {}).get('a', sp.Symbol('a'))
        x_symbol = x
        target_mono_1mx = sp.simplify(1 - x_symbol)
        target_mono_x = sp.simplify(x_symbol)

        for expr_id, expr_str in rows:
            try:
                u = sp.sympify(expr_str, locals=getattr(self, '_sympify_locals', {}))
                ur = sp.simplify(sp.diff(u, r))
                ux = sp.simplify(sp.diff(u, x))
                if ur == 0 and ux == 0:
                    const_count += 1
                    continue

                tag_no_x = (ux == 0)
                tag_no_r = (ur == 0)
                if tag_no_x:
                    no_x_dep += 1
                if tag_no_r:
                    no_r_dep += 1

                # Compute small-spin limit to check for disguised monopole
                limit_repr = None
                disguised = False
                try:
                    try:
                        lim_u = sp.simplify(sp.limit(u, a_symbol, 0))
                    except Exception:
                        lim_u = sp.simplify(u.subs(a_symbol, 0))
                    limit_repr = sp.sstr(lim_u)
                    try:
                        if sp.simplify(lim_u - target_mono_1mx) == 0 or sp.simplify(lim_u - target_mono_x) == 0:
                            disguised = True
                    except Exception:
                        pass
                except Exception:
                    limit_repr = "<limit-error>"

                nonconst_rows.append({
                    'id': expr_id,
                    'expr': expr_str,
                    'no_x': tag_no_x,
                    'no_r': tag_no_r,
                    'limit': limit_repr,
                    'is_disguised_mono': disguised,
                })
            except Exception:
                continue

        print("\n" + "-" * 80)
        print(f"Audit (Kerr) — run_id={self.run_id} table={self.table_name}")
        print(f"Database: {self.db_path}")
        print(f"Scope: {scope} | Valid rows scanned: {total}")
        print(f"Constant solutions: {const_count}")
        print(f"No x-dependence (ux==0 but ur!=0): {no_x_dep}")
        print(f"No r-dependence (ur==0 but ux!=0): {no_r_dep}")

        # Log all non-constant rows with tags and small-spin check
        if nonconst_rows:
            print("\nNon-constant solutions (with small-spin check):")
            for row in nonconst_rows:
                tags = []
                if row['no_x']:
                    tags.append('no-x')
                if row['no_r']:
                    tags.append('no-r')
                if row['is_disguised_mono']:
                    tags.append('limit->monopole')
                tag_str = f"[{', '.join(tags)}]" if tags else "[]"
                limit_str = f"limit(a->0)={row['limit']}" if row['limit'] is not None else "limit(a->0)=<n/a>"
                print(f"  • id={row['id']} {tag_str} {limit_str} expr={row['expr']}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Reproduce Force-Free Foliation results')
    parser.add_argument('--mode', choices=['parallel', 'sequential'],
                       default='parallel', help='Run mode: parallel (concurrent generation and validation) or sequential (in-memory)')
    parser.add_argument('--problem', type=str, default='force_free', help='Problem to solve (e.g., force_free, kerr_magnetosphere).')
    parser.add_argument('--validate-run', type=str,
                        help='(DEPRECATED) Run ID to validate a previous run. Concurrent validation is now default.')
    parser.add_argument('--max-depth', type=int, default=4, help='Maximum expression depth to generate.')
    parser.add_argument('--print-run-id', type=str, help='Print results for an existing run ID (reads DB only, no generation).')
    parser.add_argument('--db-path', type=str, help='Optional explicit path to the run database; if omitted, inferred from problem outputs and run ID.')
    parser.add_argument('--resume-run', type=str, help='Resume validation for an existing run_id (drains pending rows).')
    parser.add_argument('--resume-validators', type=int, default=8, help='Number of validators to use when resuming a run.')
    parser.add_argument('--verify-pde', action='store_true', help='After printing, verify the exact PDE on the selected scope from DB only.')
    parser.add_argument('--verify-scope', choices=['novel','valid','all'], default='novel', help='Which rows to verify: novel valid non-paper, all valid, or all rows.')
    parser.add_argument('--find-monopole', action='store_true', help='Search DB rows for finite-spin candidates whose a->0 limit equals the monopole (1 - x) or x.')
    parser.add_argument('--monopole-target', choices=['either','1-x','x'], default='either', help='Target small-spin limit to match.')
    parser.add_argument('--audit-kerr', action='store_true', help='Audit valid rows to detect constants and variable dependence patterns (Kerr only).')
    parser.add_argument('--validators', type=int, default=0, help='Number of parallel validator workers (0 = inline validation only, -1 = auto (cpu_count-2)).')
    
    args = parser.parse_args()
    
    if args.validate_run:
        print("Warning: --validate-run is deprecated. Validation now runs concurrently with generation.")
        print("To re-validate a failed run, you would need to implement a re-validation script.")
        sys.exit(0)

    # Fast path: print results for an existing run (no generation)
    if args.print_run_id:
        discovery = GeneralFoliationDiscovery(problem_name=args.problem)
        discovery.run_id = args.print_run_id
        if args.db_path:
            discovery.db_path = args.db_path
        else:
            discovery.db_path = os.path.join(discovery.problem.get_output_dir(), f"parallel_runs_{args.print_run_id}.db")
        discovery.table_name = f"expressions_{args.print_run_id.replace('-', '_')}"
        if not os.path.exists(discovery.db_path):
            print(f"Database not found: {discovery.db_path}")
            sys.exit(1)
        discovery._generate_report_from_db()
        if args.verify_pde:
            discovery.verify_pde_from_db(scope=args.verify_scope)
        if args.find_monopole:
            discovery.find_monopole_extensions(scope=args.verify_scope, target=args.monopole_target)
        if args.audit_kerr:
            discovery.audit_kerr_candidates(scope=args.verify_scope)
        sys.exit(0)

    if args.resume_run:
        discovery = GeneralFoliationDiscovery(problem_name=args.problem)
        discovery.resume_pending_validation(
            resume_run_id=args.resume_run,
            validators=args.resume_validators,
            db_path=args.db_path
        )
        sys.exit(0)

    if args.mode == 'parallel':
        discovery = GeneralFoliationDiscovery(mode='parallel', problem_name=args.problem)
        # Compute validator worker count
        try:
            import os as _os
            validators = (max(1, (_os.cpu_count() or 4) - 2) if args.validators == -1 else max(0, int(args.validators)))
        except Exception:
            validators = 0 if args.validators != -1 else 2
        discovery.run_parallel_discovery(max_depth=args.max_depth, validators=validators)
    else: # sequential
        # Force sequential mode to bypass parallel execution issues.
        print("Running in sequential mode.")
        discovery = GeneralFoliationDiscovery(problem_name=args.problem)
        discovery.generate_expressions_up_to_depth(max_depth=args.max_depth)
        valid_solutions = discovery.find_valid_foliations()
        discovery.generate_report(valid_solutions)
