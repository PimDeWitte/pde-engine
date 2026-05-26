import sympy as sp
import unittest

from tools.export_kerr_paper_target_gate import (
    CORRECTION_ANGULAR_FACTORS,
    CORRECTION_COEFFICIENTS,
    CORRECTION_RADIAL_FACTORS,
    FULL_TARGET_REJECTION,
    PAIR_CORRECTION_ANGULAR_FACTORS,
    PAIR_CORRECTION_COEFFICIENTS,
    PAIR_CORRECTION_RADIAL_FACTORS,
    assess_expression,
    correction_basis,
    exact_zero,
    generate_anchor_pair_correction_rows,
    full_kerr_split_monopole_residual,
    generate_anchor_correction_rows,
    pair_correction_basis,
    sympify_locals,
)


class KerrPaperTargetGateTest(unittest.TestCase):
    def test_degenerate_expression_is_rejected(self):
        assessment = assess_expression("1/(1 - 1)")

        self.assertFalse(assessment["paper_admissible"])
        self.assertIn("non-finite symbolic expression", assessment["paper_rejection_reasons"])
        self.assertEqual(assessment["full_target"]["reason"], "skipped because strict prechecks failed")

    def test_known_small_spin_anchor_is_not_a_paper_candidate(self):
        assessment = assess_expression("1 - x")

        self.assertFalse(assessment["paper_admissible"])
        self.assertEqual(assessment["small_spin_limit"], "1 - x")
        self.assertTrue(assessment["small_spin_anchor_matches"])
        self.assertTrue(assessment["equivalent_to_known_anchor"])
        self.assertIn("missing required variables: a, r", assessment["paper_rejection_reasons"])
        self.assertIn("equivalent to known small-spin anchor", assessment["paper_rejection_reasons"])

    def test_schwarzschild_monopole_anchor_solves_a_zero_limit(self):
        locals_map = sympify_locals()
        a = locals_map["a"]
        x = locals_map["x"]

        residual = full_kerr_split_monopole_residual(1 - x)

        self.assertTrue(exact_zero(sp.factor(sp.cancel(sp.together(residual.subs(a, 0))))))

    def test_nontrivial_anchor_like_expression_fails_full_target(self):
        assessment = assess_expression("1 - x + a**2*r*x")

        self.assertFalse(assessment["paper_admissible"])
        self.assertEqual(assessment["variables"], ["a", "r", "x"])
        self.assertEqual(assessment["small_spin_limit"], "1 - x")
        self.assertTrue(assessment["small_spin_anchor_matches"])
        self.assertIn(FULL_TARGET_REJECTION, assessment["paper_rejection_reasons"])
        self.assertFalse(assessment["full_target"]["exact_zero"])

    def test_anchor_correction_grammar_preserves_small_spin_anchor(self):
        locals_map = sympify_locals()
        a = locals_map["a"]
        x = locals_map["x"]

        rows = generate_anchor_correction_rows()

        self.assertEqual(
            len(rows),
            len(CORRECTION_ANGULAR_FACTORS)
            * len(CORRECTION_RADIAL_FACTORS)
            * len(CORRECTION_COEFFICIENTS),
        )
        self.assertEqual(
            len(correction_basis()),
            len(CORRECTION_ANGULAR_FACTORS) * len(CORRECTION_RADIAL_FACTORS),
        )
        self.assertEqual(len({row.expression for row in rows}), len(rows))
        for row in rows:
            expr = sp.sympify(row.expression, locals=locals_map)
            self.assertEqual(row.source, "anchor_correction_grammar")
            self.assertIn("expanded anchor correction grammar:", row.validation_reason or "")
            self.assertEqual(sp.simplify(sp.limit(expr, a, 0) - (1 - x)), 0)
            self.assertNotEqual(sp.simplify(expr - (1 - x)), 0)
            self.assertTrue(
                {"a", "r", "x"}.issubset({str(symbol) for symbol in expr.free_symbols})
            )

    def test_pair_correction_grammar_preserves_small_spin_anchor(self):
        locals_map = sympify_locals()
        a = locals_map["a"]
        x = locals_map["x"]

        rows = generate_anchor_pair_correction_rows()
        pair_basis_count = len(PAIR_CORRECTION_ANGULAR_FACTORS) * len(
            PAIR_CORRECTION_RADIAL_FACTORS
        )
        expected_pair_count = (
            pair_basis_count
            * (pair_basis_count - 1)
            // 2
            * len(PAIR_CORRECTION_COEFFICIENTS) ** 2
        )

        self.assertEqual(len(pair_correction_basis()), pair_basis_count)
        self.assertEqual(len(rows), expected_pair_count)
        self.assertEqual(len({row.expression for row in rows}), len(rows))
        for row in rows:
            expr = sp.sympify(row.expression, locals=locals_map)
            self.assertEqual(row.source, "anchor_pair_correction_grammar")
            self.assertIn("two-term anchor correction grammar:", row.validation_reason or "")
            self.assertEqual(sp.simplify(sp.limit(expr, a, 0) - (1 - x)), 0)
            self.assertNotEqual(sp.simplify(expr - (1 - x)), 0)
            self.assertTrue(
                {"a", "r", "x"}.issubset({str(symbol) for symbol in expr.free_symbols})
            )


if __name__ == "__main__":
    unittest.main()
