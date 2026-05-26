import sympy as sp
import unittest

from tools.export_kerr_paper_target_gate import (
    FULL_TARGET_REJECTION,
    assess_expression,
    exact_zero,
    full_kerr_split_monopole_residual,
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


if __name__ == "__main__":
    unittest.main()
