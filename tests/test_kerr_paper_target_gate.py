import unittest

from tools.export_kerr_paper_target_gate import (
    TARGET_BLOCKER,
    assess_expression,
)


class KerrPaperTargetGateTest(unittest.TestCase):
    def test_degenerate_expression_is_rejected(self):
        assessment = assess_expression("1/(1 - 1)")

        self.assertFalse(assessment["paper_admissible"])
        self.assertIn("non-finite symbolic expression", assessment["paper_rejection_reasons"])
        self.assertIn(TARGET_BLOCKER, assessment["paper_rejection_reasons"])

    def test_known_small_spin_anchor_is_not_a_paper_candidate(self):
        assessment = assess_expression("1 - x")

        self.assertFalse(assessment["paper_admissible"])
        self.assertEqual(assessment["small_spin_limit"], "1 - x")
        self.assertTrue(assessment["small_spin_anchor_matches"])
        self.assertTrue(assessment["equivalent_to_known_anchor"])
        self.assertIn("missing required variables: a, r", assessment["paper_rejection_reasons"])
        self.assertIn("equivalent to known small-spin anchor", assessment["paper_rejection_reasons"])

    def test_nontrivial_anchor_like_expression_still_waits_for_full_target(self):
        assessment = assess_expression("1 - x + a**2*r*x")

        self.assertFalse(assessment["paper_admissible"])
        self.assertEqual(assessment["variables"], ["a", "r", "x"])
        self.assertEqual(assessment["small_spin_limit"], "1 - x")
        self.assertTrue(assessment["small_spin_anchor_matches"])
        self.assertIn(TARGET_BLOCKER, assessment["paper_rejection_reasons"])


if __name__ == "__main__":
    unittest.main()
