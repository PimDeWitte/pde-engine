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
    build_artifact,
    correction_basis,
    exact_zero,
    generate_anchor_pair_correction_rows,
    full_kerr_split_monopole_residual,
    generate_anchor_correction_rows,
    literature_slow_rotation_anchor_metadata,
    pair_correction_basis,
    solve_leading_order_coefficient_rows,
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

    def test_leading_order_coefficient_solve_finds_no_one_term_solution(self):
        rows, metadata = solve_leading_order_coefficient_rows()

        self.assertEqual(rows, [])
        self.assertTrue(metadata["enabled"])
        self.assertEqual(metadata["status"], "no_leading_order_solution")
        self.assertEqual(metadata["unknown_count"], len(correction_basis()))
        self.assertEqual(metadata["equation_count"], 66)
        self.assertEqual(metadata["matrix_shape"], [66, 48])
        self.assertEqual(metadata["linsolve_result"], "EmptySet")
        self.assertEqual(metadata["generated_candidates"], 0)

    def test_artifact_records_implemented_and_missing_paper_criteria(self):
        artifact = build_artifact(
            None,
            [],
            {
                "total_generated": 0,
                "total_completed": 0,
                "total_valid_rows": 0,
                "known_solution_rows": 0,
                "rows_loaded_for_full_target_gate": 0,
            },
            include_probes=True,
            include_corrections=False,
            include_pair_corrections=False,
            include_coefficient_solve=False,
            include_literature_anchor=False,
        )
        criteria = {item["id"]: item for item in artifact["criteria_status"]}

        self.assertEqual(criteria["full_residual"]["status"], "implemented_negative")
        self.assertEqual(criteria["literature_perturbative_anchor"]["status"], "not_run")
        self.assertEqual(
            criteria["global_regularities"]["status"],
            "not_implemented_for_positive_claim",
        )
        self.assertIn("rational safe points", criteria["global_regularities"]["evidence"])
        self.assertEqual(
            criteria["equivalence_filters"]["status"],
            "partial_not_sufficient_for_positive_claim",
        )

    def test_literature_slow_rotation_anchor_passes_leading_order_residual(self):
        metadata = literature_slow_rotation_anchor_metadata()

        self.assertTrue(metadata["enabled"])
        self.assertEqual(metadata["status"], "passes_leading_order_anchor")
        self.assertTrue(metadata["leading_residual_exact_zero"])
        self.assertEqual(metadata["leading_residual_simplified"], "0")
        self.assertFalse(metadata["finite_spin_exact_candidate"])
        self.assertIn("perturbative", metadata["claim_boundary"])
        self.assertTrue(all(item["zero"] for item in metadata["point_checks"]))


if __name__ == "__main__":
    unittest.main()
