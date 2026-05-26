import unittest

import sympy as sp

from tools.export_force_free_novel_discoveries import (
    CandidateRow,
    equivalent_to_known,
    exact_zero,
    force_free_determinant,
    known_solution_map,
    select_discoveries,
    sympify_locals,
)


class ForceFreeNovelDiscoveriesTest(unittest.TestCase):
    def test_independent_determinant_separates_solution_from_control(self):
        locals_map = sympify_locals()
        rho = locals_map["rho"]
        z = locals_map["z"]

        self.assertTrue(exact_zero(force_free_determinant(rho + z)))
        self.assertEqual(
            sp.factor(sp.cancel(sp.together(force_free_determinant(rho * z)))),
            16 * rho * z,
        )

    def test_registered_known_filter_keeps_engine_novel_candidate(self):
        locals_map = sympify_locals()
        known = known_solution_map()

        self.assertTrue(equivalent_to_known(sp.sympify("rho**2", locals=locals_map), known, locals_map))
        self.assertFalse(equivalent_to_known(sp.sympify("rho + z", locals=locals_map), known, locals_map))

    def test_discovery_selection_prefers_representative_candidates(self):
        rows = [
            CandidateRow(1, "rho", 1, "Valid"),
            CandidateRow(2, "rho + z", 2, "Valid"),
            CandidateRow(3, "rho/z", 2, "Valid"),
            CandidateRow(4, "rho**2 + z**2", 1, "Valid"),
        ]

        selected = select_discoveries(rows)

        self.assertEqual([row.expression for row in selected], ["rho + z", "rho/z", "rho**2 + z**2"])


if __name__ == "__main__":
    unittest.main()
