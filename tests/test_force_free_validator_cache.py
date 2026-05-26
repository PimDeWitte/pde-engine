import tempfile
import unittest
from pathlib import Path

import sympy as sp

from problems.force_free.validator import PreciseFoliationValidator


class ForceFreeValidatorCacheTests(unittest.TestCase):
    def test_cache_key_separates_rotation_context(self):
        rho = sp.Symbol("rho", real=True, positive=True)
        z = sp.Symbol("z", real=True)
        bent = rho**2 * sp.exp(-2 * z)
        cache_db = str(Path(tempfile.mkdtemp()) / "validator_cache.db")

        nonrotating = PreciseFoliationValidator(
            cache_db=cache_db,
            use_lean=False,
            Omega=0,
        )
        rotating = PreciseFoliationValidator(
            cache_db=cache_db,
            use_lean=False,
            Omega=1,
        )

        nonrotating_valid, nonrotating_reason = nonrotating.validate(
            bent,
            check_regularity=True,
            fast_point_only=True,
        )
        rotating_valid, rotating_reason = rotating.validate(
            bent,
            check_regularity=True,
            fast_point_only=True,
        )

        self.assertTrue(nonrotating_valid, nonrotating_reason)
        self.assertFalse(rotating_valid, rotating_reason)
        self.assertNotEqual(
            nonrotating._get_expr_hash(
                bent,
                check_regularity=True,
                fast_point_only=True,
            ),
            rotating._get_expr_hash(
                bent,
                check_regularity=True,
                fast_point_only=True,
            ),
        )


if __name__ == "__main__":
    unittest.main()
