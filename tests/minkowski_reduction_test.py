import unittest

from sage.quadratic_forms.all import random_quadraticform_with_conditions, QuadraticForm
from sage.rings.all import ZZ, QQ
from sage.matrix.all import identity_matrix
from ..hecke_module import _minkowski_reduction


def is_reduced(b1, b2, b3, S):

    def _norm(a):
        return a * S * a
    n2 = _norm(b2)
    n3 = _norm(b3)
    return (_norm(b1) <= n2 <= n3 and
            all(_norm(b2 + x1 * b1) >= n2 for x1 in range(-3, 3)) and
            all(_norm(b3 + x2 * b2 + x1 * b1) >= n3 for x1 in range(-3, 3) for x2 in range(-3, 3)))


class MinkowskiReductionTest(unittest.TestCase):

    def test_reduction(self):
        for _ in range(100):
            q = random_quadraticform_with_conditions(
                ZZ, 3,
                [QuadraticForm.is_positive_definite], [-20, 20])
            S = q.Gram_matrix_rational()
            print S.list()
            b1, b2, b3 = identity_matrix(QQ, 3).columns()
            b1, b2, b3 = _minkowski_reduction(b1, b2, b3, S)
            self.assertTrue(is_reduced(b1, b2, b3, S))

suite = unittest.TestLoader().loadTestsFromTestCase(MinkowskiReductionTest)
unittest.TextTestRunner(verbosity=2).run(suite)

# Local Variables:
# compile-command: "cd ..; make test-minkowski"
# End:
