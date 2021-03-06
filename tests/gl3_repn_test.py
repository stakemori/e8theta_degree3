import unittest

from sage.matrix.all import diagonal_matrix, random_matrix
from sage.modules.all import random_vector
from sage.rings.all import QQ

from ..gl3_repn import gl3_repn_module


class Gl3RepnTest(unittest.TestCase):

    def test_to_vector(self):
        M = gl3_repn_module((3, 2, 1))
        d = M.dimension()
        for _ in range(30):
            v = random_vector(QQ, d)
            pol = M.to_pol(v)
            self.assertEqual(M.to_pol(M.to_vector(pol)), pol)

    def test_matrix_repn(self):
        M = gl3_repn_module((5, 3, 1))
        for _ in range(30):
            while True:
                g1 = random_matrix(QQ, 3)
                if not g1.is_singular():
                    break
            while True:
                g2 = random_matrix(QQ, 3)
                if not g2.is_singular():
                    break
            self.assertEqual(M.matrix_representaion(g1) * M.matrix_representaion(g2),
                             M.matrix_representaion(g1 * g2))

    def test_left_action(self):
        M = gl3_repn_module((12, 12, 12))
        v = M((3, ))
        self.assertEqual(v.left_action(diagonal_matrix([5, 7, 11]))[0],
                         3 * (5 * 7 * 11)**12)
        M1 = gl3_repn_module((1, 0, 0))
        v1 = M1((1, 2, 3))
        self.assertEqual(list(v1.left_action(diagonal_matrix([5, 7, 11])).vector),
                         [5, 2 * 7, 3 * 11])

    def test_sort_basis(self):

        def assert_sort(wt):
            M = gl3_repn_module(wt)
            M1 = gl3_repn_module(tuple([a - wt[-1] for a in wt]))
            for b1, b2 in zip(M.basis(), M1.basis()):
                self.assertTrue(all(a1[wt[-1]:] == a2 for a1, a2 in
                                    zip(b1.right_tableau.row_numbers,
                                        b2.right_tableau.row_numbers)))
        assert_sort((10, 8, 6))
        assert_sort((23, 4, 3))
        assert_sort((12, 8, 2))

suite = unittest.TestLoader().loadTestsFromTestCase(Gl3RepnTest)
unittest.TextTestRunner(verbosity=2).run(suite)
