import unittest
from sage.modules.all import random_vector
from sage.rings.all import QQ
from sage.matrix.all import random_matrix

from ..gl3_repn import GL3RepnModule


class Gl3RepnTest(unittest.TestCase):

    def test_to_vector(self):
        M = GL3RepnModule([3, 2, 1])
        d = M.dimension()
        for _ in range(30):
            v = random_vector(QQ, d)
            pol = M.to_pol(v)
            self.assertEqual(M.to_pol(M.to_vector(pol)), pol)

    def test_matrix_repn(self):
        M = GL3RepnModule([5, 3, 1])
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


suite = unittest.TestLoader().loadTestsFromTestCase(Gl3RepnTest)
unittest.TextTestRunner(verbosity=2).run(suite)
