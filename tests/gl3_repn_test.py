import unittest
from sage.modules.all import random_vector
from sage.rings.all import QQ
from ..deg3_repn import GL3RepnModule


class Gl3RepnTest(unittest.TestCase):

    def test_to_vector(self):
        M = GL3RepnModule([3, 2, 1])
        d = M.dimension()
        for _ in range(30):
            v = random_vector(QQ, d)
            pol = M.to_pol(v)
            self.assertEqual(M.to_pol(M.to_vector(pol)), pol)


suite = unittest.TestLoader().loadTestsFromTestCase(Gl3RepnTest)
unittest.TextTestRunner(verbosity=2).run(suite)
