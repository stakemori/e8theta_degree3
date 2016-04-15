import unittest
# from unittest import skip
from ..hecke_module import _index_of_gamma_0_gl_n, _gl3_coset_gamma0


class CosetTest(unittest.TestCase):

    def assert_index(self, alphas, p):
        self.assertEqual(len(_gl3_coset_gamma0(alphas, p)),
                         _index_of_gamma_0_gl_n(alphas, p))

    def test_indices(self):
        self.assert_index([0, 1, 2], 2)

suite = unittest.TestLoader().loadTestsFromTestCase(CosetTest)
unittest.TextTestRunner(verbosity=2).run(suite)
