import unittest

from sage.rings.integer_ring import ZZ

# from unittest import skip
from ..hecke_module import _gl3_coset_gamma0, _index_of_gamma_0_gl_n


def alphas_ls():
    return [[0, 1, 2], [0, 1, 3], [0, 2, 3], [0, 0, 1], [0, 0, 2], [0, 1, 1], [0, 2, 2]]


class CosetTest(unittest.TestCase):

    def assert_index(self, alphas, p):
        self.assertEqual(len(_gl3_coset_gamma0(alphas, p)),
                         _index_of_gamma_0_gl_n(alphas, p))

    def test_indices(self):
        for alphas in alphas_ls():
            self.assert_index(alphas, 2)
            self.assert_index(alphas, 3)

    def check_distinct_cosets(self, alphas, p, gs):
        res = True
        p = ZZ(p)
        a0, a1, a2 = alphas
        for i, g in enumerate(gs):
            for j, h in enumerate(gs):
                if i < j:
                    m = g * h**(-1)
                    in_gam0 = ((ZZ(m[0, 1]) / p**(a1 - a0)) in ZZ and
                               (ZZ(m[0, 2]) / p**(a2 - a0)) in ZZ and
                               (ZZ(m[1, 2]) / p**(a2 - a1)) in ZZ)
                    res = res and not in_gam0
                    if not res:
                        print g
                        print h
                        return False
        return True

    def test_distinct_cosets(self):
        for p in [ZZ(2), ZZ(3)]:
            for alphas in alphas_ls():
                print alphas, p
                gs = _gl3_coset_gamma0(alphas, p)
                self.assertTrue(self.check_distinct_cosets(alphas, p, gs))


suite = unittest.TestLoader().loadTestsFromTestCase(CosetTest)
unittest.TextTestRunner(verbosity=2).run(suite)
