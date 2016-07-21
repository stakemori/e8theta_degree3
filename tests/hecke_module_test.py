import unittest

from sage.matrix.all import matrix
from sage.rings.integer_ring import ZZ

from ..gl3_repn import GL3RepnElement
# from unittest import skip
from ..hecke_module import (HalfIntMatElement, _gl3_coset_gamma0,
                            _index_of_gamma_0_gl_n, spinor_l_euler_factor)


def alphas_ls():
    return [[0, 1, 2], [0, 1, 3], [0, 2, 3], [0, 0, 1], [0, 0, 2], [0, 1, 1], [0, 2, 2]]


class CosetTest(unittest.TestCase):

    def assert_index(self, alphas, p):
        self.assertEqual(len(_gl3_coset_gamma0(alphas, p)),
                         _index_of_gamma_0_gl_n(alphas, p))

    # @skip("OK")
    def test_indices(self):
        for alphas in alphas_ls():
            self.assert_index(alphas, 2)
            self.assert_index(alphas, 3)

    # @skip("OK")
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

    # @skip("OK")
    def test_distinct_cosets(self):
        for p in [ZZ(2), ZZ(3)]:
            for alphas in alphas_ls():
                print alphas, p
                gs = _gl3_coset_gamma0(alphas, p)
                self.assertTrue(self.check_distinct_cosets(alphas, p, gs))

    def test_spinor_l_miyawaki_at_2(self):
        miyawaki_fc = [([ZZ(1), ZZ(1) / ZZ(2), ZZ(1) / ZZ(2), ZZ(1) / ZZ(2),
                         ZZ(1), ZZ(1) / ZZ(2), ZZ(1) / ZZ(2), ZZ(1) / ZZ(2), ZZ(3)],
                        ZZ(70495764480)),
                       ([ZZ(2), ZZ(1), ZZ(1), ZZ(1), ZZ(2), ZZ(1), ZZ(1), ZZ(1), ZZ(2)],
                        -ZZ(6995218268160)),
                       ([ZZ(1), ZZ(0), ZZ(0), ZZ(0), ZZ(1), ZZ(0), ZZ(0), ZZ(0), ZZ(1)],
                        ZZ(8705802240)),
                       ([ZZ(1), ZZ(1) / ZZ(2), ZZ(1) / ZZ(2), ZZ(1) / ZZ(2), ZZ(1), ZZ(1) / ZZ(2),
                         ZZ(1) / ZZ(2), ZZ(1) / ZZ(2), ZZ(1)], ZZ(53084160)),
                       ([ZZ(2), ZZ(0), ZZ(0), ZZ(0), ZZ(2), ZZ(0), ZZ(0), ZZ(0), ZZ(2)],
                        -ZZ(361848813649920)),
                       ([ZZ(1), ZZ(0), ZZ(0), ZZ(0), ZZ(1), ZZ(0), ZZ(0), ZZ(0), ZZ(2)],
                        -ZZ(53508833280)),
                       ([ZZ(1), ZZ(0), ZZ(0), ZZ(0), ZZ(3), ZZ(1), ZZ(0), ZZ(1), ZZ(3)],
                        ZZ(93047614341120))]
        d = {HalfIntMatElement(matrix(3, a)): GL3RepnElement([b], (12, 12, 12))
             for a, b in miyawaki_fc}
        pl = spinor_l_euler_factor(2, d)
        t = pl.parent().gens()[0]
        delta_pl = 1 + 24 * t + 2**11 * t**2
        self.assertEqual(pl, (delta_pl.subs({t: 2**9 * t}) *
                              delta_pl.subs({t: 2**10 * t}) *
                              (1 + 2**6 * 171 * t - 2**17 * 10831 * t**2 +
                               2**36 * 171 * t**3 + 2**60 * t**4)))


suite = unittest.TestLoader().loadTestsFromTestCase(CosetTest)
unittest.TextTestRunner(verbosity=2).run(suite)
