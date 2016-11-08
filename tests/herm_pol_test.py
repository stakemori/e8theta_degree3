import unittest

from e8theta_degree3.binding.theta_code_format import (_bideterminant_prime_factors_dict,
                                                       _s_t_u_ring, matrix_var,
                                                       _pol_basis_as_polof_factors)
from sage.matrix.all import matrix
from sage.rings.all import QuadraticField
from sage.misc.all import mul


class TestHermonicPol(unittest.TestCase):

    def test_bi_determinant_sub(self):
        R = _s_t_u_ring()
        stu_mt = matrix(R, 3, R.gens())
        stu_mt = matrix(stu_mt.rows()).transpose()
        K = QuadraticField(-1, name='i')
        i = K.gen()
        mat = matrix([[1, 0, 0, i, 0, 0, 0, 0],
                      [0, 1, 0, 0, i, 0, 0, 0],
                      [0, 0, 1, 0, 0, i, 0, 0]])
        subs_dct = dict(zip(matrix_var().list(), (mat * stu_mt).list()))
        d = _bideterminant_prime_factors_dict(mat, (4, 2, 2))
        for bd, (rl, im) in d.items():
            self.assertEqual(bd.subs(subs_dct), rl + im * i)

    def test_pol_basis_as_polof_factors(self):
        K = QuadraticField(-1, name='i')
        d, subs_d = _pol_basis_as_polof_factors((6, 2, 1), K)
        subs_d = {k: rl + im * K.gen() for k, (rl, im) in subs_d.items()}

        def _sub(f):
            return mul(subs_d[p] ** i for p, i in f.factor().pols) * f.factor().const

        for f, (rl, im) in d.items():
            self.assertEqual(_sub(f), rl + im * K.gen())

suite = unittest.TestLoader().loadTestsFromTestCase(TestHermonicPol)
unittest.TextTestRunner(verbosity=2).run(suite)
