import os
import unittest
# from unittest import skip

from e8theta_degree3.gl3_repn import gl3_repn_module, matrix_var_right_mul_dict, GL3RepnElement
from e8theta_degree3.results.data.data_utils import data_dir
from e8theta_degree3.hecke_module import tp_action_fc_alist, HalfIntMatElement
from sage.all import ZZ, QQ, cached_function, load, random_matrix, matrix


@cached_function
def trans_mat_dict():
    return load(os.path.join(data_dir(), "trans_mats.sobj"))


def left_action(m, g, wt):
    v = m.vector
    d = matrix_var_right_mul_dict(g)
    M = gl3_repn_module(wt)
    return sum(a * b.subs_and_compute_pol(d) for a, b in zip(v, M.basis()))


class TransMatTest(unittest.TestCase):

    def trans_mat_test_gen(self, wt, mat):
        wt = tuple([a - wt[-1] for a in wt])
        M = gl3_repn_module(wt)
        for _ in range(10):
            g = random_matrix(QQ, 3, algorithm='unimodular')
            m = (M.matrix_representaion(g.transpose())).transpose()
            n = M.matrix_representaion(g)
            self.assertEqual(m * mat, mat * n)

    # @skip("OK")
    def test_trans_mat(self):
        for k, b in trans_mat_dict().items():
            self.trans_mat_test_gen(k, b)

    def assert_hecke_eigen(self, f, wt):
        B = trans_mat_dict()[wt]
        M = gl3_repn_module(wt)
        f_transed = {k: M(B * v.vector) for k, v in f.items()}
        T0 = HalfIntMatElement(matrix(QQ, [[2, 1, 1], [1, 2, 1], [1, 1, 2]])/QQ(2))
        alst_t2 = tp_action_fc_alist(ZZ(2), T0)
        pl = sum(left_action(f_transed[S], g, wt) * a for S, a, g in alst_t2)
        pl1 = M.to_pol(f_transed[T0].vector)
        a = pl1//pl
        self.assertEqual(pl1, pl * a)

    def test_hecke_eigen(self):
        l = load(os.path.join(data_dir(), "wt14_13_5.sobj"))
        f = {HalfIntMatElement(t): GL3RepnElement(v, (14, 13, 5)) for t, v in l}
        self.assert_hecke_eigen(f, (14, 13, 5))

suite = unittest.TestLoader().loadTestsFromTestCase(TransMatTest)
unittest.TextTestRunner(verbosity=2).run(suite)
