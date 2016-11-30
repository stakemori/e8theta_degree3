import os
import unittest
from unittest import skip

from e8theta_degree3.gl3_repn import (GL3RepnElement, gl3_repn_module,
                                      matrix_var, matrix_var_right_mul_dict)
from e8theta_degree3.hecke_module import HalfIntMatElement, tp_action_fc_alist
from e8theta_degree3.results.data.data_utils import data_dir, dict_sum
from sage.all import (QQ, SR, ZZ, NumberField, cached_function, load, matrix,
                      random_matrix)


@cached_function
def trans_mat_dict():
    return load(os.path.join(data_dir(), "trans_mats.sobj"))


def left_action(m, g, wt):
    w = tuple([a - wt[-1] for a in wt])
    v = m.vector
    d = matrix_var_right_mul_dict(g)
    M = gl3_repn_module(w)
    return (sum(a * b.subs_and_compute_pol(d) for a, b in zip(v, M.basis()))
            * (matrix_var().det()**wt[-1]) * g.det()**wt[-1])


class TransMatTest(unittest.TestCase):

    def trans_mat_test_gen(self, wt, mat):
        wt = tuple([a - wt[-1] for a in wt])
        M = gl3_repn_module(wt)
        for _ in range(10):
            g = random_matrix(QQ, 3, algorithm='unimodular')
            m = (M.matrix_representaion(g.transpose())).transpose()
            n = M.matrix_representaion(g)
            self.assertEqual(m * mat, mat * n)

    @skip("OK")
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
        a = pl1.lc()//pl.lc()
        self.assertEqual(pl1, pl * a)

    def test_hecke_eigen_14_13_5(self):
        l = load(os.path.join(data_dir(), "wt14_13_5.sobj"))
        f = {HalfIntMatElement(t): GL3RepnElement(v, (14, 13, 5)) for t, v in l}
        self.assert_hecke_eigen(f, (14, 13, 5))

    def test_hecke_eigen_16_13_7(self):
        dcts16_13_7 = load(os.path.join(data_dir(), "wt16_13_7_dicts.sobj"))
        f0 = dict_sum([ZZ(32539262976), ZZ(252907706779)], dcts16_13_7)
        f1 = dict_sum([ZZ(441367252992), ZZ(1352603388073)], dcts16_13_7)
        self.assert_hecke_eigen(f0, (16, 13, 7))
        self.assert_hecke_eigen(f1, (16, 13, 7))

    def test_hecke_eigen_16_16_14(self):
        x = SR.var("x")
        K = NumberField(x**2 - ZZ(8640)*x - ZZ(454569984), names="a")
        a = K.gen()
        fc_dct1, fc_dct2 = load(os.path.join(data_dir(), "wt16_16_14.sobj"))
        F = {k: ((-3280338755058*a + 321567813112379328) * fc_dct1[k] -
                 1900007009496106555 * fc_dct2[k]) for k in fc_dct1}
        self.assert_hecke_eigen(F, (16, 16, 14))

    def test_hecke_eigen_18_13_5(self):
        bs = load(os.path.join(data_dir(), "wt18_13_5_dicts.sobj"))
        f0 = dict_sum((ZZ(10149020898562), -ZZ(390753)), bs)
        f1 = dict_sum((ZZ(30563021287421), ZZ(321984)), bs)
        self.assert_hecke_eigen(f0, (18, 13, 5))
        self.assert_hecke_eigen(f1, (18, 13, 5))

    def test_hecke_eigen_18_17_5(self):
        bs = load(os.path.join(data_dir(), "wt18_17_5_dicts.sobj"))
        f0 = dict_sum((-ZZ(27411940570812415758258813595452421519092205),
                       ZZ(426844081323506188348315458383164538880),
                       ZZ(1108258366195376623065993820933390921728)), bs)
        x = SR.var("x")
        K = NumberField(x**2 - 122112*x - 6040977408, names="a")
        f1 = dict_sum((115755188162636478306143242492773391279012500590415297895*K.gen() -
                       280871869782212043390032985770730420043620758653959144573711360,
                       -3520271549093870111699822711292893500461035796234240*K.gen() +
                       1792446151032915388102568677460511010394423622703062712320,
                       -23730541121480814714846359164016657918248298142879744*K.gen() +
                       8813953574896657677285798006856812369943698972051863764992),
                      bs)
        self.assert_hecke_eigen(f0, (18, 17, 5))
        self.assert_hecke_eigen(f1, (18, 17, 5))


suite = unittest.TestLoader().loadTestsFromTestCase(TransMatTest)
unittest.TextTestRunner(verbosity=2).run(suite)
