import unittest
from os.path import join
from sage.matrix.all import matrix
from sage.misc.all import load
from sage.rings.all import ZZ
from sage.quadratic_forms.all import QuadraticForm
from e8theta_degree3.results.data.data_utils import data_dir
from e8theta_degree3.hecke_module import HalfIntMatElement
from e8theta_degree3.gl3_repn import GL3RepnElement


class FourierCoeffTest(unittest.TestCase):

    def assert_unimodular_invariance(self, d):
        for t, v in d.items():
            q = QuadraticForm(ZZ, 2 * t.T)
            gens = q.automorphism_group().gens()
            for g in gens:
                g = matrix(3, g.list())
                self.assertEqual(t.right_action(g), t)
                self.assertEqual(v.left_action(g.transpose()).vector, v.vector)

    def test_wt14_13_5(self):
        l = load(join(data_dir(), "wt14_13_5.sobj"))
        d = dict([(HalfIntMatElement(t), GL3RepnElement(v, (14, 13, 5))) for t, v in l])
        self.assert_unimodular_invariance(d)

    def test_wt16_16_14(self):
        ds = load(join(data_dir(), "wt16_16_14.sobj"))
        for d in ds:
            self.assert_unimodular_invariance(d)

    def test_wt16_13_7(self):
        ds = load(join(data_dir(), "wt16_13_7_dicts.sobj"))
        for d in ds:
            self.assert_unimodular_invariance(d)

    def test_wt18_13_5(self):
        ds = load(join(data_dir(), "wt18_13_5_dicts.sobj"))
        for d in ds:
            self.assert_unimodular_invariance(d)

    def test_wt18_17_5(self):
        ds = load(join(data_dir(), "wt18_17_5_dicts.sobj"))
        for d in ds:
            self.assert_unimodular_invariance(d)


suite = unittest.TestLoader().loadTestsFromTestCase(FourierCoeffTest)
unittest.TextTestRunner(verbosity=2).run(suite)
