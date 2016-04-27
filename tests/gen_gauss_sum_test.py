import unittest
from sage.matrix.all import diagonal_matrix
from sage.rings.all import FiniteField, ZZ
from ..hecke_module import _gen_gauss_sum_direct_way, _generalized_gauss_sum


class GenGaussSum(unittest.TestCase):

    def assert_gauss_sums(self, diag, p, r, verbose=False):
        N = diagonal_matrix(diag)
        N = N.change_ring(FiniteField(p))
        if verbose:
            print "diag:", diag
            print "r: ", r
        self.assertEqual(_generalized_gauss_sum(N, p, r),
                         _gen_gauss_sum_direct_way(N, p, r))

    def test_gauss_sums(self):
        N = diagonal_matrix([1, 2, 0])
        p = ZZ(3)
        N = N.change_ring(FiniteField(p))
        for diag in [[1, 2, 0], [1, 0, 0], [2, 0, 0], [1, 0, 0], [0, 0, 0],
                     [1, 1, 1], [1, 1, 2]]:
            for r in [1, 2, 3]:
                self.assert_gauss_sums(diag, ZZ(3), r, verbose=True)

        for diag in [[1, 2], [1, 1], [2, 0], [1, 0], [0, 0]]:
            for r in [1, 2]:
                self.assert_gauss_sums(diag, ZZ(3), r, verbose=True)


suite = unittest.TestLoader().loadTestsFromTestCase(GenGaussSum)
unittest.TextTestRunner(verbosity=2).run(suite)

# Local Variables:
# compile-command: "cd ..; make test-gen-gauss-sum"
# End:
