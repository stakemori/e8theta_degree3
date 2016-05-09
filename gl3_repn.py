import operator

from sage.rings.all import PolynomialRing, QQ
from sage.matrix.all import matrix
from sage.misc.all import cached_function, cached_method
from sage.modules.all import vector
from e8theta_degree3.young_tableau import YoungTableu, semistandard_young_tableaux, poly_repn_dim
from e8theta_degree3.repn import ReplSpaceElement
from e8theta_degree3.utils import find_linearly_indep_indices


@cached_function
def matrix_var(base_field=QQ):
    R = PolynomialRing(base_field, names=[
                       'x%s%s' % (i, j) for i in range(3) for j in range(3)])
    return matrix(3, R.gens())


def matrix_var_right_mul_dict(g):
    '''
    Return a dictionary which sends matrix_var() => matrix_var() * g.
    '''
    m = matrix_var()
    m1 = m * g
    return {m[i, j]: m1[i, j] for i in range(3) for j in range(3)}


def left_action_as_pol(pol, g):
    d = matrix_var_right_mul_dict(g)
    return pol.subs(d)


def _bideterminant(a, b):
    '''
    a, b: an instance of young_tableau.YoungTableu
    '''
    m = matrix_var()
    res = 1
    for l1, l2 in zip(a.col_numbers, b.col_numbers):
        res *= m.matrix_from_rows_and_columns(
            [i - 1 for i in l1], [j - 1 for j in l2]).det()
    return res


def _t_lambda(wt):
    return YoungTableu(n=3, row_numbers=[[i + 1 for _ in range(a)]
                                         for i, a in enumerate(wt)])


class GL3RepnModule(object):

    def __init__(self, wt):
        self._wt = wt

    @property
    def wt(self):
        return self._wt

    @cached_method
    def dimension(self):
        return poly_repn_dim(self.wt)

    @cached_method
    def basis_as_pol(self):
        t = _t_lambda(self.wt)
        return [_bideterminant(t, a) for a in semistandard_young_tableaux(3, self.wt)]

    @cached_method
    def linearly_indep_tpls(self):
        kys = reduce(operator.add, [a.dict().keys()
                                    for a in self.basis_as_pol()], [])
        kys = list(set(kys))
        vecs = [[a[t] for a in self.basis_as_pol()] for t in kys]
        return [kys[i] for i in find_linearly_indep_indices(vecs, self.dimension())]

    @cached_method
    def _transform_mat(self):
        return matrix([[b[t] for t in self.linearly_indep_tpls()] for b in self.basis_as_pol()])

    def to_vector(self, a):
        '''
        an element of the parent of matrix_var().
        Return vector corresponding a.
        '''
        v = vector([a[t] for t in self.linearly_indep_tpls()])
        m = self._transform_mat()
        return v * m ** (-1)

    def to_pol(self, v):
        return sum(a * b for a, b in zip(v, self.basis_as_pol()))

    def matrix_representaion(self, g):
        '''
        g: matrix of size 3.
        Return matrix representation of the left action of g by self.basis_as_pol().
        '''
        d = matrix_var_right_mul_dict(g)
        bs_acted = [a.subs(d) for a in self.basis_as_pol()]
        return matrix([self.to_vector(a) for a in bs_acted]).transpose()


def element_constructor(wt):

    M = GL3RepnModule(wt)

    class GL3RepnElement(ReplSpaceElement):

        def left_action(self, g):
            pol = M.to_pol(self.vector)
            pol1 = left_action_as_pol(pol, g)
            return GL3RepnElement(M.to_vector(pol1))

        def parent(self):
            return M

    return GL3RepnElement
