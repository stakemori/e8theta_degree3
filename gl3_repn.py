import itertools
import operator

from e8theta_degree3.repn import ReplSpaceElement
from e8theta_degree3.utils import find_linearly_indep_indices
from e8theta_degree3.young_tableau import (YoungTableu, poly_repn_dim,
                                           semistandard_young_tableaux)
from sage.matrix.all import matrix
from sage.misc.all import cached_function, cached_method, mul
from sage.modules.all import vector
from sage.rings.all import QQ, PolynomialRing


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


class BiDeterminant(object):

    def __init__(self, a, b):
        self._a = a
        self._b = b

    def __repr__(self):
        return repr(self.factor())

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @cached_method
    def determinants(self):
        m = matrix_var()
        res = []
        for l1, l2 in zip(self.a.col_numbers, self.b.col_numbers):
            res.append(m.matrix_from_rows_and_columns(
                [i - 1 for i in l1], [j - 1 for j in l2]).det())
        return res

    def factor(self):
        const = QQ(1)
        dets = []
        for d in self.determinants():
            const = const * d.lc()
            dets.append(d * d.lc() ** (-1))
        l = [(k, len(list(v)))
             for k, v in itertools.groupby(sorted(dets), key=lambda x: x)]
        return NormFactorELt(const, l)

    def as_pol(self):
        return mul(a for a in self.determinants())

    def __hash__(self):
        return hash(('BiDeterminant', self.as_pol()))

    def __eq__(self, other):
        if isinstance(other, BiDeterminant):
            return self.as_pol() == other.as_pol()
        else:
            return False


def _t_lambda(wt):
    return YoungTableu(n=3, row_numbers=[[i + 1 for _ in range(a)]
                                         for i, a in enumerate(wt)])


@cached_function
def gl3_repn_module(wt):
    return GL3RepnModule(wt)


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
        return [b.as_pol() for b in self.basis()]

    @cached_method
    def basis(self):
        t = _t_lambda(self.wt)
        return [BiDeterminant(t, a) for a in semistandard_young_tableaux(3, self.wt)]

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
    '''
    wt: a list/tuple of non-increasing integers of length 3.
    Returns a sub class of ReplSpaceElement of given wt.
    Used for Hecke operators.
    '''

    M = GL3RepnModule(wt)

    class GL3RepnElement(ReplSpaceElement):

        def left_action(self, g):
            pol = M.to_pol(self.vector)
            pol1 = left_action_as_pol(pol, g)
            return GL3RepnElement(M.to_vector(pol1))

        def parent(self):
            return M

        @classmethod
        def dimension(cls):
            return M.dimension()

        @classmethod
        def weight(cls):
            return M.wt

        @classmethod
        def zero(cls):
            return cls(vector([QQ(0) for _ in range(M.dimension())]))

    return GL3RepnElement


class NormFactorELt(object):

    def __init__(self, c, facs):
        self._c = c
        self._facs = facs

    @property
    def const(self):
        return self._c

    @property
    def pols(self):
        return self._facs

    def __iter__(self):
        yield (self.const, 1)
        for a in self.pols:
            yield a

    def __repr__(self):
        l = [self.const] + self.pols
        return repr(l)
