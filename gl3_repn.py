import itertools
import operator

from e8theta_degree3.repn import ReplSpaceElement
from e8theta_degree3.young_tableau import (YoungTableu, poly_repn_dim,
                                           semistandard_young_tableaux)
from sage.matrix.all import identity_matrix, matrix
from sage.misc.all import cached_function, cached_method, flatten, mul
from sage.modules.all import vector
from sage.rings.all import QQ, PolynomialRing
from sage.parallel.all import fork
from degree2.utils import find_linearly_indep_indices


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
    def left_tableau(self):
        return self._a

    @property
    def right_tableau(self):
        return self._b

    @cached_method
    def determinants(self):
        m = matrix_var()
        res = []
        for l1, l2 in zip(self.left_tableau.col_numbers, self.right_tableau.col_numbers):
            if l1 and l2:
                res.append(m.matrix_from_rows_and_columns(
                    [i - 1 for i in l1], [j - 1 for j in l2]).det())
            else:
                # empty matrix
                res.append(matrix(QQ, []))
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

    def subs_and_compute_pol(self, d):
        return mul(a.subs(d) for a in self.determinants())

    def __hash__(self):
        return hash(('BiDeterminant', self.as_pol()))

    def __eq__(self, other):
        if isinstance(other, BiDeterminant):
            return self.as_pol() == other.as_pol()
        else:
            return False


class BiDetAsstoSSYT(BiDeterminant):

    def __init__(self, b, wt):
        t = _t_lambda(wt)
        super(BiDetAsstoSSYT, self).__init__(t, b)

    def element_weight(self):
        d = {i + 1: v for i, v in enumerate(identity_matrix(QQ, 3).columns())}
        return tuple(sum(d[a] for a in flatten(self.right_tableau.col_numbers)))


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
        R = matrix_var().base_ring()
        return [R(b.as_pol()) for b in self.basis()]

    @cached_method
    def basis(self):
        d = {i + 1: v for i, v in enumerate(identity_matrix(QQ, 3).columns())}
        ssyt = list(semistandard_young_tableaux(3, self.wt))
        if self.wt[0] != 0:
            ssyt = list(reversed(
                sorted(ssyt,
                       key=lambda x: (list(sum(d[a] for a in flatten(x.col_numbers))) +
                                      sum([a[self.wt[-1]:] for a in x.row_numbers], [])))))
        return [BiDetAsstoSSYT(a, self.wt) for a in ssyt]

    @cached_method
    def linearly_indep_tpls(self):
        kys = reduce(operator.add, [a.dict().keys()
                                    for a in self.basis_as_pol()], [])
        kys = list(set(kys))
        vecs = [[a[t] for a in self.basis_as_pol()] for t in kys]
        return [kys[i] for i in find_linearly_indep_indices(vecs, self.dimension())]

    @cached_method
    def _transform_mat_inv(self):
        m = matrix([[b[t] for t in self.linearly_indep_tpls()] for b in self.basis_as_pol()])
        return m**(-1)

    def to_vector(self, a):
        '''
        a: an element of the parent of matrix_var().
        Return vector corresponding to a.
        '''
        R = matrix_var().base_ring()
        a = R(a)
        v = vector([a[t] for t in self.linearly_indep_tpls()])
        m = self._transform_mat_inv()
        return v * m

    def to_pol(self, v):
        return sum(a * b for a, b in zip(v, self.basis_as_pol()))

    def matrix_representaion(self, g):
        '''
        g: matrix of size 3.
        Return matrix representation of the left action of g by self.basis_as_pol().
        Here we take the matrix representation rho(g) so that
        (b1(Xg), ..., bm(Xg)) = (b1(X), ..., bm(X)) rho(g.transpose()).transpose(),
        where b1, ..., bm are basis as polynomials.
        '''
        d = matrix_var_right_mul_dict(g.transpose())
        bs_acted = [a.subs_and_compute_pol(d) for a in self.basis()]
        m = matrix([[a[t] for t in self.linearly_indep_tpls()] for a in bs_acted])
        return m * self._transform_mat_inv()

    def __call__(self, v):
        '''
        v: vector with length self.dimension()
        Create an element which corresponds to GL3RepnElement
        '''
        if len(v) != self.dimension():
            raise ValueError
        return GL3RepnElement(v, self.wt)


@fork                           # to reduce memory usage
def matrix_repn(M, g):
    d = matrix_var_right_mul_dict(g.transpose())
    bs_acted = [a.subs_and_compute_pol(d) for a in M.basis()]
    m = matrix([[a[t] for t in M.linearly_indep_tpls()] for a in bs_acted])
    return m * M._transform_mat_inv()


class GL3RepnElement(ReplSpaceElement):

    def __init__(self, v, wt):
        super(GL3RepnElement, self).__init__(v, wt)
        self._parent = gl3_repn_module(wt)

    def left_action(self, g):
        det_wt = self.weight[-1]
        non_det_wt = tuple([a - det_wt for a in self.weight])
        M = gl3_repn_module(non_det_wt)
        m = M.matrix_representaion(g)
        return GL3RepnElement(m * self.vector * g.det()**det_wt, self.weight)

    def parent(self):
        return self._parent

    def dimension(self):
        return self.parent().dimension()


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
