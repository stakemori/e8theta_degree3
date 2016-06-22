import operator
import itertools

from sage.rings.all import PolynomialRing, QQ
from sage.matrix.all import matrix
from sage.misc.all import cached_function, cached_method, mul
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

    return GL3RepnElement


@cached_function
def _s_t_u_ring(base_ring=None):
    if base_ring is None:
        base_ring = QQ
    R = PolynomialRing(
        base_ring, names=("s0, s1, s2, s3, s4, s5, s6, s7,"
                          "t0, t1, t2, t3, t4, t5, t6, t7,"
                          "u0, u1, u2, u3, u4, u5, u6, u7"))
    return R


def _conj(pol):
    return pol.map_coefficients(lambda c: c.conjugate())


def _rl_part(pol):
    return ((pol + _conj(pol)) / QQ(2)).change_ring(QQ)


def _im_part(pol):
    K = pol.base_ring()
    i = K.gen()
    return ((pol - _conj(pol)) / (QQ(2) * i)).change_ring(QQ)


def _normalized_factor(pol):
    '''
    Prime factors of pol may differ by constant.
    '''
    l = [(a / a.lc(), b) for a, b in pol.factor() if not a.is_constant()]
    a = pol.lc() / mul(a ** b for a, b in l).lc()
    return NormFactorELt(a, l)


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


@cached_function
def _pol_basis_factor_dct_and_ls(wt):
    '''
    wt: a list/tuple of non-increasing integers of length 3.
    Let M a corresponding repn of GL3.
    return a pair (d, l)
    l: a list of polynomials obtained from prime factors of M.basis_as_pol().
    d: dict s.t whose keys are M.basis_as_pol() and
    f => [const, (b1, t1), (b2, t2), ... ],
    where f = const b1^t1 * b2^t2 * ..., and b1, b2, ... in l.
    '''
    M = gl3_repn_module(wt)
    basis = M.basis_as_pol()
    facs = [(pol, _normalized_factor(pol)) for pol in basis]
    l = list(
        set(itertools.chain(*[[a for a, _ in fc.pols] for _, fc in facs])))
    return [{pol: [fc.const] + fc.pols for pol, fc in facs}, l]


@cached_function
def euclidean_basis():
    basis_vecs = [(QQ(1) / QQ(2), QQ(1) / QQ(2), QQ(1) / QQ(2),
                   QQ(1) / QQ(2), QQ(1) / QQ(2), QQ(1) / QQ(2),
                   QQ(1) / QQ(2), QQ(1) / QQ(2)),
                  (0, 1, 0, 0, 0, 0, 0, 1),
                  (0, 0, 1, 0, 0, 0, 0, 1),
                  (0, 0, 0, 1, 0, 0, 0, 1),
                  (0, 0, 0, 0, 1, 0, 0, 1),
                  (0, 0, 0, 0, 0, 1, 0, 1),
                  (0, 0, 0, 0, 0, 0, 1, 1),
                  (0, 0, 0, 0, 0, 0, 0, 2)]
    return [vector([QQ(a) for a in v1]) for v1 in basis_vecs]


def to_eulidian_vec(t):
    return sum([a * b for a, b in zip(t, euclidean_basis())])


def _bideterminants_dict(mat, wt):
    '''
    wt: a list/tuple of non-increasing integers of length 3
    mat: 3 * 8 matrix with mat * mat.transpose() = 0 with coefficients in
    an imaginary quadratic field.
    Return a dict
    a: (real_part, imag_part) as polynomials of _s_t_u_ring(QQ),
    where a is a prime factor of polynomial basis
    '''
    R = _s_t_u_ring()
    stu_mt = matrix(R, 3, R.gens())
    stu_mt = matrix([to_eulidian_vec(v) for v in stu_mt.rows()]).transpose()
    subs_dct = dict(zip(matrix_var().list(), (mat * stu_mt).list()))
    _, l = _pol_basis_factor_dct_and_ls(wt)
    d = {a: a.subs(subs_dct) for a in l}
    return {a: (_rl_part(b), _im_part(b)) for a, b in d.items()}


def _pol_basis_as_polof_factors(wt, imag_quad, names_base=('rl', 'im')):
    '''
    wt: weight of repn of GL3
    imag_quad: imaginary quadratic field
    Returns a pair of dicts
    '''
    d, l = _pol_basis_factor_dct_and_ls(wt)
    n = len(l)
    names = itertools.chain(
        *[[names_base[0] + str(a), names_base[1] + str(a)] for a in range(n)])
    R = PolynomialRing(imag_quad, names=list(names))
    omega = imag_quad.gen()
    r_gens = R.gens()
    subs_dct = {fc: gns[0] + omega * gns[1]
                for gns, fc in zip([r_gens[a:a + 2] for a in range(0, n, 2)], l)}

    def _subs(ls):
        return ls[0] * mul(subs_dct[a] ** b for a, b in ls[1:])

    res = {k: _subs(v) for k, v in d.items()}
    return ({k: (_rl_part(v), _im_part(v)) for k, v in res.items()},
            subs_dct)
