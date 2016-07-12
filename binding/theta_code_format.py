import itertools
from sage.misc.all import cached_function, mul
from sage.matrix.all import matrix
from sage.rings.all import QQ, PolynomialRing
from e8theta_degree3.gl3_repn import gl3_repn_module, matrix_var
from sage.modules.all import vector


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
    return pol.map_coefficients(lambda x: x.list()[0]).change_ring(QQ)


def _im_part(pol):
    return pol.map_coefficients(lambda x: x.list()[1]).change_ring(QQ)


@cached_function
def _pol_basis_factor_dct_and_ls(wt):
    '''
    wt: a list/tuple of non-increasing integers of length 3.
    Let M a corresponding repn of GL3.
    return a pair (d, l)
    l: a list of polynomials obtained from prime factors of M.basis_as_pol().
    d: dict s.t whose keys are M.basis() and
    f => [const, (b1, t1), (b2, t2), ... ],
    where f = const b1^t1 * b2^t2 * ..., and b1, b2, ... in l.
    '''
    M = gl3_repn_module(wt)
    basis = M.basis()
    facs = [(pol, pol.factor()) for pol in basis]
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


def _bideterminant_prime_factors_dict(mat, wt):
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
    Returns a pair of dicts (d, subs_dct).
    subs_dct is a dict s.t.
    a => (rli, imi)
    where a is a factorization of baisis_as_pol,
    rli, imi are variables and omega is the gen of imag_quad.
    The set of keys of d is basis of the corresponding representation.
    Its value at x is a pair (f, g) of polynomials of rl0, im0, rl1, im1, ... s.t.
    x.subs({a => rli + omega * imi}) = f + g * omega.
    '''
    d, l = _pol_basis_factor_dct_and_ls(wt)
    n = len(l)
    names = itertools.chain(
        *[[names_base[0] + str(a), names_base[1] + str(a)] for a in range(n)])
    R = PolynomialRing(imag_quad, names=list(names))
    omega = imag_quad.gen()
    r_gens = R.gens()
    subs_dct = {fc: gns[0] + omega * gns[1]
                for gns, fc in zip([r_gens[a:a + 2] for a in range(0, 2 * n, 2)], l)}

    def _subs(ls):
        return ls[0] * mul(subs_dct[a] ** b for a, b in ls[1:])

    res = {k: _subs(v) for k, v in d.items()}
    return ({k: (_rl_part(v), _im_part(v)) for k, v in res.items()},
            {k: (_rl_part(v), _im_part(v)) for k, v in subs_dct.items()})


def _init_code(wt, mat):
    pass


def code_format(wt, mat, real_part=True):
    '''
    wt: non-increasing list/tuple of non-negative integers of length 3.
    mat: 3 * 8 matrix with mat * mat.transpose() = 0 with coefficients in
    an imaginary quadratic field.
    real_part: boolian. If False, the imaginary part of the Fourier coefficient
    will be computed.
    Here for alpha = a + b * omega in K (an imaginary quadratic field with the
    generator omega), the real part and imaginary part of alpha are
    a and b respectively.
    '''

    header = '''
#include "e8vectors.h"

inline int inner_prod(int s[8], int t[8])
{
  return ((2*s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7]) * t[0] +
          (s[0] + 2*s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + 2*s[7]) * t[1] +
          (s[0] + s[1] + 2*s[2] + s[3] + s[4] + s[5] + s[6] + 2*s[7]) * t[2] +
          (s[0] + s[1] + s[2] + 2*s[3] + s[4] + s[5] + s[6] + 2*s[7]) * t[3] +
          (s[0] + s[1] + s[2] + s[3] + 2*s[4] + s[5] + s[6] + 2*s[7]) * t[4] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + 2*s[5] + s[6] + 2*s[7]) * t[5] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + 2*s[6] + 2*s[7]) * t[6] +
          (s[0] + 2*s[1] + 2*s[2] + 2*s[3] + 2*s[4] + 2*s[5] + 2*s[6] + 4*s[7]) * t[7]);
}

'''


def code_format(wt, mat):
    '''
    mat: 3 * 8 matrix with mat * mat.transpose() = 0 with coefficients in
    an imaginary quadratic field.
    '''
