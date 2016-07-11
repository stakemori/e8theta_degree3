# -*- coding: utf-8 -*-
from itertools import groupby
from sage.all import mul
from sage.arith.all import kronecker_symbol
from sage.matrix.all import (diagonal_matrix, matrix, block_diagonal_matrix,
                             identity_matrix, block_matrix)
from sage.misc.all import cached_function
from sage.rings.all import FiniteField, CyclotomicField, ZZ, QQ, PolynomialRing
from sage.quadratic_forms.all import least_quadratic_nonresidue, QuadraticForm
from sage.functions.all import sgn, floor, ceil
import itertools


def _index_of_gamma_0_gl_n(alpha, p):
    '''
    Returns delta(a1, ..., an) defined in Shimura, Euler products and Eisenstein
    series, pp 118, (15.1.7).
    '''
    if p in ZZ:
        p = ZZ(p)

    def _bn(n):
        return mul(1 - p ** (-i) for i in xrange(1, n + 1))

    e_r_ls = [(k, len(list(v)))
              for k, v in groupby(sorted(alpha), lambda x: x)]
    res = _bn(len(alpha)) / mul(_bn(r) for _, r in e_r_ls)
    for i, (ei, ri) in enumerate(e_r_ls):
        for j, (ej, rj) in enumerate(e_r_ls):
            if i < j:
                res *= p ** ((ej - ei) * ri * rj)
    return res


def _gl2_coset_gamma0(a, p):
    w = matrix([[0, -1],
                [1, 0]])
    for m12 in range(p ** a):
        yield matrix([[1, m12],
                      [0, 1]])
    for m21 in range(p ** (a - 1)):
        m = matrix([[1, 0],
                    [p * m21, 1]])
        yield w * m


def _gl3_coset_gamma0(alpha, p):
    r'''
    Let alpha = [a0, a1, a2] with a0 <= a1 <= a2,
    g = diag([p^a0, p^a1, p^a2]), and Gamma0 = g^(-1) GL3(Z) g ∧ GL3(Z).
    Return a complete set Gamma0 \ GL3(Z).
    '''
    if p in ZZ:
        p = ZZ(p)
    a0, a1, a2 = alpha
    if a0 < a1 < a2:
        return list(__gl3_coset_gamma0_distinct(a0, a1, a2, p))
    elif a0 == a1 and a1 < a2:
        return list(__gl3_coset_gamma0_2_1(a0, a2, p))
    elif a0 < a1 and a1 == a2:
        return list(__gl3_coset_gamma0_1_2(a0, a2, p))
    elif a0 == a1 == a2:
        return [identity_matrix(ZZ, 3)]
    else:
        raise ValueError


def __gl3_coset_gamma0_2_1(a1, a3, p):
    w23 = matrix([[1, 0, 0],
                  [0, 0, 1],
                  [0, 1, 0]])
    for m13 in range(p ** (a3 - a1 - 1)):
        for m23 in range(p ** (a3 - a1 - 1)):
            m = matrix([[1, 0, p * m13],
                        [0, 1, p * m23],
                        [0, 0, 1]])
            yield m

    for m32 in range(p ** (a3 - a1)):
        m = matrix([[1, 0, 0],
                    [0, 1, 0],
                    [0, m32, 1]])
        for g in _gl2_coset_gamma0(a3 - a1, p):
            n = block_diagonal_matrix(g, matrix([[1]]))
            yield w23 * m * n


def __gl3_coset_gamma0_1_2(a1, a2, p):
    w12 = matrix([[0, 1, 0],
                  [1, 0, 0],
                  [0, 0, 1]])

    for m12 in range(p ** (a2 - a1 - 1)):
        for m13 in range(p ** (a2 - a1 - 1)):
            m = matrix([[1, p * m12, p * m13],
                        [0, 1, 0],
                        [0, 0, 1]])
            yield m
    for m21 in range(p ** (a2 - a1)):
        m = matrix([[1, 0, 0],
                    [m21, 1, 0],
                    [0, 0, 1]])
        for g in _gl2_coset_gamma0(a2 - a1, p):
            n = block_diagonal_matrix(matrix([[1]]), g)
            yield w12 * m * n


def __gl3_coset_gamma0_distinct(a1, a2, a3, p):

    w12 = matrix([[0, 1, 0],
                  [1, 0, 0],
                  [0, 0, 1]])

    w23 = matrix([[1, 0, 0],
                  [0, 0, 1],
                  [0, 1, 0]])

    w13 = matrix([[0, 0, 1],
                  [0, 1, 0],
                  [1, 0, 0]])

    w123 = matrix([[0, 1, 0],
                   [0, 0, 1],
                   [1, 0, 0]])

    w132 = matrix([[0, 0, 1],
                   [1, 0, 0],
                   [0, 1, 0]])

    # w = 1
    for m12 in range(p ** (a2 - a1 - 1)):
        for m13 in range(p ** (a3 - a1 - 1)):
            for m23 in range(p ** (a3 - a2 - 1)):
                yield matrix([[1, p * m12, p * m13],
                              [0, 1, p * m23],
                              [0, 0, 1]])
    # w = (12)
    for m13 in range(p ** (a3 - a2 - 1)):
        for m21 in range(p ** (a2 - a1)):
            for m23 in range(p ** (a3 - a1 - 1)):
                m = matrix([[1, 0, p * m13],
                            [m21, 1, p * m23],
                            [0, 0, 1]])
                yield w12 * m
    # w = (23)
    for m12 in range(p ** (a3 - a1 - 1)):
        for m13 in range(p ** (a2 - a1 - 1)):
            for m32 in range(p ** (a3 - a2)):
                m = matrix([[1, p * m12, p * m13],
                            [0, 1, 0],
                            [0, m32, 1]])
                yield w23 * m

    # w = (13)
    for m21 in range(p ** (a3 - a2)):
        for m31 in range(p ** (a3 - a1)):
            for m32 in range(p ** (a2 - a1)):
                m = matrix([[1, 0, 0],
                            [m21, 1, 0],
                            [m31, m32, 1]])
                yield w13 * m

    # w = (123)
    for m21 in range(p ** (a3 - a1)):
        for m23 in range(p ** (a2 - a1 - 1)):
            for m31 in range(p ** (a3 - a2)):
                m = matrix([[1, 0, 0],
                            [m21, 1, p * m23],
                            [m31, 0, 1]])
                yield w123 * m
    # w = (132)
    for m12 in range(p ** (a3 - a2 - 1)):
        for m31 in range(p ** (a2 - a1)):
            for m32 in range(p ** (a3 - a1)):
                m = matrix([[1, p * m12, 0],
                            [0, 1, 0],
                            [m31, m32, 1]])
                yield w132 * m


class HalfIntMatElement(object):

    def __init__(self, T):
        '''
        :params T: half integral matrix of size 3 or a list
        '''
        if isinstance(T, list):
            a, b, c, d, e, f = [ZZ(x) for x in T]
            mat = matrix([[a, f / 2, e / 2],
                          [f / 2, b, d / 2],
                          [e / 2, d / 2, c]])
        else:
            mat = T
        self.__entries = tuple(mat.list())

    def __eq__(self, other):
        if isinstance(other, HalfIntMatElement):
            return self.__entries == other.__entries
        else:
            raise NotImplementedError

    def __repr__(self):
        return self.T.__repr__()

    def __hash__(self):
        return hash(self.__entries)

    @property
    def T(self):
        return matrix(3, self.__entries)

    def right_action(self, g):
        '''
        :param g: matrix of size n
        return self[g] (Siegel's notation)
        '''
        S = g.transpose() * self.T * g
        return HalfIntMatElement(S)

    def satisfy_cong_condition_tp(self, p, alpha):
        '''
        Test if sum_{B mod D} exp(2pi T B D^(-1)) is zero, where D = diag(p^a1, p^a2, a^a3),
        a1, a2, a3 = alpha.
        '''
        return (all(ZZ(self.T[i, i]) % p ** alpha[i] == 0 for i in range(3)) and
                all(ZZ(self.T[i, j] * 2) % p ** alpha[i] == 0
                    for i in range(3) for j in range(i + 1, 3)))

    def is_divisible_by(self, m):
        '''
        Test if self is divisible by m
        :param m: integer
        '''
        return _half_int_mat_is_div_by(self.T, m)

    def __floordiv__(self, other):
        S = matrix(QQ, 3)
        for i in range(3):
            S[i, i] = ZZ(self.T[i, i]) // other
        for i in range(3):
            for j in range(i + 1, 3):
                S[i, j] = S[j, i] = (ZZ(self.T[i, j] * 2) // other) / 2
        return HalfIntMatElement(S)


def alpha_list(dl):
    '''
    Return a list of (a0, a1, a2) with 0 <= a0 <= a1 <= a2 <= dl
    '''
    return [(a0, a1, a2) for a0 in range(dl + 1)
            for a1 in range(a0, dl + 1) for a2 in range(a1, dl + 1)]


def tp_action_fourier_coeff(p, T, F):
    '''
    Return the Tth Fourier coefficient of F|T(p), where F is a modular form.
    :param p: a prime number
    :param T: a half integral matrix or an instance of HalfIntMatElement
    :param F: a dictionary or a Siegel modular form of degree 3
    '''
    p = ZZ(p)
    return _action_fc_base(tp_action_fc_alist(p, T), F, T)


def tp2_action_fourier_coeff(p, i, T, F):
    '''
    Similar to tp_action_fourier_coeff for T_i(p^2).
    '''
    p = ZZ(p)
    return _action_fc_base(tp2_action_fc_alist(p, T, i), F, T)


def _action_fc_base(ls, F, T):
    if not isinstance(T, HalfIntMatElement):
        T = HalfIntMatElement(T)
    res = F[T].zero()
    for s, a, g in ls:
        res = a * F[s].left_action(g) + res
    return res


def hecke_eigenvalue_tp(p, F, T=None):
    '''
    p, F, T: same as aruments of tp_action_fourier_coeff.
    Assuming F is an eigenform, return the eigenvalue for T(p),
    T is used for the computation of Fourier coefficients.
    If T is omitted, T will be set to
    matrix([[1, 1/2, 1/2], [1/2, 1, 1/2], [1/2, 1/2, 1]]).
    '''
    return _hecke_eigenvalue_base(lambda s: tp_action_fourier_coeff(p, s, F), F, T=T)


def hecke_eigenvalue_tp2(p, i, F, T=None):
    '''
    Similar to hecke_eigenvalue_tp for T(p^2).
    '''
    return _hecke_eigenvalue_base(lambda s: tp2_action_fourier_coeff(p, i, s, F), F, T=T)


def spinor_l_euler_factor(p, F, t=None, T=None):
    '''
    F: a dict or Siegel modular form of degree 3.
    Return a polynomial G(t) of degree 8, s.t.
    G(p^(-s))^(-1) is the p-Euler factor of the spinor L function of F.
    '''
    p = ZZ(p)
    if t is None:
        t = PolynomialRing(ZZ, names='t').gen()
    c = {}
    tp = hecke_eigenvalue_tp(p, F, T=T)
    tpp1, tpp2, tpp3 = [hecke_eigenvalue_tp2(p, i, F, T=T) for i in [1, 2, 3]]
    c[0] = ZZ(1)
    c[1] = tp
    c[2] = p * (tpp1 + (p**2 + 1) * tpp2 + (p**2 + 1)**2 * tpp3)
    c[3] = p**3 * tp * (tpp2 + tpp3)
    c[4] = p**6 * (tp**2 * tpp3 + tpp2**2 - 2 * p * tpp1 * tpp3 -
                   2 * (p - 1) * tpp2 * tpp3 -
                   (p**6 + 2 * p**5 + 2 * p**3 + 2 * p - 1) * tpp3**2)
    c[5] = p**6 * tpp3 * c[3]
    c[6] = p**12 * tpp3 ** 2 * c[2]
    c[7] = p**18 * tpp3 ** 3 * c[1]
    c[8] = p**24 * tpp3 ** 4
    return sum((-1)**k * v * t**k for k, v in c.items())


def _hecke_eigenvalue_base(fc_func, F, T=None):
    if T is None:
        T = HalfIntMatElement(matrix([[ZZ(1), ZZ(1) / ZZ(2), ZZ(1) / ZZ(2)],
                                      [ZZ(1) / ZZ(2), ZZ(1), ZZ(1) / ZZ(2)],
                                      [ZZ(1) / ZZ(2), ZZ(1) / ZZ(2), ZZ(1)]]))
    if not isinstance(T, HalfIntMatElement):
        T = HalfIntMatElement(T)
    v1 = fc_func(T).vector
    v = F[T].vector
    if v == 0:
        raise ZeroDivisionError
    else:
        i = next(i for i in range(len(v)) if v[i] != 0)
        return v1[i] / v[i]


@cached_function
def tp_action_fc_alist(p, T):
    '''
    return a list of tuples (S, a, g) s.t.
    S: an instance of HalfIntMatElement
    a: integer
    g: 3 by 3 matrix s.t.
    F|T(p) = sum(a rho(g) F[S] | (a, g, S)).
    '''
    res1 = []
    for alpha in alpha_list(1):
        D = diagonal_matrix([p ** a for a in alpha])
        for V in _gl3_coset_gamma0(alpha, p):
            M = D * V
            S = T.right_action(M.transpose())
            if S.is_divisible_by(p):
                S = S // p
                if S.satisfy_cong_condition_tp(p, alpha):
                    # p**(-6) and p in the third item are for normalization.
                    res1.append(
                        (S, p ** (-6) * mul(p ** alpha[i] for i in range(3) for j in range(i, 3)),
                         M ** (-1) * p))
    return __convert_reduced_nonisom_matrices(res1)


def __convert_reduced_nonisom_matrices(alst):
    red_res = []
    for s, a, g in alst:
        u = _minkowski_reduction_transform_matrix(s.T)
        t = s.right_action(u)
        red_res.append((t, a, g * u.transpose() ** (-1)))

    non_isoms = []

    for s, a, g in red_res:
        q = QuadraticForm(ZZ, 2 * s.T)
        u = None
        for t, _, _ in non_isoms:
            q1 = QuadraticForm(ZZ, 2 * t.T)
            if q.det() == q1.det():
                u = q.is_globally_equivalent_to(q1, return_matrix=True)
                if u:
                    break
        if u:
            non_isoms.append((s.right_action(u), a, g * u.transpose() ** (-1)))
        else:
            non_isoms.append((s, a, g))
    return non_isoms


@cached_function
def tp2_action_fc_alist(p, T, i):
    '''
    similar to tp_action_fc_alist for T_i(p^2) for i = 0, 1, 2, 3.
    '''
    res1 = []

    for alpha in alpha_list(2):
        D = diagonal_matrix([p ** a for a in alpha])
        for V in _gl3_coset_gamma0(alpha, p):
            M = D * V
            S = T.right_action(M.transpose())
            if S.is_divisible_by(p ** 2):
                S = S // (p ** 2)
                res1.append((S, p ** (-12) * _expt_sum(S, p, alpha, i),
                             M ** (-1) * p ** 2))

    return __convert_reduced_nonisom_matrices([(a, b, c) for a, b, c in res1 if b != 0])


def _nearest_integer(x):
    r = floor(x)
    if x - r > 0.5:
        return r + 1
    else:
        return r


def _gaussian_reduction(b1, b2, S):
    '''
    b1, b2: vectors of length 3
    S: symmetric matrix of size 3
    '''
    while True:
        nb1 = b1 * S * b1
        nb2 = b2 * S * b2
        if nb2 < nb1:
            b1, b2 = b2, b1
        x = (b2 * S * b1) / (b1 * S * b1)
        r = _nearest_integer(x)
        a = b2 - r * b1
        if a * S * a >= b2 * S * b2:
            return (b1, b2)
        else:
            b1, b2 = a, b1


def _sym_mat_gen(p, n):
    if n == 1:
        for a in range(p):
            yield matrix([[a]])
    else:
        for s in _sym_mat_gen(p, n - 1):
            ls = [range(p) for _ in range(n)]
            for a in itertools.product(*ls):
                v = matrix([a[:-1]])
                yield block_matrix([[s, v.transpose()], [v, matrix([[a[-1]]])]])


def _gen_gauss_sum_direct_way(N, p, r):
    res = 0
    K = CyclotomicField(p)
    zeta = K.gen()
    for S in _sym_mat_gen(p, N.ncols()):
        if S.change_ring(FiniteField(p)).rank() == r:
            res += zeta ** ((N * S).trace())
    try:
        return QQ(res)
    except TypeError:
        return res


def _generalized_gauss_sum(N, p, r):
    if r == 0:
        return 1
    if p == 2:
        return _gen_gauss_sum_direct_way(N, p, r)
    else:
        N_mp = N.change_ring(FiniteField(p))
        d, _, v = N_mp.smith_form()
        t = d.rank()
        N1 = (v.transpose() * N_mp *
              v).matrix_from_rows_and_columns(range(t), range(t))
        eps = kronecker_symbol(N1.det(), p)
        return _gen_gauss_sum_non_dyadic(p, eps, N.ncols(), t, r)


def _half_int_mat_is_div_by(S, m):
    n = S.ncols()
    return (all(ZZ(S[i, i]) % m == 0 for i in range(n)) and
            all(ZZ(2 * S[i, j]) % m == 0 for i in range(n) for j in range(i + 1, n)))


@cached_function
def _gen_gauss_sum_non_dyadic(p, eps, n, t, r):
    '''
    cf. H. Saito, a generalization of Gauss sums
    '''

    def parenthesis_prod(a, b, m):
        if m == 0:
            return 1
        else:
            return mul(1 - a * b ** i for i in range(m))

    if (n - t) % 2 == 0:
        m = (n - t) // 2
    else:
        m = (n - t + 1) // 2

    if n == r:
        if n % 2 == 1:
            return ((-1) ** ((n - 2 * m + 1) // 2) * p ** ((n ** 2 + (2 * m) ** 2 - 1) // 4) *
                    parenthesis_prod(p ** (-1), p ** (-2), m))
        elif n % 2 == t % 2 == 0:
            return ((-kronecker_symbol(-1, p)) ** ((n - 2 * m) // 2) *
                    eps * p ** ((n ** 2 + (2 * m + 1) ** 2 - 1) // 4) *
                    parenthesis_prod(p ** (-1), p ** (-2), m))
        else:
            return 0
    else:
        diag = [1 for _ in range(t)]
        if eps == -1:
            diag[-1] = least_quadratic_nonresidue(p)
        diag = diag + [0 for _ in range(n - t)]
        N = diagonal_matrix(diag).change_ring(FiniteField(p))
        return _gen_gauss_sum_direct_way(N, p, r)


def _expt_sum(S, p, alpha, i):
    '''
    Return the exponential sum in Miyawaki's paper, where alpha[-1] <= 2, for T_i(p^2).
    '''
    a, b, c = [alpha.count(_i) for _i in range(3)]
    S33 = S.T.matrix_from_rows_and_columns(range(a + b, 3), range(a + b, 3))
    S22 = S.T.matrix_from_rows_and_columns(range(a, a + b), range(a, a + b))
    S32 = S.T.matrix_from_rows_and_columns(range(a + b, 3), range(a))

    if c > 0 and not _half_int_mat_is_div_by(S33, p ** 2):
        return 0
    if c > 0 and b > 0 and any(x % p != 0 for x in (S32 * ZZ(2)).change_ring(ZZ).list()):
        return 0

    if b == 0 and a + c == 3 - i:
        return p ** (c * (c + 1))
    elif b == 0:
        return 0
    else:
        return p ** (c * (c + 1)) * p ** (b * c) * _generalized_gauss_sum(S22, p, b - i)


def _minkowski_reduction(b1, b2, b3, S):

    def inner_prod(x, y):
        return x * S * y

    while True:
        b1, b2, b3 = sorted([b1, b2, b3], key=lambda b: b * S * b)

        b1, b2 = _gaussian_reduction(b1, b2, S)

        b11 = inner_prod(b1, b1)
        b12 = inner_prod(b1, b2)
        b13 = inner_prod(b1, b3)
        b22 = inner_prod(b2, b2)
        b23 = inner_prod(b2, b3)
        b33 = inner_prod(b3, b3)

        y1 = - (b13 / b11 - b12 * b23 / (b11 * b22)) / \
            (1 - b12 ** 2 / (b11 * b22))
        y2 = - (b23 / b22 - b12 * b13 / (b11 * b22)) / \
            (1 - b12 ** 2 / (b11 * b22))

        # Find integers x1, x2 so that norm(b3 + x2 * b2 + x1 * b1) is minimal.
        a_norms_alst = []

        for x1 in [floor(y1), ceil(y1)]:
            for x2 in [floor(y2), ceil(y2)]:
                a = b3 + x2 * b2 + x1 * b1
                a_norms_alst.append((x1, x2, a, inner_prod(a, a)))
        _inner_prod_a = min(x[-1] for x in a_norms_alst)
        x1, x2, a, _ = next(x for x in a_norms_alst if x[-1] == _inner_prod_a)

        if _inner_prod_a >= b33:
            # Change sings of b1, b2, b3 and terminate the alogrithm
            sngs = [sgn(b12), sgn(b13), sgn(b23)]
            bs = [b1, b2, b3]
            try:
                # If b12, b13 or b23 is zero, change sgns of b1, b2, b3 so that
                # b12, b13, b23 >= 0.
                zero_i = sngs.index(0)
                set_ls = [set([1, 2]), set([1, 3]), set([2, 3])]
                t = set_ls[zero_i]
                _other = [x for x in [1, 2, 3] if x not in t][0]
                for x in t:
                    i = set_ls.index(set([x, _other]))
                    if sngs[i] < 0:
                        bs[x - 1] *= -1
                b1, b2, b3 = bs
            except ValueError:
                # Else change sgns so that b12, b13 > 0
                if b12 < 0:
                    b2 = -b2
                if b13 < 0:
                    b3 = -b3
            return (b1, b2, b3)
        else:
            b3 = a


def _minkowski_reduction_transform_matrix(S):
    '''
    Return a unimodular matrix u such that u^t * S * u is reduced in Minkowski's sense.
    '''
    b1, b2, b3 = identity_matrix(QQ, 3).columns()
    c1, c2, c3 = _minkowski_reduction(b1, b2, b3, S)
    return matrix([c1, c2, c3]).transpose()
