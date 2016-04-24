# -*- coding: utf-8 -*-
from itertools import groupby
from sage.all import mul
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.functions.other import floor
from sage.matrix.constructor import diagonal_matrix, matrix, block_diagonal_matrix, identity_matrix
from sage.rings.number_field.number_field import CyclotomicField
from sage.misc.cachefunc import cached_function


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
    for m12 in range(p**a):
        yield matrix([[1, m12],
                      [0, 1]])
    for m21 in range(p**(a - 1)):
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
    for m13 in range(p**(a3 - a1 - 1)):
        for m23 in range(p**(a3 - a1 - 1)):
            m = matrix([[1, 0, p * m13],
                        [0, 1, p * m23],
                        [0, 0, 1]])
            yield m

    for m32 in range(p**(a3 - a1)):
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

    for m12 in range(p**(a2 - a1 - 1)):
        for m13 in range(p**(a2 - a1 - 1)):
            m = matrix([[1, p * m12, p * m13],
                        [0, 1, 0],
                        [0, 0, 1]])
            yield m
    for m21 in range(p**(a2 - a1)):
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
    for m12 in range(p**(a2 - a1 - 1)):
        for m13 in range(p**(a3 - a1 - 1)):
            for m23 in range(p**(a3 - a2 - 1)):
                yield matrix([[1, p * m12, p * m13],
                              [0, 1, p * m23],
                              [0, 0, 1]])
    # w = (12)
    for m13 in range(p**(a3 - a2 - 1)):
        for m21 in range(p**(a2 - a1)):
            for m23 in range(p**(a3 - a1 - 1)):
                m = matrix([[1, 0, p * m13],
                            [m21, 1, p * m23],
                            [0, 0, 1]])
                yield w12 * m
    # w = (23)
    for m12 in range(p**(a3 - a1 - 1)):
        for m13 in range(p**(a2 - a1 - 1)):
            for m32 in range(p**(a3 - a2)):
                m = matrix([[1, p * m12, p * m13],
                            [0, 1, 0],
                            [0, m32, 1]])
                yield w23 * m

    # w = (13)
    for m21 in range(p**(a3 - a2)):
        for m31 in range(p**(a3 - a1)):
            for m32 in range(p**(a2 - a1)):
                m = matrix([[1, 0, 0],
                            [m21, 1, 0],
                            [m31, m32, 1]])
                yield w13 * m

    # w = (123)
    for m21 in range(p**(a3 - a1)):
        for m23 in range(p**(a2 - a1 - 1)):
            for m31 in range(p**(a3 - a2)):
                m = matrix([[1, 0, 0],
                            [m21, 1, p * m23],
                            [m31, 0, 1]])
                yield w123 * m
    # w = (132)
    for m12 in range(p**(a3 - a2 - 1)):
        for m31 in range(p**(a2 - a1)):
            for m32 in range(p**(a3 - a1)):
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
        return (all(ZZ(self.T[i, i]) % p**alpha[i] == 0 for i in range(3)) and
                all(ZZ(self.T[i, j] * 2) % p ** alpha[i] == 0
                    for i in range(3) for j in range(i + 1, 3)))

    def is_divisible_by(self, m):
        '''
        Test if self is divisible by m
        :param m: integer
        '''
        return (all(ZZ(self.T[i, i]) % m == 0 for i in range(3)) and
                all(ZZ(self.T[i, j] * 2) % m == 0
                    for i in range(3) for j in range(i + 1, 3)))

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
    res = 0
    p = ZZ(p)
    if not isinstance(T, HalfIntMatElement):
        T = HalfIntMatElement(T)


@cached_function
def __tp_action_fc_alist(p, T):
    res1 = []
    for alpha in alpha_list(1):
        D = diagonal_matrix([p**a for a in alpha])
        for V in _gl3_coset_gamma0(alpha, p):
            M = D * V
            S = T.right_action(M.transpose())
            if S.is_divisible_by(p):
                S = S // p
                if S.satisfy_cong_condition_tp(p, alpha):
                    # p**(-6) and p in the third item are for normalization.
                    res1.append(
                        (S, p**(-6) * mul(p**alpha[i] for i in range(3) for j in range(i, 3)),
                         M**(-1) * p))
    res = []
    # Use reduced quadratic forms
    for s, a, g in res1:
        u = _minkowski_reduction_transform_matrix(s.T)
        t = s.right_action(u)
        res.append((t, a, g * u.transpose()))
    return res


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
        if a * S * a >= b1 * S * b1:
            return (b1, b2)
        else:
            b1, b2 = a, b1


def _B_genrator(p, alpha):
    '''
    Yield B s.t. BD^(-1) symmetric and BD^(-1) mod Sym_3(Z)
    '''
    a1, a2, a3 = alpha
    for b11 in range(p**a1):
        for b22 in range(p**a2):
            for b33 in range(p**a3):
                for b21 in range(p**a1):
                    for b31 in range(p**a1):
                        for b32 in range(p**a2):
                            yield matrix([[b11, b21 * p**(a2 - a1), b31 * p**(a3 - a1)],
                                          [b21, b22, b32 * p**(a3 - a2)],
                                          [b31, b32, b33]])


def _expt_sum(S, p, alpha, D, pred=None):
    '''
    Return the exponential sum in Miyawaki's paper, where alpha[-1] <= 2.
    '''
    p = ZZ(p)
    K = CyclotomicField(p**2)
    a = K.gens()[0]
    res = 0
    for B in _B_genrator(p, alpha):
        if pred is None or pred(B, D, p):
            res += a ** ZZ((S.T * B * D**(-1)).trace() * p ** 2)
    try:
        return QQ(res)
    except TypeError:
        return res


def __minkowski_reduction(b1, b2, b3, S):

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
        x1 = _nearest_integer(y1)
        x2 = _nearest_integer(y2)

        a = b3 + x2 * b2 + x1 * b1
        if inner_prod(a, a) >= b33:
            return (b1, b2, b3)
        else:
            b3 = a


def _minkowski_reduction_transform_matrix(S):
    '''
    Return a unimodular matrix u such that u * S * u^t is reduced in Minkowski's sense.
    '''
    b1, b2, b3 = identity_matrix(QQ, 3).columns()
    c1, c2, c3 = __minkowski_reduction(b1, b2, b3, S)
    return matrix([c1, c2, c3]).transpose()
