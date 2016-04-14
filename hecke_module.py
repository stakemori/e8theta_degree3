# -*- coding: utf-8 -*-

from sage.all import ZZ, mul, matrix, block_diagonal_matrix
from itertools import groupby


def _index_of_gamma_0_gl_n(alphas, p):
    '''
    Returns delta(a1, ..., an) defined in Shimura, Euler products and Eisenstein
    series, pp 118, (15.1.7).
    '''
    if p in ZZ:
        p = ZZ(p)

    def _bn(n):
        return mul(1 - p ** (-i) for i in xrange(1, n + 1))

    e_r_ls = [(k, len(list(v)))
              for k, v in groupby(sorted(alphas), lambda x: x)]
    res = _bn(len(alphas)) / mul(_bn(r) for _, r in e_r_ls)
    for i, (ei, ri) in enumerate(e_r_ls):
        for j, (ej, rj) in enumerate(e_r_ls):
            if i < j:
                res *= p ** ((ej - ei) * ri * rj)
    return res


def _gl2_coset_gamma_subscript0(a, p):
    w = matrix([[0, -1],
                [1, 0]])
    for m21 in range(p**a):
        yield matrix([[1, 0],
                      [m21, 1]])
    for m12 in range(p**(a - 1)):
        m = matrix([[1, p * m12],
                    [0, 1]])
        yield w * m


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


def _gl3_coset_gamma0(alphas, p):
    r'''
    Let alphas = [a0, a1, a2] with a0 <= a1 <= a2,
    g = diag([p^a0, p^a1, p^a2]), and Gamma0 = g^(-1) GL3(Z) g ∧ GL3(Z).
    Return a complete set Gamma0 \ GL3(Z).
    '''
    if p in ZZ:
        p = ZZ(p)
    a0, a1, a2 = alphas
    if a0 < a1 < a2:
        return list(__gl3_coset_gamma0_distinct(a0, a1, a2, p))
    elif a0 == a1 and a1 < a2:
        return list(__gl3_coset_gamma0_2_1(a0, a2, p))
    elif a0 < a1 and a1 == a2:
        return list(__gl3_coset_gamma0_1_2(a0, a2, p))
    else:
        raise ValueError


def __gl3_coset_gamma0_2_1(a1, a3, p):
    w23 = matrix([[1, 0, 0],
                  [0, 0, 1],
                  [0, 1, 0]])
    for m13 in range(p**(a3 - a1 - 1)):
        for m23 in range(p**(a3 - a1 - 1)):
            m = matrix([[1, 0, m13],
                        [0, 1, m23],
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
    for m23 in range(p**(a2 - a1)):
        m = matrix([[1, 0, 0],
                    [0, 1, m23],
                    [0, 0, 1]])
        for g in _gl2_coset_gamma_subscript0(a2 - a1, p):
            n = block_diagonal_matrix(g, matrix([[1]]))
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
                            [m13, m12, 1]])
                yield w132 * m
