# -*- coding: utf-8 -*-

from sage.all import ZZ, mul
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


def _gl3_coset_gamma0(alphas, p):
    r'''
    Let alphas = [a0, a1, a2] with a0 <= a1 <= a2,
    g = diag([p^a0, p^a1, p^a2]), and Gamma0 = g^(-1) GL3(Z) g ∧ GL3(Z).
    Return a complete set Gamma0 \ GL3(Z).
    '''
    a0, a1, a2 = alphas
    if a0 < a1 < a2:
        return __gl3_coset_gamma0_distinct(a0, a1, a2, p)
    elif a0 == a1 and a1 < a2:
        pass
    elif a0 < a1 and a1 == a2:
        pass
    else:
        raise ValueError


def __gl3_coset_gamma0_distinct(a0, a1, a2, p):
    pass
