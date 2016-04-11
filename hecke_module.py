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
