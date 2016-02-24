from sage.all import PolynomialRing, QQ, matrix, Combinations


def schur_polynomial(n, wt):
    R = PolynomialRing(QQ, names=["x%s" % i for i in xrange(n)])
    l = [i + a for i, a in zip(reversed(xrange(n)), wt)]
    return R(_vandermonde(n, l) / _vandermonde(n, reversed(range(n))))


class YoungTableu(object):

    def __init__(self, n=None, numbers=None):
        self._n = n
        self._numbers = numbers

    @property
    def n(self):
        return self._n

    @property
    def numbers(self):
        return self._numbers

    def weight(self):
        return [len(l) for l in self.numbers]

    def __repr__(self):
        return "\n".join(str(a) for a in self.numbers)


def _increasing_nums(n, m, lower_bds=None):
    '''
    n, m: positive integers,
    lower_bds: list of integers of length m.
    Returns a generator of m numbers (a0, a1, ..., a_(m-1)),
    such that a0, a1, .. in [1, ..., n] and
    a0 >= lower_bds[0], a1 >= lower_bds[1] .. and a_(m-1) >= lower_bds[m-1].
    '''
    cmbs = Combinations(range(1, n + 1), m)
    cmbs = (list(sorted(a)) for a in cmbs)
    if lower_bds is None:
        for a in cmbs:
            yield a
    else:
        for a in cmbs:
            if all(x >= y for x, y in zip(a, lower_bds)):
                yield a


def semistandard_young_tableaux(n, wt):
    '''Returns a list of semistandard Young tableaux.
    '''
    pass


def _vandermonde(n, wt):
    R = PolynomialRing(QQ, names=["x%s" % i for i in range(n)])
    vrs = R.gens()
    m = matrix([[x ** a for x in vrs] for a in wt])
    return m.det()
