from sage.all import PolynomialRing, QQ, matrix


def schur_polynomial(n, wt):
    R = PolynomialRing(QQ, names=["x%s" % i for i in range(n)])
    l = [i + a for i, a in zip(reversed(range(n)), wt)]
    return R(_vandermonde(n, l) / _vandermonde(n, reversed(range(n))))


def _vandermonde(n, wt):
    R = PolynomialRing(QQ, names=["x%s" % i for i in range(n)])
    vrs = R.gens()
    m = matrix([[x ** a for x in vrs] for a in wt])
    return m.det()
