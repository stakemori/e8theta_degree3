from os.path import dirname, join
import e8theta_degree3
from sage.rings.all import QQ
from sage.misc.all import cached_function
from sage.all import gcd, latex, ZZ
from sage.structure.factorization import Factorization


@cached_function
def data_dir():
    return join(dirname(e8theta_degree3.__file__), "results", "data")


def half_int_mat_to_list(t):
    return ([t.T[a, a] for a in range(3)] +
            [2 * t.T[i, j] for i, j in [(1, 2), (0, 2), (0, 1)]])


def sort_ts(ts):
    def key(t):
        T = t.T
        return (T.det(), T[0, 0], T[1, 1], T[2, 2], T[0, 1], T[0, 2], T[1, 2])
    return list(sorted(ts, key=key))


def dict_sum(l, dcts):
    kys = list(dcts[0].keys())
    res = {k: 0 for k in kys}
    for k in kys:
        for a, d in zip(l, dcts):
            res[k] += a * d[k]
    return res


def modulo_p(alpha, p, a):
    return QQ(sum(c * a ** i for i, c in enumerate(alpha.list()))) % p


def gcd_of_dict_vals(d):
    return gcd([gcd(d[k].vector) for k in d])


def factorization_normalized(pl):
    '''
    Assume f[0] == 1
    '''
    assert pl[0] == 1
    return Factorization([(a/a[0], b) for a, b in pl.factor()])


def _latex_pol(pl):
    x = pl.parent().gens()[0]

    def _e(v, k):
        if k == 0:
            return ""
        elif v < 0:
            return "-"
        else:
            return "+"

    def _term(v, k):
        if k > 0:
            if v.abs() == 1:
                coeff = ""
            else:
                coeff = latex(v.abs().factor())
            return "%s %s %s" % (_e(v, k), coeff, latex(x**k))
        else:
            return "1"

    def _key(x):
        k = x[0]
        if k in ZZ:
            return k
        else:
            return k[0]

    return " ".join([_term(v, k if k in ZZ else k[0])
                     for k, v in sorted(list(pl.dict().items()),
                                        key=_key)])


def _latex_expt(a, b):
    if b == 1:
        return r"\left(%s\right)" % a
    else:
        return r"\left(%s\right)^{%s}" % (a, b)


def factor_latex(pl):
    '''
    Assume f[0] == 1
    '''
    assert pl[0] == 1
    l = [(a/a[0], b) for a, b in pl.factor()]
    return "".join([_latex_expt(_latex_pol(a), b) for a, b in l])


def _to_allow_break(v):
    return str(v).replace(",", r", \allowbreak")


def _print_dmath(alst, name):
    for a, b in alst:
        print r"\begin{dmath*}"
        print "  b(%s, %s) = %s," % (a, name, _to_allow_break(b))
        print r"\end{dmath*}"
