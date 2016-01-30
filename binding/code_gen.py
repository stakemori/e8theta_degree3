# -*- coding: utf-8; mode: sage -*-
from sage.all import PolynomialRing, mul, ZZ, uniq


class Add(object):

    def __init__(self, l, r):
        self.l = l
        self.r = r

    def __repr__(self):
        return "Add(%s, %s)" % (self.l, self.r)


def set_mul(a, b, c):
    return "fmpz_mul(%s, %s, %s)" % (a, b, c)


def set_add_mul(a, b, c):
    return "fmpz_addmul(%s, %s, %s)" % (a, b, c)


class Mul(object):

    def __init__(self, l, r, ngens):
        self.l = l
        self.r = r
        self.ngens = ngens

    def __repr__(self):
        return "Mul(%s, %s)" % (self.l, self.r)

    def codes1(self, tmp_var):
        cds = self.l.codes1(tmp_var)
        cds.append(set_mul(tmp_var, tmp_var, self.r))
        return cds

    def codes(self, tmp_vars):
        '''
        Return list of codes, result tmp var and used tmp vars.
        '''
        if self.ngens == 1:
            return (self.codes1(tmp_vars[0]), tmp_vars[0], [tmp_vars[0]])
        cd, v, vrs = self.l.codes(tmp_vars)
        cd.append(set_mul(v, v, self.r))
        vrs.append(v)
        return (cd, v, vrs)


class AddMul(object):

    def __init__(self, a, b, c, ngens):
        # a + b * c
        self.a = a
        self.b = b
        self.c = c
        self.ngens = ngens

    def __repr__(self):
        return "AddMul(%s, %s, %s)" % (self.a, self.b, self.c)

    def codes1(self, tmp_var):
        m = Mul(self.b, self.c, 1)
        cds = m.codes1(tmp_var)
        a = ZZ(self.a.pl)
        if a > 0:
            cds.append("fmpz_add_ui(%s, %s, %s)" % (tmp_var, tmp_var, a))
        elif a < 0:
            cds.append("fmpz_sub_ui(%s, %s, %s)" % (tmp_var, tmp_var, -a))
        return cds

    def codes(self, tmp_vars):
        if self.ngens == 1:
            return (self.codes1(tmp_vars[0]), tmp_vars[0], [tmp_vars[0]])
        cds, v, vrs = self.b.codes(tmp_vars)
        cds1, v1, vrs1 = self.a.codes([a for a in tmp_vars if a != v])
        cds.extend(cds1)
        cds.append(set_add_mul(v1, v, self.c))
        vrs.extend(vrs1)
        vrs.extend([v, v1])
        return (cds, v1, vrs)


class Deg1Pol(object):

    def __init__(self, pl):
        self.pl = pl

    def __repr__(self):
        return "Deg1Pol(%s)" % (self.pl, )

    def codes(self, tmp_vars):
        v = tmp_vars[0]
        if self.pl.constant_coefficient() != 0:
            codes = ["fmpz_set_si(%s, %s)" %
                     (v, self.pl.constant_coefficient())]
        else:
            codes = ["fmpz_zero(%s)" % (v, )]
        pl = self.pl
        if pl.parent().ngens() == 1:
            x = pl.parent().gen()
            cd = _admul_code(v, x, pl[1])
            if cd:
                codes.append(cd)
        else:
            gens = pl.parent().gens()
            for k, a in pl.dict().iteritems():
                if sum(k) == 1:
                    x = _expt(k, gens)
                    cd = _admul_code(v, x, a)
                    if cd:
                        codes.append(cd)
        return (codes, v, [v])

    def codes1(self, tmp_var):
        return self.codes([tmp_var])[0]


def _admul_code(v, x, a):
    res = None
    if a > 1:
        res = "fmpz_addmul_ui(%s, %s, %s)" % (v, x, a)
    elif a == 1:
        res = "fmpz_add(%s, %s, %s)" % (v, v, x)
    elif a == -1:
        res = "fmpz_sub(%s, %s, %s)" % (v, v, x)
    elif a < -1:
        res = "fmpz_submul_ui(%s, %s, %s)" % (v, x, -a)
    return res


def _count(ts, i):
    return sum(1 for t in ts if t[i] > 0)


def _expt(t, vrs):
    return mul(x ** e for e, x in zip(t, vrs))


def _to_expr1(pl):
    A = pl.parent()
    x = A.gen()
    if pl.degree() <= 1:
        return Deg1Pol(pl)
    elif pl[0] == 0:
        return Mul(_to_expr1(A(pl / x)), x, 1)
    else:
        return AddMul(Deg1Pol(A(pl[0])), _to_expr1(A((pl - pl[0]) / x)), x, 1)


def _to_expr(pl):
    if pl.degree() <= 1:
        return Deg1Pol(pl)
    A = pl.parent()
    n = A.ngens()
    if n == 1:
        return _to_expr1(pl)
    kys = [tuple(t) for t, _ in pl.dict().iteritems()]
    idx = list(sorted([(_count(kys, i), i)
                       for i in xrange(n)], key=lambda x: x[0]))[-1][1]
    _gns = A.gens()
    x = _gns[idx]
    B = PolynomialRing(A.base_ring(), names=[str(a) for a in _gns if a != x])
    f1 = A(sum(_expt(t, _gns) * a for t,
               a in pl.dict().iteritems() if t[idx] > 0) / x)
    f2 = B(sum(_expt(t, _gns) * a for t,
               a in pl.dict().iteritems() if t[idx] == 0))
    if f2 == 0:
        return Mul(_to_expr(f1), x, n)
    else:
        return AddMul(_to_expr(f2), _to_expr(f1), x, n)


def _expr_to_pol(expr):
    if isinstance(expr, Deg1Pol):
        return expr.pl
    elif isinstance(expr, AddMul):
        return _expr_to_pol(expr.a) + _expr_to_pol(expr.b) * _expr_to_pol(expr.c)
    elif isinstance(expr, Mul):
        return _expr_to_pol(expr.l) * _expr_to_pol(expr.r)
    else:
        return expr


def pol_to_fmpz_code_and_result_var_horner(pl, name):
    n = pl.parent().ngens()
    vrs = [name + str(a) for a in range(n)]
    e = _to_expr(pl)
    codes, v, vrs = e.codes(vrs)
    return ("\n".join(a + ";" for a in codes), v, uniq(vrs))


def pol_to_fmpz_code_and_result_var(pl, name, res_var_name, algorithm=None):
    if algorithm == "horner":
        n = pl.parent().ngens()
        vrs = [name + str(a) for a in range(n)]
        e = _to_expr(pl)
        codes, v, vrs = e.codes(vrs)
        codes.append("fmpz_set(%s, %s)" % (res_var_name, v))
    else:
        v = name + "0"
        if pl.constant_coefficient() == 0:
            codes = ["fmpz_zero(%s)" % (res_var_name, )]
        else:
            codes = ["fmpz_set_si(%s, %s)" %
                     (res_var_name, pl.constant_coefficient())]
        vrs = [v]
        if pl.parent().ngens() == 1:
            x = pl.parent().gen()
            for e, cf in pl.dict().items():
                if e > 1:
                    codes.append("fmpz_pow_ui(%s, %s, %s)" % (v, x, e))
                    codes.append(_admul_code(res_var_name, v, cf))
                elif e == 1:
                    codes.append(_admul_code(res_var_name, x, cf))
        else:
            gns = pl.parent().gens()
            for t, cf in pl.dict().items():
                if sum(t) > 1:
                    codes.append(_monom_codes(t, v, gns))
                    codes.append(_admul_code(res_var_name, v, cf))
                elif sum(t) == 1:
                    codes.append(_admul_code(res_var_name, _expt(t, gns), cf))

        return ("\n".join(a + ";" for a in codes), uniq(vrs))


def _monom_codes(t, v, gns):
    '''t: monomial, v: variable name.
    Return codes.
    '''
    _tv = [(e, x) for e, x in zip(t, v) if e > 0]
    if all(e == 1 for e, _ in _tv):
        pass
