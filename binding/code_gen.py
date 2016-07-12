# -*- coding: utf-8 -*-
from sage.all import PolynomialRing, mul, ZZ, uniq, SR, factor, is_prime_power
from abc import ABCMeta, abstractmethod, abstractproperty


class CodeStyle(object):

    __metaclass__ = ABCMeta

    @abstractproperty
    def sep(self):
        pass

    @abstractmethod
    def zero_z(self, a):
        pass

    @abstractmethod
    def set_mul(self, a, b, c):
        pass

    @abstractmethod
    def set_add_mul(self, a, b, c):
        pass

    @abstractmethod
    def set_si(self, a, b):
        pass

    @abstractmethod
    def set_z(self, a, b):
        pass

    @abstractmethod
    def add_ui(self, a, b, c):
        pass

    @abstractmethod
    def sub_ui(self, a, b, c):
        pass

    @abstractmethod
    def pow_ui(self, a, b, c):
        pass

    @abstractmethod
    def add_mul_ui(self, a, b, c):
        pass

    @abstractmethod
    def add_z(self, a, b, c):
        pass

    @abstractmethod
    def sub_mul_ui(self, a, b, c):
        pass

    @abstractmethod
    def sub_z(self, a, b, c):
        pass

    @abstractmethod
    def mul_ui(self, a, b, c):
        pass

    @abstractmethod
    def neg(self, a, b):
        pass

    @abstractmethod
    def mul_2exp(self, a, b, e):
        pass

    def mul_si(self, a, b, c):
        c = ZZ(c)
        res_dct_special_case = {0: self.zero_z(a),
                                1: self.set_z(a, b),
                                -1: self.neg(a, b)}
        if c in res_dct_special_case:
            return res_dct_special_case[c]

        tpw = False
        e = None
        fc = factor(c.abs())
        if is_prime_power(c.abs()) and fc[0][0] == 2:
            # The case when c is a power of 2
            tpw = True
            e = fc[0][1]

        if c > 0:
            return self._mul_si_pos(a, b, c, e, tpw)
        else:
            return "%s%s %s" % (self._mul_si_pos(a, b, -c, e, tpw),
                                self.sep,
                                self.neg(a, a))

    def _mul_si_pos(self, a, b, c, e, tpw):
        '''
        Assume c > 1.
        '''
        if tpw:
            return self.mul_2exp(a, b, e)
        else:
            return self.mul_ui(a, b, c)


class FmpzStyle(CodeStyle):

    @property
    def sep(self):
        return ";"

    def zero_z(self, a):
        return "fmpz_zero(%s)" % a

    def set_mul(self, a, b, c):
        return "fmpz_mul(%s, %s, %s)" % (a, b, c)

    def set_add_mul(self, a, b, c):
        return "fmpz_addmul(%s, %s, %s)" % (a, b, c)

    def set_si(self, a, b):
        return "fmpz_set_si(%s, %s)" % (a, b)

    def set_z(self, a, b):
        return "fmpz_set(%s, %s)" % (a, b)

    def add_ui(self, a, b, c):
        return "fmpz_add_ui(%s, %s, %s)" % (a, b, c)

    def sub_ui(self, a, b, c):
        return "fmpz_sub_ui(%s, %s, %s)" % (a, b, c)

    def pow_ui(self, a, b, c):
        return "fmpz_pow_ui(%s, %s, %s)" % (a, b, c)

    def add_mul_ui(self, a, b, c):
        return "fmpz_addmul_ui(%s, %s, %s)" % (a, b, c)

    def add_z(self, a, b, c):
        return "fmpz_add(%s, %s, %s)" % (a, b, c)

    def sub_mul_ui(self, a, b, c):
        return "fmpz_submul_ui(%s, %s, %s)" % (a, b, c)

    def sub_z(self, a, b, c):
        return "fmpz_sub(%s, %s, %s)" % (a, b, c)

    def neg(self, a, b):
        return "fmpz_neg(%s, %s)" % (a, b)

    def mul_2exp(self, a, b, e):
        return "fmpz_mul_2exp(%s, %s, %s)" % (a, b, e)

    def mul_ui(self, a, b, c):
        return "fmpz_mul_ui(%s, %s, %s)" % (a, b, c)


class PythonStyle(CodeStyle):

    @property
    def sep(self):
        return ";"

    def zero_z(self, a):
        return "%s = 0" % a

    def set_mul(self, a, b, c):
        return "%s = %s * %s" % (a, b, c)

    def set_add_mul(self, a, b, c):
        return "%s = %s + %s * %s" % (a, a, b, c)

    def set_si(self, a, b):
        return "%s = %s" % (a, b)

    def set_z(self, a, b):
        return self.set_si(a, b)

    def add_ui(self, a, b, c):
        return "%s = %s + %s" % (a, b, c)

    def sub_ui(self, a, b, c):
        return "%s = %s - %s" % (a, b, c)

    def pow_ui(self, a, b, c):
        return "%s = %s ** %s" % (a, b, c)

    def add_mul_ui(self, a, b, c):
        return self.set_add_mul(a, b, c)

    def add_z(self, a, b, c):
        return self.add_ui(a, b, c)

    def sub_mul_ui(self, a, b, c):
        return "%s = %s - %s * %s" % (a, a, b, c)

    def sub_z(self, a, b, c):
        return "%s = %s - %s" % (a, b, c)

    def mul_ui(self, a, b, c):
        return "%s = %s * (%s)" % (a, b, c)

    def neg(self, a, b):
        return "%s = -%s" % (a, b)

    def mul_2exp(self, a, b, e):
        return "%s = %s * 2**(%s)" % (a, b, e)


class Mul(object):

    def __init__(self, l, r, ngens):
        self.l = l
        self.r = r
        self.ngens = ngens

    def __repr__(self):
        return "Mul(%s, %s)" % (self.l, self.r)

    def codes1(self, tmp_var, sty):
        cds = self.l.codes1(tmp_var, sty)
        cds.append(sty.set_mul(tmp_var, tmp_var, self.r))
        return cds

    def codes(self, tmp_vars, sty):
        '''
        Return list of codes, result tmp var and used tmp vars.
        '''
        if self.ngens == 1:
            return (self.codes1(tmp_vars[0], sty), tmp_vars[0], [tmp_vars[0]])
        cd, v, vrs = self.l.codes(tmp_vars, sty)
        cd.append(sty.set_mul(v, v, self.r))
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

    def codes1(self, tmp_var, sty):
        m = Mul(self.b, self.c, 1)
        cds = m.codes1(tmp_var, sty)
        a = ZZ(self.a.pl)
        if a > 0:
            cds.append(sty.add_ui(tmp_var, tmp_var, a))
        elif a < 0:
            cds.append(sty.sub_ui(tmp_var, tmp_var, -a))
        else:
            raise ValueError("a must be non zero")
        return cds

    def codes(self, tmp_vars, sty):
        if self.ngens == 1:
            return (self.codes1(tmp_vars[0], sty), tmp_vars[0], [tmp_vars[0]])
        cds, v, vrs = self.b.codes(tmp_vars, sty)
        cds1, v1, vrs1 = self.a.codes([a for a in tmp_vars if a != v], sty)
        cds.extend(cds1)
        cds.append(sty.set_add_mul(v1, v, self.c))
        vrs.extend(vrs1)
        vrs.extend([v, v1])
        return (cds, v1, vrs)


class Deg1Pol(object):

    def __init__(self, pl):
        self.pl = pl

    def __repr__(self):
        return "Deg1Pol(%s)" % (self.pl, )

    def codes(self, tmp_vars, sty):
        v = tmp_vars[0]
        pl = self.pl
        gens = pl.parent().gens()

        if pl.constant_coefficient() != 0:
            codes = [sty.set_si(v, self.pl.constant_coefficient())]
        else:
            _ls = list(pl.dict().items())
            _t, a = _ls[0]
            codes = [sty.mul_si(v, _expt(_t, gens), a)]
            pl = pl.parent()(sum(_expt(_t, gens) * a for _t, a in _ls[1:]))
        if pl.parent().ngens() == 1:
            x = pl.parent().gen()
            cd = _admul_code(v, x, pl[1], sty)
            if cd:
                codes.append(cd)
        else:
            gens = pl.parent().gens()
            for k, a in pl.dict().iteritems():
                if sum(k) == 1:
                    x = _expt(k, gens)
                    cd = _admul_code(v, x, a, sty)
                    if cd:
                        codes.append(cd)
        return (codes, v, [v])

    def codes1(self, tmp_var, sty):
        return self.codes([tmp_var], sty)[0]


def _admul_code(v, x, a, sty):
    '''
    If a is 0, it return None
    '''
    res = None
    if a > 1:
        res = sty.add_mul_ui(v, x, a)
    elif a == 1:
        res = sty.add_z(v, v, x)
    elif a == -1:
        res = sty.sub_z(v, v, x)
    elif a < -1:
        res = sty.sub_mul_ui(v, x, -a)
    return res


def _count(ts, i):
    return sum(1 for t in ts if t[i] > 0)


def _expt(t, vrs):
    if t in ZZ:
        return vrs[0] ** t
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


def _pol_to_codes_and_res_var_pow(pl, name, res_var_name, sty):
    v = name + "0"
    if pl.constant_coefficient() == 0:
        codes = [sty.zero_z(res_var_name)]
    else:
        codes = [sty.set_si(res_var_name, pl.constant_coefficient())]
    vrs = [v]
    v1 = None
    if pl.parent().ngens() == 1:
        x = pl.parent().gen()
        for e, cf in pl.dict().items():
            if e > 1:
                codes.append(sty.pow_ui(v, x, e))
                codes.append(_admul_code(res_var_name, v, cf, sty))
            elif e == 1:
                codes.append(_admul_code(res_var_name, x, cf, sty))
    else:
        gns = pl.parent().gens()
        for t, cf in pl.dict().items():
            if sum(t) > 1:
                if any(a > 1 for a in t):
                    v1 = name + "1"
                codes.extend(_monom_codes(t, v, v1, gns, sty))
                codes.append(_admul_code(res_var_name, v, cf, sty))
            elif sum(t) == 1:
                codes.append(_admul_code(res_var_name, _expt(t, gns), cf, sty))
        if v1 is not None:
            vrs.append(v1)
    return (codes, uniq(vrs))


def pol_to_fmpz_codes_and_result_var(pl, name, res_var_name, sty=None, algorithm=None):
    '''
    pl: polynomial
    name: name for tmp vars
    res_var_name: variable name for storing result
    Return (codes, tmp_var_names),
    where codes is a list of strings for computing pl and
    tmp_var_names is a list of used temp variable names.
    '''
    pl = pl.change_ring(ZZ)
    if sty is None:
        sty = FmpzStyle()

    # If pl is equal to a single variable, just set the result
    if pl.degree() == 1 and pl.dict().values()[0] == 1:
        return ([sty.set_z(res_var_name, str(pl))], [])

    if algorithm == "horner":
        n = pl.parent().ngens()
        vrs = [name + str(a) for a in range(n)]
        e = _to_expr(pl)
        codes, v, vrs = e.codes(vrs, sty)
        codes.append(sty.set_z(res_var_name, v))
        return (codes, uniq(vrs))
    elif algorithm is None:
        return _pol_to_codes_and_res_var_pow(pl, name, res_var_name, sty)
    else:
        raise ValueError


def pol_factor_to_code_and_result_var(pl, name, res_var_name, sty=None, algorithm=None):
    '''
    Similar to pol_to_fmpz_codes_and_result_var. But factor pl before computing code.
    '''
    pl = pl.change_ring(ZZ)
    if sty is None:
        sty = FmpzStyle()
    pl = pl.change_ring(ZZ)
    R = pl.parent()
    factors = [(f, e) for f, e in pl.factor() if R(f).degree() > 0]
    if len(factors) == 1:
        return pol_to_fmpz_codes_and_result_var(pl, name, res_var_name,
                                                sty=sty, algorithm=algorithm)

    _pl = mul(f**e for f, e in factors)
    if pl.parent().ngens() == 1:
        const_cf = ZZ(pl.leading_coefficient() / _pl.leading_coefficient())
    else:
        const_cf = ZZ(pl.lc() / _pl.lc())
    assert const_cf * _pl == pl

    # Compute tmp vars
    tvars = []
    factor_codes_ls = []
    for f, e in factors:
        _cds, _tvars = pol_to_fmpz_codes_and_result_var(
            f, name, "{res}", sty=sty, algorithm=algorithm)
        tvars.extend(_tvars)
        factor_codes_ls.append(_cds)
    tvars = uniq(tvars)
    tvar0 = name + str(len(tvars))
    assert tvar0 not in tvars
    tvars.append(tvar0)

    def _replace_res_in_codes(v, codes):
        return [c.format(res=v) for c in codes]

    codes = []
    if all(e == 1 for _, e in factors):  # just multiplication
        codes.extend(_replace_res_in_codes(res_var_name, factor_codes_ls[0]))
        for factor_codes in factor_codes_ls[1:]:
            codes.extend(_replace_res_in_codes(tvar0, factor_codes))
            codes.append(sty.set_mul(res_var_name, res_var_name, tvar0))
    else:
        fcodes0 = factor_codes_ls[0]
        e0 = factors[0][1]
        # set res = f0**e0
        codes.extend(_replace_res_in_codes(res_var_name, fcodes0))
        if e0 > 1:
            codes.append(sty.pow_ui(res_var_name, res_var_name, e0))

        for (_, e), fcodes in zip(factors[1:], factor_codes_ls[1:]):
            # set tvar0 = f**e
            codes.extend(_replace_res_in_codes(tvar0, fcodes))
            if e > 1:
                codes.append(sty.pow_ui(tvar0, tvar0, e))
            # set res = res * tvar0
            codes.append(sty.set_mul(res_var_name, res_var_name, tvar0))
    if const_cf != 1:
        codes.append(sty.mul_si(res_var_name, res_var_name, const_cf))
    return (codes, tvars)


def fmpz_init_clear(vrs, sep=" "):
    if isinstance(vrs, str):
        vrs = SR.var(vrs)
    res = sep.join(
        ["fmpz_t %s;" % v for v in vrs] + ["fmpz_init(%s);" % v for v in vrs])
    res = res + "\n\n"
    res = res + sep.join("fmpz_clear(%s);" % v for v in vrs)
    return res


def _monom_codes(t, v, v1, gns, sty):
    '''t: monomial, v: variable name for result.
    v1: tmp var name.
    Return codes.
    '''
    codes = []
    _tv = [(e, x) for e, x in zip(t, gns) if e > 0]
    if all(e == 1 for e, x in _tv):
        codes.append(sty.set_z(v, _tv[0][1]))
        for _, x in _tv[1:]:
            codes.append(sty.set_mul(v, v, x))
    else:
        _tv = list(sorted(_tv, key=lambda x: -x[1]))
        e0, x0 = _tv[0]
        codes.append(sty.pow_ui(v, x0, e0))
        for e, x in _tv[1:]:
            if e == 1:
                codes.append(sty.set_mul(v, v, x))
            else:
                codes.append(sty.pow_ui(v1, x, e))
                codes.append(sty.set_mul(v, v, v1))
    return codes


def _test():
    gi = globals()["get_ipython"]
    ip = gi()

    def check(R):
        for _ in range(100):
            if R.ngens() == 1:
                while True:
                    f = R.random_element(degree=(-1, 10))
                    h = R.random_element(degree=(-1, 10))
                    g = R.random_element(degree=(-1, 10))
                    if all(a != 0 for a in [f, g, h]):
                        break
            else:
                while True:
                    f = R.random_element(degree=10)
                    g = R.random_element(degree=10)
                    h = R.random_element(degree=10)
                    if all(a != 0 for a in [f, g, h]):
                        break

            for alg in [None, 'horner']:
                print alg
                c = "; ".join(pol_to_fmpz_codes_and_result_var(
                    f, "a", "res", sty=PythonStyle(), algorithm=alg)[0])
                ip.run_cell(c)
                res = globals()["res"]
                assert res == f

            for alg in [None, 'horner']:
                print alg, "factor"
                for f in [f, f * g, f * g * h, f**2, f * g**2, f**2 * g**2 * h]:
                    for g in [f, ZZ(10) * f]:
                        c = "; ".join(pol_factor_to_code_and_result_var(
                            g, "a", "res", sty=PythonStyle(), algorithm=alg)[0])
                        ip.run_cell(c)
                        res = globals()["res"]
                        assert res == g, str(g)

    R1 = PolynomialRing(ZZ, names="x")
    R2 = PolynomialRing(ZZ, names="x, y")
    R3 = PolynomialRing(ZZ, names="x, y, z")
    print "check1"
    ip.run_cell("x = PolynomialRing(ZZ, names='x').gen()")
    check(R1)
    print "check2"
    ip.run_cell("x, y = PolynomialRing(ZZ, names='x, y').gens()")
    check(R2)
    print "check3"
    ip.run_cell("x, y, z = PolynomialRing(ZZ, names='x, y, z').gens()")
    check(R3)
