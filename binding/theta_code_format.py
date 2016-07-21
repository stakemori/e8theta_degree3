import itertools

from e8theta_degree3.binding.code_gen import (FmpzStyle,
                                              pol_to_fmpz_codes_and_result_var)
from e8theta_degree3.gl3_repn import gl3_repn_module, matrix_var
from sage.arith.all import lcm
from sage.matrix.all import matrix
from sage.misc.all import cached_function, mul
from sage.modules.all import vector
from sage.rings.all import QQ, PolynomialRing


@cached_function
def _s_t_u_ring(base_ring=None):
    if base_ring is None:
        base_ring = QQ
    R = PolynomialRing(
        base_ring, names=("s0, s1, s2, s3, s4, s5, s6, s7,"
                          "t0, t1, t2, t3, t4, t5, t6, t7,"
                          "u0, u1, u2, u3, u4, u5, u6, u7"))
    return R


def _conj(pol):
    return pol.map_coefficients(lambda c: c.conjugate())


def _rl_part(pol):
    return pol.map_coefficients(lambda x: x.list()[0]).change_ring(QQ)


def _im_part(pol):
    return pol.map_coefficients(lambda x: x.list()[1]).change_ring(QQ)


@cached_function
def _pol_basis_factor_dct_and_ls(wt):
    '''
    wt: a list/tuple of non-increasing integers of length 3.
    Let M a corresponding repn of GL3.
    return a pair (d, l)
    l: a list of polynomials obtained from prime factors of M.basis_as_pol().
    d: dict s.t whose keys are M.basis() and
    f => [const, (b1, t1), (b2, t2), ... ],
    where f = const b1^t1 * b2^t2 * ..., and b1, b2, ... in l.
    '''
    M = gl3_repn_module(wt)
    basis = M.basis()
    facs = [(pol, pol.factor()) for pol in basis]
    l = list(
        set(itertools.chain(*[[a for a, _ in fc.pols] for _, fc in facs])))
    return [{pol: [fc.const] + fc.pols for pol, fc in facs}, l]


@cached_function
def euclidean_basis():
    basis_vecs = [(QQ(1) / QQ(2), QQ(1) / QQ(2), QQ(1) / QQ(2),
                   QQ(1) / QQ(2), QQ(1) / QQ(2), QQ(1) / QQ(2),
                   QQ(1) / QQ(2), QQ(1) / QQ(2)),
                  (0, 1, 0, 0, 0, 0, 0, 1),
                  (0, 0, 1, 0, 0, 0, 0, 1),
                  (0, 0, 0, 1, 0, 0, 0, 1),
                  (0, 0, 0, 0, 1, 0, 0, 1),
                  (0, 0, 0, 0, 0, 1, 0, 1),
                  (0, 0, 0, 0, 0, 0, 1, 1),
                  (0, 0, 0, 0, 0, 0, 0, 2)]
    return [vector([QQ(a) for a in v1]) for v1 in basis_vecs]


def to_eulidian_vec(t):
    return sum([a * b for a, b in zip(t, euclidean_basis())])


def _bideterminant_prime_factors_dict(mat, wt):
    '''
    wt: a list/tuple of non-increasing integers of length 3
    mat: 3 * 8 matrix with mat * mat.transpose() = 0 with coefficients in
    an imaginary quadratic field.
    Return a dict
    a: (real_part, imag_part) as polynomials of _s_t_u_ring(QQ),
    where a is a prime factor of polynomial basis
    '''
    R = _s_t_u_ring()
    stu_mt = matrix(R, 3, R.gens())
    stu_mt = matrix([to_eulidian_vec(v) for v in stu_mt.rows()]).transpose()
    subs_dct = dict(zip(matrix_var().list(), (mat * stu_mt).list()))
    _, l = _pol_basis_factor_dct_and_ls(wt)
    d = {a: a.subs(subs_dct) for a in l}
    return {a: (_rl_part(b), _im_part(b)) for a, b in d.items()}


def _pol_basis_as_polof_factors(wt, imag_quad, names_base=('real_part', 'imag_part')):
    '''
    wt: weight of repn of GL3
    imag_quad: imaginary quadratic field
    Returns a pair of dicts (d, subs_dct).
    subs_dct is a dict s.t.
    a => (rli, imi)
    where a is a prime factor of baisis_as_pol,
    rli, imi are variables and omega is the gen of imag_quad.
    The set of keys of d is basis of the corresponding representation.
    Its value at x is a pair (f, g) of polynomials of rl0, im0, rl1, im1, ... s.t.
    x.subs({a => rli + omega * imi}) = f + g * omega.
    '''
    d, l = _pol_basis_factor_dct_and_ls(wt)
    n = len(l)
    names = itertools.chain(
        *[[names_base[0] + str(a), names_base[1] + str(a)] for a in range(n)])
    R = PolynomialRing(imag_quad, names=list(names))
    omega = imag_quad.gen()
    r_gens = R.gens()
    subs_dct = {fc: gns[0] + omega * gns[1]
                for gns, fc in zip([r_gens[a:a + 2] for a in range(0, 2 * n, 2)], l)}

    def _subs(ls):
        return ls[0] * mul(subs_dct[a] ** b for a, b in ls[1:])

    res = {k: _subs(v) for k, v in d.items()}
    return ({k: (_rl_part(v), _im_part(v)) for k, v in res.items()},
            {k: (_rl_part(v), _im_part(v)) for k, v in subs_dct.items()})


def _init_code(variables, res_str):
    indent = "  "
    res = ";\n".join([indent + "fmpz_t %s" % v for v in variables]) + ";\n"
    res = res + "\n"
    res = res + ";\n".join([indent + "fmpz_init(%s)" % v for v in variables]) + ";\n"
    res = res + "\n" + indent + "char * {res_str};".format(res_str=res_str)

    return res


def _bi_det_factors_code(rl_alst, im_alst, denom_lcm, var_dct):
    indent = " " * 26
    lines_format = ("{indent}/* Computation of {var}" +
                    " = {part} part of {lcm} * ({bi_det}) */\n{indent}{code}")

    def _code(alst, idx, part):
        return (";\n\n").join([lines_format.format(indent=indent, var=var_dct[k][idx],
                                                   bi_det=k,
                                                   lcm=denom_lcm,
                                                   code=(";\n" + indent).join(codes),
                                                   part=part) for k, codes in alst])
    rl_code = _code(rl_alst, 0, "real")
    im_code = _code(im_alst, 1, "imaginary")
    return rl_code + ";\n\n\n" + im_code + ";"


def _coeffs_code(pol_code_alst, res_vars):
    indent = " " * 26
    comment_line = "{indent}/* Computation of {var} = {pol} */\n{indent}"

    code_blocks = [(comment_line.format(var=v, pol=str(pl), indent=indent) +
                    (";\n" + indent).join(codes) + ";")
                   for (pl, codes), v in zip(pol_code_alst, res_vars)]
    return "\n\n".join(code_blocks)


def _str_concat_code(res_vars, res_str):
    indent = " " * 2
    _format_str = '"' + ",".join(['%s' for _ in range(len(res_vars))]) + '"'
    get_str_code = ", ".join(["fmpz_get_str(NULL, 10, %s)" % (v, ) for v in res_vars])

    _code = '''
int buf_size = asprintf(&{res_str}, {_format_str}, {get_str_code});

if (buf_size == -1) {{
  exit(1);
}}
    '''.format(get_str_code=get_str_code, _format_str=_format_str, res_str=res_str)
    return indent + _code.replace("\n", "\n" + indent)


def _cleanup_code(variables):
    indent = "  "
    return ";\n".join(indent + "fmpz_clear(%s)" % (v,) for v in variables) + ";"


def code_format(func_name, wt, mat, real_part=True):
    '''
    wt: non-increasing list/tuple of non-negative integers of length 3.
    mat: 3 * 8 matrix with mat * mat.transpose() = 0 with coefficients in
    an imaginary quadratic field.
    real_part: boolian. If False, the imaginary part of the Fourier coefficient
    will be computed.
    Here for alpha = a + b * omega in K (an imaginary quadratic field with the
    generator omega), the real part and imaginary part of alpha are
    a and b respectively.
    '''
    tmp_var_name = "a"
    sum_tmp_var_name = "tmp"
    res_str_name = "res_str"
    sty = FmpzStyle()

    bdt_facs = _bideterminant_prime_factors_dict(mat, wt)
    _facs_pols = itertools.chain(*([a, b] for a, b in bdt_facs.values()))
    _facs_pols_lcm = lcm([b.denominator()
                          for b in itertools.chain(*(a.dict().values() for a in _facs_pols))])
    # Remove denominators
    bdt_facs = {k: (a * _facs_pols_lcm, b * _facs_pols_lcm) for k, (a, b) in bdt_facs.items()}

    Vrho = gl3_repn_module(wt)
    bs_pl_dct, bdt_var_dct = _pol_basis_as_polof_factors(wt, mat.base_ring())
    coef_pol_pairs = [bs_pl_dct[b] for b in Vrho.basis()]

    if real_part:
        coef_pols = [a for a, _ in coef_pol_pairs]
    else:
        coef_pols = [a for _, a in coef_pol_pairs]

    def _facs_pols_code_and_vars(i):
        return {k: pol_to_fmpz_codes_and_result_var(v[i], tmp_var_name, bdt_var_dct[k][i])
                for k, v in bdt_facs.items()}

    facs_pols_code_rl_dct = _facs_pols_code_and_vars(0)
    facs_pols_code_im_dct = _facs_pols_code_and_vars(1)

    res_vars = ["res" + str(i) for i, _ in enumerate(coef_pols)]

    coefs_pol_code_alst = [(pl, pol_to_fmpz_codes_and_result_var(pl, tmp_var_name,
                                                                 sum_tmp_var_name))
                           for pl in coef_pols]

    tmp_vars = list(set(list(itertools.chain(*(l for _, l in facs_pols_code_rl_dct.values()))) +
                        list(itertools.chain(*(l for _, l in facs_pols_code_im_dct.values()))) +
                        list(itertools.chain(*(l for _, (_, l) in coefs_pol_code_alst)))))
    tmp_vars = sorted(tmp_vars)
    tmp_vars.append(sum_tmp_var_name)

    def _key_fun(x):
        return str(bdt_var_dct[x[0]])

    facs_pols_code_rl_alst = sorted([(k, codes) for k, (codes, _) in facs_pols_code_rl_dct.items()],
                                    key=_key_fun)
    facs_pols_code_im_alst = sorted([(k, codes) for k, (codes, _) in facs_pols_code_im_dct.items()],
                                    key=_key_fun)

    # Remove tmp vars and add summing
    coefs_pol_code_alst1 = [(pl, codes + [sty.add_z(v, v, sum_tmp_var_name)])
                            for (pl, (codes, _)), v in zip(coefs_pol_code_alst, res_vars)]

    _vrs = (tmp_vars + res_vars + sorted([str(r) for r, _ in bdt_var_dct.values()]) +
            sorted([str(im) for _, im in bdt_var_dct.values()]) +
            list(itertools.chain(*[[s + str(i) for i in range(8)] for s in ["s", "t", "u"]])))

    init_code_str = _init_code(_vrs, res_str_name)

    bi_det_factors_code_str = _bi_det_factors_code(facs_pols_code_rl_alst,
                                                   facs_pols_code_im_alst, _facs_pols_lcm,
                                                   bdt_var_dct)
    coeffs_code_str = _coeffs_code(coefs_pol_code_alst1, res_vars)

    str_concat_code_str = _str_concat_code(res_vars, res_str_name)

    cleanup_code_str = _cleanup_code(_vrs)

    header = '''#define _GNU_SOURCE             /* for asprintf */
#include <stdio.h>
#include "e8vectors.h"

inline int inner_prod(int s[8], int t[8])
{
  return ((2*s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7]) * t[0] +
          (s[0] + 2*s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + 2*s[7]) * t[1] +
          (s[0] + s[1] + 2*s[2] + s[3] + s[4] + s[5] + s[6] + 2*s[7]) * t[2] +
          (s[0] + s[1] + s[2] + 2*s[3] + s[4] + s[5] + s[6] + 2*s[7]) * t[3] +
          (s[0] + s[1] + s[2] + s[3] + 2*s[4] + s[5] + s[6] + 2*s[7]) * t[4] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + 2*s[5] + s[6] + 2*s[7]) * t[5] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + 2*s[6] + 2*s[7]) * t[6] +
          (s[0] + 2*s[1] + 2*s[2] + 2*s[3] + 2*s[4] + 2*s[5] + 2*s[6] + 4*s[7]) * t[7]);
}

'''

    code = '''{header}

char * {func_name}(int a, int b, int c, int d, int e, int f)
{{
  cache_vectors();
  /* Use static to avoid segmentation fault */
  static int vs1[MAX_NM_OF_VECTORS][8];
  static int vs2[MAX_NM_OF_VECTORS][8];
  static int vs3[MAX_NM_OF_VECTORS][8];

  _set_vs3(vs1, vs2, vs3, a, b, c);

{init_code}

  for (int i = 0; i < num_of_vectors[a]; i++)
    {{
      for (int j = 0; j < num_of_vectors[b]; j++)
        {{
          for (int k = 0; k < num_of_vectors[c]; k++)
            {{
              if (inner_prod(vs1[i], vs2[j]) == f)
                {{
                  if (inner_prod(vs1[i], vs3[k]) == e)
                    {{
                      if (inner_prod(vs2[j], vs3[k]) == d)
                        {{
                          fmpz_set_si(s0, vs1[i][0]);
                          fmpz_set_si(s1, vs1[i][1]);
                          fmpz_set_si(s2, vs1[i][2]);
                          fmpz_set_si(s3, vs1[i][3]);
                          fmpz_set_si(s4, vs1[i][4]);
                          fmpz_set_si(s5, vs1[i][5]);
                          fmpz_set_si(s6, vs1[i][6]);
                          fmpz_set_si(s7, vs1[i][7]);

                          fmpz_set_si(t0, vs2[j][0]);
                          fmpz_set_si(t1, vs2[j][1]);
                          fmpz_set_si(t2, vs2[j][2]);
                          fmpz_set_si(t3, vs2[j][3]);
                          fmpz_set_si(t4, vs2[j][4]);
                          fmpz_set_si(t5, vs2[j][5]);
                          fmpz_set_si(t6, vs2[j][6]);
                          fmpz_set_si(t7, vs2[j][7]);

                          fmpz_set_si(u0, vs3[k][0]);
                          fmpz_set_si(u1, vs3[k][1]);
                          fmpz_set_si(u2, vs3[k][2]);
                          fmpz_set_si(u3, vs3[k][3]);
                          fmpz_set_si(u4, vs3[k][4]);
                          fmpz_set_si(u5, vs3[k][5]);
                          fmpz_set_si(u6, vs3[k][6]);
                          fmpz_set_si(u7, vs3[k][7]);

{bi_det_factors_code}

{coeffs_code}

                        }}
                    }}
                }}
            }}
        }}
    }}

{str_concat_code}
{cleanup_code}

  /* {res_str_name} must be freed later. */
  return {res_str_name};
}}
    '''.format(init_code=init_code_str,
               bi_det_factors_code=bi_det_factors_code_str,
               coeffs_code=coeffs_code_str,
               str_concat_code=str_concat_code_str,
               cleanup_code=cleanup_code_str,
               res_str_name=res_str_name,
               header=header, func_name=func_name)
    return code