import itertools
import hashlib
import os
import re

from e8theta_degree3.binding.code_gen import (MpirStyle,
                                              pol_to_mpz_codes_and_result_var)
from e8theta_degree3.gl3_repn import gl3_repn_module, matrix_var
from sage.arith.all import lcm
from sage.matrix.all import matrix
from sage.misc.all import cached_function, mul
from sage.rings.all import QQ, PolynomialRing


class CCode(object):

    def __init__(self, code, s_codes=None):
        self.code = code
        self.s_codes = s_codes


class SepCCode(object):

    def __init__(self, func_name, code, args):
        self.func_name = func_name
        self.code = code
        self.args = args


def save_code_to_file(directory, fname_base, func_names, wt, mats, real_parts=None,
                      overwrite=False, sty=None, num_of_procs=None, is_sparse_mat=False,
                      separate_code=False):
    '''
    fname_base: string
    This save code for theta series to fname.h and fname.c
    If separate_code, return file names for separated code, otherwise return None.
    '''
    hf = header_file_format(fname_base, func_names)
    code_obj = code_format(func_names, wt, mats, real_parts=real_parts,
                           sty=sty, num_of_procs=num_of_procs, is_sparse_mat=is_sparse_mat,
                           separate_code=separate_code)
    fnameh = os.path.join(directory, fname_base + ".h")
    fnamec = os.path.join(directory, fname_base + ".c")
    if (not overwrite) and (os.path.exists(fnameh) or os.path.exists(fnamec)):
        raise IOError("file aleadly exists.")
    with open(fnameh, "w") as fp:
        fp.write(hf)

    with open(fnamec, "w") as fp:
        fp.write(code_obj.code)
    if separate_code:
        sep_code_fs = []
        for s_c in code_obj.s_codes:
            fname_scode = os.path.join(directory, s_c.func_name + '.c')
            sep_code_fs.append(fname_scode)
            with open(fname_scode, 'w') as fp:
                fp.write(s_c.code)
        return sep_code_fs
    else:
        return


def generate_cython_and_build_scripts(directory,
                                      cython_fname_base,
                                      c_fname_base,
                                      func_names,
                                      c_func_names, wt, mats, real_parts=None,
                                      overwrite=False, sty=None, num_of_procs=1,
                                      is_sparse_mat=False, separate_code=False):
    '''
    directory(string): must be a subdirectory of binding
    Generate c source, cython source ,build scripts and a makefile in directory.
    If is_sparse_mat is true, then every (i, j)th element of matrices given should be zero
    if j > 5 (in case when 3 * 8 matrix) j > 6 (in case when 3 * 16 matrix).
    Can compile by "make compile-cython" in that directory if
    required libs were compiled by "make compile-theta_vectors".
    '''
    vec_len = mats[0].ncols()
    for m in mats:
        assert m * m.transpose() == 0
    if is_sparse_mat:
        for m in mats:
            assert all(m[i, j] == 0 for i in range(3)
                       for j in range(6 if vec_len == 8 else 7, vec_len))

    c_fname = c_fname_base + "_c"
    _cython_code = _cython_format(c_fname, c_func_names, func_names, num_of_procs)
    _setup_py_code = _setup_py_format(cython_fname_base, c_fname)
    sep_code_fs_maybe = save_code_to_file(directory, c_fname, c_func_names, wt, mats,
                                          real_parts=real_parts, overwrite=overwrite, sty=sty,
                                          num_of_procs=num_of_procs, is_sparse_mat=is_sparse_mat,
                                          separate_code=separate_code)
    _makefile_code = _makefile_format(c_fname, sep_code_fs_maybe)

    def _fname(f):
        return os.path.join(directory, f)

    cython_file = _fname(cython_fname_base + ".pyx")
    setup_file = _fname("setup.py")
    make_file = _fname("Makefile")
    _init_file = _fname("__init__.py")
    file_code_alst = [(cython_file, _cython_code),
                      (setup_file, _setup_py_code),
                      (make_file, _makefile_code),
                      (_init_file, "\n")]
    if (not overwrite) and any(os.path.exists(a) for a, _ in file_code_alst):
        raise IOError("file already exists.")
    for f, c in file_code_alst:
        with open(f, "w") as fp:
            fp.write(c)


def _setup_py_format(cython_src_file, c_lib_name):
    fmt = '''from os.path import abspath, curdir, join, pardir
from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize

setup(
    name="{cython_src_file_upp}",
    ext_modules=cythonize(
        Extension("{cython_src_file}",
                  sources=["{cython_src_file}.pyx"],
                  include_dirs=[abspath(curdir), abspath(pardir)],
                  library_dirs=[join(abspath(pardir), "lib")],
                  libraries=["{c_lib_name}", "e8vectors", "rank16_vectors"]),
    ),
)
'''.format(cython_src_file=cython_src_file,
           cython_src_file_upp=cython_src_file.upper(),
           c_lib_name=c_lib_name)
    return fmt


def _makefile_format(c_src_file, sep_code_fnames):
    if sep_code_fnames:
        code = " ".join(os.path.basename(f) for f in sep_code_fnames)
    else:
        code = ""
    _fmt = '''current_dir = $(shell pwd)
parent_dir = $(shell dirname "$(current_dir)")
DEBUGOPT = -Wall -g -Og -std=c11
PATHOPT = -L$(parent_dir)/lib -I$(current_dir) -I$(parent_dir)
LIBOPTBASE = -lmpir
OPT = -O3 -std=c11 -Wall -Wextra
SHARED = -shared -fPIC
CC = gcc

compile-c-lib:
\t$(CC) {c_src_file}.c {code} -o $(parent_dir)/lib/lib{c_src_file}.so $(PATHOPT) $(OPT) \\
\t$(LIBOPTBASE) $(SHARED)

compile-cython: compile-c-lib
\tsage -c 'sh.eval("python setup.py build_ext -i")'
'''.format(c_src_file=c_src_file, code=code)
    return _fmt


def _cython_format(c_header_file, c_func_names, cython_func_names, num_of_procs):
    ext_fmt = "    cpdef char * {c_func_name}(int, int, int, int, int, int, int)"
    ext_code = "\n".join([ext_fmt.format(c_func_name=c_func_name) for c_func_name in c_func_names])
    ext_code = ('cdef extern from "{c_header_file}.h":\n'.format(
        c_header_file=c_header_file) + ext_code)
    funcs_code = "\n".join([_cython_format_each(c_fcn, cy_fcn, num_of_procs)
                            for c_fcn, cy_fcn in zip(c_func_names, cython_func_names)])

    code = '''import itertools
from multiprocessing import Pool
from sage.rings.all import Integer, QQ, ZZ
from sage.misc.all import cached_function
from sage.modules.all import vector
from sage.matrix.all import MatrixSpace
from libc.stdlib cimport free
from e8theta_degree3.gl3_repn import GL3RepnElement
include "cysignals/signals.pxi"

{ext_code}

{funcs_code}
'''.format(ext_code=ext_code, funcs_code=funcs_code)
    return code


def _cython_format_each(c_func_name, cython_func_name, num_of_procs):
    if num_of_procs == 1:
        theta_body = 'return _{cfn}_part((0, m))\n'.format(cfn=cython_func_name)
    else:
        theta_body = '''p = Pool(processes={npcs})
    try:
        res = sum(p.map(_{cfn}_part, zip(range({npcs}), itertools.repeat(m, {npcs}))))
    except KeyboardInterrupt:
        p.terminate()
        p.join()
    finally:
        p.close()
        p.join()
    return res'''.format(cfn=cython_func_name, npcs=num_of_procs)
    _fmt = '''
def _{cfn}_part(i_red_m):
    i_red, m = i_red_m

    if not (m in MatrixSpace(QQ, 3) and (2 * m in MatrixSpace(ZZ, 3)) and
            (m[a, a] in ZZ for a in range(3)) and m.transpose() == m):
        raise ValueError("m must be a half integral matrix of size 3.")
    l = [m[i, i] for i in range(3)] + [2*m[t] for t in [(1, 2), (0, 2), (0, 1)]]
    a, b, c, d, e, f = [int(x) for x in l]
    if max([a, b, c]) > 7:
        raise ValueError("Diagonal elements are too large.")
    sig_on()
    cdef char* c_str = {c_func_name}(i_red, a, b, c, d, e, f)
    cdef bytes py_str
    try:
        py_str = c_str
    finally:
        py_strs = py_str.split(",")
        # If py_str is empty, asprintf was not called
        if len(py_strs) > 1:
            free(c_str)
    sig_off()
    assert len(py_strs) > 1, "MAX_NM_REPRS_RK16 is too small."
    res = [Integer(a) for a in py_strs]
    normalizing_num = ZZ(res[0])
    return vector(res[1:]) / normalizing_num


@cached_function
def {cfn}(m):
    {theta_body}
'''.format(c_func_name=c_func_name, cfn=cython_func_name, theta_body=theta_body)
    return _fmt


@cached_function
def _s_t_u_ring(base_ring=None, vec_len=8):
    if base_ring is None:
        base_ring = QQ
    R = PolynomialRing(
        base_ring, names=(["s" + str(i) for i in range(vec_len)] +
                          ["t" + str(i) for i in range(vec_len)] +
                          ["u" + str(i) for i in range(vec_len)]))
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


def _bideterminant_prime_factors_dict(mat, wt):
    '''
    wt: a list/tuple of non-increasing integers of length 3
    mat: 3 * 8 or (3 * 16) matrix with mat * mat.transpose() = 0 with coefficients in
    an imaginary quadratic field.
    Return a dict
    a: (real_part, imag_part) as polynomials of _s_t_u_ring(QQ),
    where a is a prime factor of polynomial basis
    '''
    vec_len = mat.ncols()
    R = _s_t_u_ring(vec_len=vec_len, base_ring=mat.base_ring())
    stu_mt = matrix(R, 3, R.gens()).transpose()
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


def _init_code(variables, res_str, vec_len,
               normalizing_num, normalizing_num_var, wt_sum,
               sty=None, is_sparse_mat=False):
    indent = "  "
    res = ";\n".join([indent + sty.z_type_decl(v) for v in variables]) + ";\n"
    res = res + "\n"
    res = res + ";\n".join([indent + sty.init(str(v)) for v in variables]) + ";\n"
    res = res + "\n" + indent + "char * {res_str};".format(res_str=res_str) + "\n"
    res = res + indent + sty.set_str(normalizing_num_var,
                                     '"%s"' % (normalizing_num, ), "10") + ";\n"
    res = res + indent + sty.mul_2exp(normalizing_num_var,
                                      normalizing_num_var, wt_sum) + ";\n"

    res = res + _vec_k_init_code(is_sparse_mat, vec_len)
    return res


def _vec_k_init_code(is_sparse_mat, vec_len):
    if is_sparse_mat and vec_len == 16:
        res = '''
  static int reprs_k[MAX_NM_REPRS_RK16][16];
  static int reprs_j[MAX_NM_REPRS_RK16][16];
  static unsigned int num_of_classes_k[MAX_NM_REPRS_RK16];
  static unsigned int num_of_classes_j[MAX_NM_REPRS_RK16];

  int num_of_reprs_k = repr_modulo_autom_rk16(c, reprs_k, num_of_classes_k);


  if (num_of_reprs_k + 1 > MAX_NM_REPRS_RK16)
    {
      return "";
    }

  static int vecs_j[MAX_NM_OF_VECTORS_RK16][16];

  static int reprs_i[MAX_NM_REPRS_RK16][16];
  static unsigned int num_of_classes_i[MAX_NM_REPRS_RK16];
  static int vecs_i[MAX_NM_OF_VECTORS_RK16][16];
'''
    elif (not is_sparse_mat) and vec_len == 16:
        raise NotImplementedError
    elif (not is_sparse_mat) and vec_len == 8:
        res = '''
  int * vec_k = cached_vectors_ptr[c];
  vec_k += i_red * 8;
'''
    elif is_sparse_mat and vec_len == 8:
        res = '''

  static int reprs_k[MAX_NM_REPRS][8];
  static unsigned int num_of_classes_k[MAX_NM_REPRS];

  static int reprs_j[MAX_NM_REPRS][8];
  static unsigned int num_of_classes_j[MAX_NM_REPRS];
  static int vecs_j[MAX_NM_OF_VECTORS][8];
  int num_of_vecs_j;

  static int reprs_i[MAX_NM_REPRS][8];
  static unsigned int num_of_classes_i[MAX_NM_REPRS];
  static int vecs_i[MAX_NM_OF_VECTORS][8];

  int num_of_reprs_k = repr_modulo_autom(c, reprs_k, num_of_classes_k);
  int num_of_reprs_j;
'''
    return res


def _bi_det_factors_code(rl_alst, im_alst, denom_lcm, var_dct, num_spaces):
    indent = " " * num_spaces
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


def _coeffs_code(pol_code_alst, res_vars, num_spaces):
    indent = " " * num_spaces
    comment_line = "{indent}/* Computation of {var} = {pol} */\n{indent}"

    code_blocks = [(comment_line.format(var=v, pol=str(pl), indent=indent) +
                    (";\n" + indent).join(codes) + ";")
                   for (pl, codes), v in zip(pol_code_alst, res_vars)]
    return "\n\n".join(code_blocks)


def _str_concat_code(res_vars, res_str, sty=None):
    indent = " " * 2
    _format_str = '"' + ",".join(['%s' for _ in range(len(res_vars))]) + '"'
    get_str_code = ", ".join([sty.get_str("NULL", "10", str(v)) for v in res_vars])

    _code = '''int buf_size = asprintf(&{res_str}, {_format_str}, {get_str_code});

if (buf_size == -1) {{
  exit(1);
}}
'''.format(get_str_code=get_str_code, _format_str=_format_str, res_str=res_str)
    _code = _code.replace("\n", "\n" + indent)
    regexp = re.compile(r"^\s*$")
    _code = indent + _code
    return "\n".join(["" if regexp.match(c) else c for c in _code.split("\n")])


def _cleanup_code(variables, sty=None):
    indent = "  "
    return ";\n".join(indent + sty.clear(str(v)) for v in variables) + ";"


def header_file_format(fname, func_names):
    decl = '''
char * {func_name}(int i_red, int a, int b, int c, int d, int e, int f);
'''
    funcs_declaration = "\n\n".join([decl.format(func_name=a) for a in func_names])

    res = '''#ifndef _{filename_name_upp_base}_H_
#define _{filename_name_upp_base}_H_

{funcs_declaration}

#endif /* _{filename_name_upp_base}_H_ */
'''.format(funcs_declaration=funcs_declaration,
           filename_name_upp_base=os.path.basename(fname.upper()))
    return res


def _inner_prod_code(vec_len):
    code8 = '''inline int inner_prod(int s[8], int t[8])
{
  return (s[0]*t[0] + s[1]*t[1] + s[2]*t[2] + s[3]*t[3] +
          s[4]*t[4] + s[5]*t[5] + s[6]*t[6] + s[7]*t[7]) >> 2;
}
'''
    code16 = '''int inner_prod(int s[16], int t[16])
{
  return (s[0]*t[0] + s[1]*t[1] + s[2]*t[2] + s[3]*t[3] + s[4]*t[4] + s[5]*t[5] +
          s[6]*t[6] + s[7]*t[7] + s[8]*t[8] + s[9]*t[9] + s[10]*t[10] + s[11]*t[11] +
          s[12]*t[12] + s[13]*t[13] + s[14]*t[14] + s[15]*t[15]) >> 2;
}
'''
    if vec_len == 8:
        return code8
    else:
        return code16


def _set_s_code(vec_len, set_si_func, vecs_dict, num_spaces):
    indent = " " * num_spaces

    def _set_s_code_each(s_name, vecs):
        return "\n".join([indent +
                          "{set_si_func}({s_name}{n}, {vecs}[{n}]);".format(
                              set_si_func=set_si_func, vecs=vecs, n=str(n),
                              s_name=s_name)
                          for n in range(vec_len)])

    return "\n\n".join(_set_s_code_each(s, vecs_dict[i])
                       for s, i in zip(["s", "t", "u"], ["i", "j", "k"]))


def code_format_header_innerprod(vec_len):
    code = '''#define _GNU_SOURCE             /* for asprintf */
#include <stdio.h>
#include "{header_file}.h"
#include <stdlib.h>
#include <mpir.h>
#include "memory.h"
#include "vector_utils.h"

{inner_prod}
'''.format(header_file="e8vectors" if vec_len == 8 else "rank16_vectors",
           inner_prod=_inner_prod_code(vec_len))
    return code


def _vec_j_normalize_code(vec_len, is_sparse_mat):
    if vec_len == 8 and is_sparse_mat:
        return '''
      int wo_sign_indices_array[8][16] = {{0}};
      int w_sign_indices[8] = {0};
      set_w_sign_indices(w_sign_indices, reprs_k[k], 8, 2);
      set_wo_sign_indices_array(wo_sign_indices_array, reprs_k[k], 8, 2);

      num_of_vecs_j = 0;
      int * cached_vec_b = cached_vectors_ptr[b];
      for (int j = 0; j < num_of_vectors[b]; j++, cached_vec_b += 8)
        {
          if (inner_prod(cached_vec_b, reprs_k[k]) == d)
            {
              memcpy(vecs_j[num_of_vecs_j++], cached_vec_b, sizeof(int) * 8);
            }
        }
      num_of_reprs_j = repr_modulo_autom_w_indices(vecs_j, num_of_vecs_j, reprs_j,
                                                   num_of_classes_j,
                                                   w_sign_indices, wo_sign_indices_array);
'''
    elif vec_len == 16 and is_sparse_mat:
        return '''
      int wo_sign_indices_array[8][16] = {{0}};
      int w_sign_indices[16] = {0};
      set_w_sign_indices(w_sign_indices, reprs_k[k], 16, 9);
      set_wo_sign_indices_array(wo_sign_indices_array, reprs_k[k], 16, 9);
      int num_of_vecs_j = 0;

      int * cached_vec_b = cached_vectors_rk16_ptr[b];

      for (int j = 0; j < num_of_vectors_rk16[b]; j++, cached_vec_b += 16)
        {
          if (inner_prod(cached_vec_b, reprs_k[k]) == d)
            {
              memcpy(vecs_j[num_of_vecs_j++], cached_vec_b, sizeof(int) * 16);
            }
        }
      int num_of_reprs_j = repr_modulo_autom_rk16_w_indices(vecs_j, num_of_vecs_j, reprs_j,
                                                            num_of_classes_j,
                                                            w_sign_indices, wo_sign_indices_array);

      if (num_of_reprs_j + 1 > MAX_NM_REPRS_RK16)
        {
          return "";
        }

'''
    elif vec_len == 8 and (not is_sparse_mat):
        return '''
      int * vec_j = cached_vectors_ptr[b];
'''
    elif vec_len == 16 and (not is_sparse_mat):
        raise NotImplementedError


def _vec_i_normalize_code(vec_len, is_sparse_mat):
    if vec_len == 8 and is_sparse_mat:
        return '''
          int wo_sign_indices_array[8][16] = {{0}};
          int w_sign_indices[8] = {0};
          set_wo_sign_indices_array2(wo_sign_indices_array, reprs_j[j], reprs_k[k], 8, 2);
          set_w_sign_indices_2(w_sign_indices, reprs_j[j], reprs_k[k], 8, 2);

          int num_of_vecs_i = 0;
          int * cached_vec_a = cached_vectors_ptr[a];
          for (int l = 0; l < num_of_vectors[a]; l++, cached_vec_a += 8)
            {
              if ((inner_prod(cached_vec_a, reprs_k[k]) == e) &&
                  (inner_prod(cached_vec_a, reprs_j[j]) == f))
                {
                  memcpy(vecs_i[num_of_vecs_i++], cached_vec_a, sizeof(int) * 8);
                }
            }
          int num_of_reprs_i = repr_modulo_autom_w_indices(vecs_i, num_of_vecs_i, reprs_i,
                                                           num_of_classes_i,
                                                           w_sign_indices,
                                                           wo_sign_indices_array);
'''
    elif vec_len == 16 and is_sparse_mat:
        return '''
          int wo_sign_indices_array[8][16] = {{0}};
          int w_sign_indices[16] = {0};
          set_wo_sign_indices_array2(wo_sign_indices_array, reprs_j[j], reprs_k[k], 16, 9);
          set_w_sign_indices_2(w_sign_indices, reprs_j[j], reprs_k[k], 16, 9);

          int num_of_vecs_i = 0;
          int * cached_vec_a = cached_vectors_rk16_ptr[a];
          for (int l = 0; l < num_of_vectors_rk16[a]; l++, cached_vec_a += 16)
            {
              if ((inner_prod(cached_vec_a, reprs_k[k]) == e) &&
                  (inner_prod(cached_vec_a, reprs_j[j]) == f))
                {
                  memcpy(vecs_i[num_of_vecs_i++], cached_vec_a, sizeof(int) * 16);
                }
            }

          int num_of_reprs_i = repr_modulo_autom_rk16_w_indices(vecs_i, num_of_vecs_i, reprs_i,
                                                                num_of_classes_i,
                                                                w_sign_indices,
                                                                wo_sign_indices_array);


          if (num_of_reprs_i + 1 > MAX_NM_REPRS_RK16)
            {
              return "";
            }

'''
    elif vec_len == 8 and (not is_sparse_mat):
        return '''
          int * vec_i = cached_vectors_ptr[a];
'''
    elif vec_len == 16 and (not is_sparse_mat):
        raise NotImplementedError


def code_format(func_names, wt, mats, real_parts=None,
                factor_pol=False, sty=None,
                num_of_procs=1, is_sparse_mat=False, separate_code=False):
    n = len(mats)
    if real_parts is None:
        real_parts = [True for _ in range(n)]
    code_objs = [code_format_theta(func_name, wt, mat,
                                   real_part=real_part,
                                   factor_pol=factor_pol,
                                   sty=sty, num_of_procs=num_of_procs,
                                   is_sparse_mat=is_sparse_mat, separate_code=separate_code)
                 for func_name, mat, real_part in zip(func_names, mats, real_parts)]
    header = code_format_header_innerprod(mats[0].ncols()) + "\n"
    if separate_code:
        for c in code_objs:
            for s_c in c.s_codes:
                header = header + _sep_code_signature(s_c) + "\n"

    code = (header + "\n" + "\n\n".join((c.code for c in code_objs)))
    if separate_code:
        return CCode(code, sum([c.s_codes for c in code_objs], []))
    else:
        return CCode(code, None)


def _sep_code_signature(s_code):
    code = "void {func_name}({args_code});".format(
        func_name=s_code.func_name,
        args_code=", ".join(("mpz_t " + a for a in s_code.args)))
    return code


def code_format_theta(func_name, wt, mat, real_part=True, factor_pol=False, sty=None,
                      num_of_procs=1, is_sparse_mat=False, separate_code=False):
    '''
    wt: non-increasing list/tuple of non-negative integers of length 3.
    mat: 3 * 8 (or 3 * 16) matrix with mat * mat.transpose() = 0 with coefficients in
    an imaginary quadratic field.
    real_part: boolian. If False, the imaginary part of the Fourier coefficient
    will be computed.
    Here for alpha = a + b * omega in K (an imaginary quadratic field with the
    generator omega), the real part and imaginary part of alpha are
    a and b respectively.
    Return an instance of CCode for computing the theta series of weight wt associated to the
    theta series and the matrix.
    num_of_procs: a positive integer. If it is greather than one, then it tries to
    compute by using multiple processes.
    separate_code: bool
    If true, separate code to sepate the main code to multiple files.
    '''
    tmp_var_name = "a"
    sum_tmp_var_name = "tmp"
    res_str_name = "res_str"
    if sty is None:
        sty = MpirStyle()
    vec_len = mat.ncols()

    if vec_len == 8:
        wt_small = tuple([a - 4 for a in wt])
        cache_vectors = "cache_vectors"
        num_of_vectors = "num_of_vectors"
    else:
        wt_small = tuple([a - 8 for a in wt])
        cache_vectors = "cache_vectors_rk16"
        num_of_vectors = "num_of_vectors_rk16"

    bdt_facs = _bideterminant_prime_factors_dict(mat, wt_small)
    _facs_pols = itertools.chain(*([a, b] for a, b in bdt_facs.values()))
    _facs_pols_lcm = lcm([b.denominator()
                          for b in itertools.chain(*(a.dict().values() for a in _facs_pols))])
    # Remove denominators
    bdt_facs = {k: (a * _facs_pols_lcm, b * _facs_pols_lcm) for k, (a, b) in bdt_facs.items()}

    Vrho = gl3_repn_module(wt_small)
    bs_pl_dct, bdt_var_dct = _pol_basis_as_polof_factors(wt_small, mat.base_ring())
    coef_pol_pairs = [bs_pl_dct[b] for b in Vrho.basis()]

    if real_part:
        coef_pols = [a for a, _ in coef_pol_pairs]
    else:
        coef_pols = [a for _, a in coef_pol_pairs]

    def _facs_pols_code_and_vars(i):
        return {k: pol_to_mpz_codes_and_result_var(v[i], tmp_var_name, bdt_var_dct[k][i],
                                                   algorithm='horner', sty=sty)
                for k, v in bdt_facs.items()}

    facs_pols_code_rl_dct = _facs_pols_code_and_vars(0)
    facs_pols_code_im_dct = _facs_pols_code_and_vars(1)

    res_vars = ["res" + str(i) for i, _ in enumerate(coef_pols)]

    coefs_pol_code_alst = [(pl, pol_to_mpz_codes_and_result_var(pl, tmp_var_name,
                                                                sum_tmp_var_name,
                                                                factor_pol=factor_pol,
                                                                algorithm=None,
                                                                sty=sty))
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

    coefs_pol_code_alst1 = [(pl, codes)
                            for (pl, (codes, _)), v in zip(coefs_pol_code_alst, res_vars)]

    _vrs = (tmp_vars + res_vars +
            sorted([str(r) for r, _ in bdt_var_dct.values()]) +
            sorted([str(im) for _, im in bdt_var_dct.values()]) +
            [str(a) for a in _s_t_u_ring(vec_len=vec_len).gens()])

    if separate_code:
        separated_codes = []
        alst = []
        for pl, codes in coefs_pol_code_alst1:
            s_cd, s_fname, s_args = _separated_code(func_name, _vrs, codes, pl)
            alst.append((pl, ["%s(%s)" % (s_fname, ", ".join(s_args))]))
            separated_codes.append(SepCCode(s_fname, s_cd, s_args))
        coefs_pol_code_alst1 = alst
    else:
        separated_codes = None

    # Add code for tmp *= num_of_classes_k[k] and res += + tmp
    for (_, codes), v in zip(coefs_pol_code_alst1, res_vars):
        if is_sparse_mat:
            codes.append(sty.mul_ui(sum_tmp_var_name, sum_tmp_var_name, "num_of_classes_i[i]"))
            codes.append(sty.mul_ui(sum_tmp_var_name, sum_tmp_var_name, "num_of_classes_j[j]"))
            codes.append(sty.mul_ui(sum_tmp_var_name, sum_tmp_var_name, "num_of_classes_k[k]"))
        codes.append(sty.add_z(v, v, sum_tmp_var_name))

    normalizing_num_var = "normalizig_num"
    _vrs.append(normalizing_num_var)

    init_code_str = _init_code(_vrs,
                               res_str_name,
                               vec_len,
                               _facs_pols_lcm ** wt_small[0],
                               normalizing_num_var,
                               sum(wt_small),
                               sty=sty, is_sparse_mat=is_sparse_mat)
    if is_sparse_mat:
        num_spaces = 14
    else:
        num_spaces = 18
    bi_det_factors_code_str = _bi_det_factors_code(facs_pols_code_rl_alst,
                                                   facs_pols_code_im_alst, _facs_pols_lcm,
                                                   bdt_var_dct, num_spaces)
    coeffs_code_str = _coeffs_code(coefs_pol_code_alst1, res_vars, num_spaces)

    str_concat_code_str = _str_concat_code([normalizing_num_var] + res_vars,
                                           res_str_name, sty=sty)

    cleanup_code_str = _cleanup_code(_vrs, sty=sty)

    vec_j_normalize_code = _vec_j_normalize_code(vec_len, is_sparse_mat)
    vec_i_normalize_code = _vec_i_normalize_code(vec_len, is_sparse_mat)
    k_inc_code = "k++" if num_of_procs == 1 else "k += %s" % (num_of_procs, )

    if is_sparse_mat:
        ith_vec_dict = {}
        for i in ["i", "j", "k"]:
            ith_vec_dict[i] = "reprs_{i}[{i}]".format(i=i)
        code = '''
char * {func_name}(int i_red, int a, int b, int c, int d, int e, int f)
{{
  /* mat: {mat_info}, quad_field: {quad_field_info}, real_part: {real_part} */
  /* young tableaux of the basis: {young_tableaux} */

  {cache_vectors}();

{init_code}

  for (int k = i_red; k < num_of_reprs_k; {k_inc_code})
    {{
{vec_j_normalize_code}
      for (int j = 0; j < num_of_reprs_j; j++)
        {{
{vec_i_normalize_code}
          for (int i = 0; i < num_of_reprs_i; i++)
            {{

{set_s_code}

{bi_det_factors_code}

{coeffs_code}
            }}
        }}
    }}

{str_concat_code}
{cleanup_code}

  /* The first element of {res_str_name} is a number for normalization and
    {res_str_name} must be freed later. */
  return {res_str_name};
}}
'''.format(init_code=init_code_str,
           bi_det_factors_code=bi_det_factors_code_str,
           coeffs_code=coeffs_code_str,
           str_concat_code=str_concat_code_str,
           cleanup_code=cleanup_code_str,
           res_str_name=res_str_name,
           func_name=func_name,

           vec_j_normalize_code=vec_j_normalize_code,
           vec_i_normalize_code=vec_i_normalize_code,

           mat_info=str(mat.list()),
           quad_field_info=str(mat.base_ring().polynomial()),
           real_part=str(real_part),
           young_tableaux=str([x.right_tableau.row_numbers for x in Vrho.basis()]),

           cache_vectors=cache_vectors,

           set_s_code=_set_s_code(vec_len, sty.set_si_func, ith_vec_dict, num_spaces),
           k_inc_code="k++" if num_of_procs == 1 else "k += %s" % (num_of_procs, ))
    else:
        # Non sparse matrix
        ith_vec_dict = {i: "vec_%s" % (i,) for i in ["i", "j", "k"]}
        code = '''
char * {func_name}(int i_red, int a, int b, int c, int d, int e, int f)
{{
  /* mat: {mat_info}, quad_field: {quad_field_info}, real_part: {real_part} */
  /* young tableaux of the basis: {young_tableaux} */

  {cache_vectors}();

{init_code}

  for (int k = i_red; k < {num_of_vectors}[c]; {k_inc_code}, vec_k += {vec_k_inc})
    {{
{vec_j_normalize_code}
      for (int j = 0; j < {num_of_vectors}[b]; j++, vec_j += {vec_len})
        {{
{vec_i_normalize_code}
          for (int i = 0; i < {num_of_vectors}[a]; i++, vec_i += {vec_len})
            {{
              if ((inner_prod({i_th_vec_having_norm_a}, {j_th_vec_having_norm_b}) == f) &&
                  (inner_prod({i_th_vec_having_norm_a}, {k_th_vec_having_norm_c}) == e) &&
                  (inner_prod({j_th_vec_having_norm_b}, {k_th_vec_having_norm_c}) == d))
                {{

{set_s_code}

{bi_det_factors_code}

{coeffs_code}
                }}
            }}
        }}
    }}

{str_concat_code}
{cleanup_code}

  /* The first element of {res_str_name} is a number for normalization and
    {res_str_name} must be freed later. */
  return {res_str_name};
}}
'''.format(init_code=init_code_str,
           bi_det_factors_code=bi_det_factors_code_str,
           coeffs_code=coeffs_code_str,
           str_concat_code=str_concat_code_str,
           cleanup_code=cleanup_code_str,
           res_str_name=res_str_name,
           func_name=func_name,

           vec_j_normalize_code=vec_j_normalize_code,
           vec_i_normalize_code=vec_i_normalize_code,

           mat_info=str(mat.list()),
           quad_field_info=str(mat.base_ring().polynomial()),
           real_part=str(real_part),
           young_tableaux=str([x.right_tableau.row_numbers for x in Vrho.basis()]),

           cache_vectors=cache_vectors,
           num_of_vectors=num_of_vectors,

           i_th_vec_having_norm_a=ith_vec_dict["i"],
           j_th_vec_having_norm_b=ith_vec_dict["j"],
           k_th_vec_having_norm_c=ith_vec_dict["k"],

           set_s_code=_set_s_code(vec_len, sty.set_si_func, ith_vec_dict, num_spaces),
           k_inc_code=k_inc_code,
           vec_k_inc=str(num_of_procs * vec_len), vec_len=vec_len)

    return CCode(code, separated_codes)


def _used_vars(vrs, codes):
    '''
    Return a sub list of vrs which are used in codes.
    '''
    code_str = " ".join(codes)
    res = []
    for a in vrs:
        p = re.compile(r"\b%s\b" % a)
        if p.search(code_str) is not None:
            res.append(a)
    return res


def _separated_code(func_name, vrs, codes, pl):
    '''
    Return code, func_name of separated code and args code.
    '''
    svrs = _used_vars(vrs, codes)
    sfunc_name = _separated_code_func_name(func_name, pl)
    args_code = ", ".join("mpz_t " + str(v) for v in svrs)
    body = "\n".join(c + ";" for c in codes)
    code = '''#include <mpir.h>

void {func_name}({args_code})
{{
{body}
}}
'''.format(func_name=sfunc_name, args_code=args_code, body=body)
    return (code, sfunc_name, svrs)


def _separated_code_func_name(func_name, pl):
    '''
    Return a unique name associated with func_name and codes.
    '''
    m = hashlib.sha256()
    m.update(func_name + str(pl))
    return func_name + m.hexdigest()
