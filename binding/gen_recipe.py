'''
Recipe for generating c and cython sources.
'''
from __future__ import print_function

import subprocess
from os.path import dirname, join

import e8theta_degree3
from e8theta_degree3.binding.theta_code_format import generate_cython_and_build_scripts
from sage.rings.all import QuadraticField, ZZ
from sage.matrix.all import matrix, diagonal_matrix
from sage.misc.all import random
from sage.functions.all import floor


def _gen_base(wt, mats, cyfns, cfns, real_parts=None, is_sparse_mat=False, num_of_procs=8,
              separate_code=False):
    assert all(mat * mat.transpose() == 0 for mat in mats)
    cyf, cf = _names(wt)
    generate_cython_and_build_scripts(_recipe_dir(wt),
                                      cyf, cf, cyfns, cfns,
                                      wt, mats, overwrite=True, num_of_procs=num_of_procs,
                                      real_parts=real_parts,
                                      is_sparse_mat=is_sparse_mat, separate_code=separate_code)


def _recipe_dir(wt):
    wt_s = "_".join(str(a) for a in wt)
    binding_dir = join(dirname(e8theta_degree3.__file__), "binding")
    return join(binding_dir, "wt{wt}/".format(wt=wt_s))


def _names(wt):
    '''
    Return cython file name, c file name.
    '''
    wt_s = "_".join(str(a) for a in wt)
    return ("theta{wt}_cython".format(wt=wt_s),
            "theta{wt}".format(wt=wt_s))


def _cython_func_name_default(wt):
    wt_s = "_".join(str(a) for a in wt)
    return "theta%s_cython" % (wt_s, )


def _c_func_name_default(wt):
    wt_s = "_".join(str(a) for a in wt)
    return "theta_c_%s" % (wt_s, )


def gen_wt12():
    '''
    Recipe for weight (12, 12, 12). This is for test.
    '''
    wt = (12, 12, 12)
    i = QuadraticField(-1, name="i").gen()
    mat = matrix(3, [1, 0, 0, i, 0, 0, 0, 0, 0, 1, 0, 0, i, 0, 0, 0, 0, 0, 1, 0, 0, i, 0, 0])
    _gen_base(wt, [mat], [_cython_func_name_default(wt)], [_c_func_name_default(wt)],
              is_sparse_mat=True)


def gen_wt14():
    wt = (14, 14, 14)
    i = QuadraticField(-1, name="i").gen()
    mat0 = matrix(3, [-5, -i, 5, 9 * i, -4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      3, 5 * i, 0, -2 * i, -2, 0, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      4, 6 * i, 0, -2 * i, -3, -i, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    _gen_base(wt, [mat0],
              [_cython_func_name_default(wt)], [_c_func_name_default(wt)],
              is_sparse_mat=True, num_of_procs=8)


def gen_wt14_13_5():
    '''
    Recipe for weight (14, 13, 5).
    '''
    wt = (14, 13, 5)
    i = QuadraticField(-1, name="i").gen()
    mat = matrix(3, [2, 3 * i, -2, -i, -1, 0, -1, 0, -2, -i, 0,
                     -i, -1, 0, -1, -2 * i, -1, 0, -1, 0, 0, i, 0, -i])
    _gen_base(wt, [mat], [_cython_func_name_default(wt)], [_c_func_name_default(wt)])


def gen_wt16_16_14():
    '''
    Recipe for weight (16, 16, 14).
    '''
    wt = (16, 16, 14)
    i = QuadraticField(-1, name="i").gen()
    mat0 = matrix(3, [-5, -i, 5, 9 * i, -4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      3, 5 * i, 0, -2 * i, -2, 0, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      4, 6 * i, 0, -2 * i, -3, -i, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    mat1 = matrix(3, [3, 7 * i, -3, i, -4, 0, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      3, 5 * i, 0, -2 * i, -2, 0, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      -5, -3 * i, -6, 8 * i, 0, 2 * i, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    _gen_base(wt, [mat0, mat1],
              [_cython_func_name_default(wt) + "_" + str(a) for a in range(2)],
              [_c_func_name_default(wt) + "_" + str(a) for a in range(2)],
              is_sparse_mat=True, num_of_procs=8)


def gen_wt16_13_7():
    '''
    Recipe for weight (16, 13, 7).
    '''
    wt = (16, 13, 7)
    i = QuadraticField(-1, name="i").gen()
    mat0 = matrix(3, [-3, -i, 2, 4*i, -2, 0, 0, 0, 2, 4*i, -2, -4*i, 5, i,
                      0, 0, 5, -3*i, 6, -2*i, 1, 7*i, 0, 0])
    mat1 = matrix(3, [1, 3*i, -2, 0, -2, 0, 0, 0, -7, -5*i, 4, 2*i, 0,
                      -6*i, 0, 0, 5, 3*i, -1, -3*i, -1, 3*i, 0, 0])
    _gen_base(wt, [mat0, mat1],
              [_cython_func_name_default(wt) + "_" + str(a) for a in range(2)],
              [_c_func_name_default(wt) + "_" + str(a) for a in range(2)],
              is_sparse_mat=True, num_of_procs=8, separate_code=True)


def gen_wt18_13_5():
    '''
    Recipe for weight (18, 13, 5).
    '''
    wt = (18, 13, 5)
    i = QuadraticField(-1, name="i").gen()
    mat = matrix(3, [1, 0, 0, i, 0, 0, 0, 0, 0, 1, 0, 0, i, 0, 0, 0, 0, 0, 1, 0, 0, i, 0, 0])
    _gen_base(wt, [mat], [_cython_func_name_default(wt)], [_c_func_name_default(wt)],
              is_sparse_mat=True, separate_code=True)


T0 = matrix([[ZZ(1), ZZ(1) / ZZ(2), ZZ(1) / ZZ(2)],
             [ZZ(1) / ZZ(2), ZZ(1), ZZ(1) / ZZ(2)],
             [ZZ(1) / ZZ(2), ZZ(1) / ZZ(2), ZZ(1)]])

T1 = diagonal_matrix([ZZ(1), ZZ(1), ZZ(1)])

Ts = [T0]


def _rank(funcs):
    m = matrix([sum([f(t).list() for t in Ts], []) for f in funcs])
    return m.rank()


def compute_rank_in_subprocess(dir_name, mod_name, func_names):
    funcs = ["e8theta_degree3.binding.{dir_name}.{mod_name}.{func}".format(
        dir_name=dir_name, mod_name=mod_name, func=func) for func in func_names]
    funcs_str = "[" + ", ".join(funcs) + "]"
    cmd = ("""sage -c 'import sys; sys.path.append("{pkg_dir}");""" +
           "import e8theta_degree3.binding.{dir_name}.{mod_name};" +
           "from e8theta_degree3.binding.gen_recipe import _rank;" +
           "print _rank({funcs})'").format(
               pkg_dir=dirname(dirname(e8theta_degree3.__file__)),
               dir_name=dir_name, mod_name=mod_name,
               funcs=funcs_str)
    return subprocess.check_output(cmd, shell=True).strip()


def _find_mat_wt(wt, mats_base, mats_total, dim,
                 cython_func_names=None, c_func_names=None,
                 is_sparse_mat=False, separate_code=False):
    ln = len(mats_total)
    if cython_func_names is None:
        cython_func_names = [_cython_func_name_default(wt) + "_" + str(i)
                             for i in range(dim)]
    if c_func_names is None:
        c_func_names = [_c_func_name_default(wt) + "_" + str(i)
                        for i in range(dim)]
    while True:
        mat = mats_total[floor(random() * ln)]
        print(mat.list())
        wt_str = "_".join([str(a) for a in wt])
        _gen_base(wt, mats_base + [mat], cython_func_names, c_func_names,
                  is_sparse_mat=is_sparse_mat, separate_code=separate_code)
        po = subprocess.Popen("make compile-cython", shell=True, cwd=_recipe_dir(wt))
        po.wait()
        print("Building done.")
        cyf, _, = _names(wt)
        a = compute_rank_in_subprocess("wt%s" % (wt_str, ), cyf, cython_func_names)
        print(a)
        a = ZZ(a)
        if a == dim:
            print(mat.list())
            print("found")
            return
