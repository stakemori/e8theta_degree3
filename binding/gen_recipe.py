'''
Recipe for generating c and cython sources.
'''
from __future__ import print_function

import itertools
import subprocess
from os.path import dirname, join

import e8theta_degree3
from e8theta_degree3.binding.theta_code_format import generate_cython_and_build_scripts
from e8theta_degree3.young_tableau import poly_repn_dim
from sage.rings.all import QuadraticField, ZZ
from sage.matrix.all import matrix, diagonal_matrix
from sage.misc.all import random
from sage.functions.all import floor


def _gen_base(wt, mat, suffix="", real_part=True):
    assert mat * mat.transpose() == 0
    cyf, cf, cyfn, cfn = _names(wt, suffix)
    generate_cython_and_build_scripts(_recipe_dir(wt),
                                      cyf, cf, cyfn, cfn,
                                      wt, mat, overwrite=True, num_of_procs=8,
                                      real_part=real_part)


def _recipe_dir(wt):
    wt_s = "_".join(str(a) for a in wt)
    binding_dir = join(dirname(e8theta_degree3.__file__), "binding")
    return join(binding_dir, "wt{wt}/".format(wt=wt_s))


def _names(wt, suffix):
    '''
    Return cython file name, c file name, cython function name, c function name.
    '''
    wt_s = "_".join(str(a) for a in wt)
    return ("theta{wt}{sfx}_cython".format(wt=wt_s, sfx=suffix),
            "theta{wt}{sfx}".format(wt=wt_s, sfx=suffix),
            "theta{sfx}".format(sfx=suffix),
            "theta_c_{wt}{sfx}".format(wt=wt_s, sfx=suffix))


def gen_wt12():
    '''
    Recipe for weight (12, 12, 12). This is for test.
    '''
    i = QuadraticField(-1, name="i").gen()
    mat = matrix(3, [1, 0, 0, i, 0, 0, 0, 0, 0, 1, 0, 0, i, 0, 0, 0, 0, 0, 1, 0, 0, i, 0, 0])
    _gen_base((12, 12, 12), mat)


def gen_wt14_13_5():
    '''
    Recipe for weight (14, 13, 5).
    '''
    i = QuadraticField(-1, name="i").gen()
    mat = matrix(3, [2, 3 * i, -2, -i, -1, 0, -1, 0, -2, -i, 0,
                     -i, -1, 0, -1, -2 * i, -1, 0, -1, 0, 0, i, 0, -i])
    _gen_base((14, 13, 5), mat)


T0 = matrix([[ZZ(1), ZZ(1) / ZZ(2), ZZ(1) / ZZ(2)],
             [ZZ(1) / ZZ(2), ZZ(1), ZZ(1) / ZZ(2)],
             [ZZ(1) / ZZ(2), ZZ(1) / ZZ(2), ZZ(1)]])

T1 = diagonal_matrix([ZZ(1), ZZ(1), ZZ(1)])


def compute_initial_fc_in_subprocess(dir_name, mod_name, func_name):
    return subprocess.check_output(
        ("""sage -c 'import sys; sys.path.append("{pkg_dir}");""" +
         "import e8theta_degree3.binding.{dir_name}.{mod_name};" +
         "from e8theta_degree3.binding.gen_recipe import T0;" +
         "print e8theta_degree3.binding.{dir_name}.{mod_name}.{func_name}(T0)'").format(
             pkg_dir=dirname(dirname(e8theta_degree3.__file__)),
             dir_name=dir_name, mod_name=mod_name, func_name=func_name), shell=True).strip()


def _find_mat_wt(wt, mats, suffix=""):
    ln = len(mats)
    while True:
        mat = mats[floor(random() * ln)]
        print(mat.list())
        wt_str = "_".join([str(a) for a in wt])
        for real_part in [True, False]:
            _gen_base(wt, mat, suffix=suffix, real_part=real_part)
            po = subprocess.Popen("make compile-cython", shell=True, cwd=_recipe_dir(wt))
            po.wait()
            print("Building done.")
            cyf, _, cyfn, _ = _names(wt, suffix)
            a = compute_initial_fc_in_subprocess(
                "wt%s" % (wt_str, ), cyf, cyfn)
            print(a)
            zero_vec_str = "(" + ", ".join(itertools.repeat("0", poly_repn_dim(wt))) + ")"
            if a != zero_vec_str:
                print(mat.list())
                print("found")
                return
