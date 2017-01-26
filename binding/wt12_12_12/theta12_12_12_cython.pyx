import itertools
from multiprocessing import Pool
from sage.rings.all import Integer, QQ, ZZ
from sage.misc.all import cached_function
from sage.modules.all import vector
from sage.matrix.all import MatrixSpace
from libc.stdlib cimport free
from e8theta_degree3.gl3_repn import GL3RepnElement
include "cysignals/signals.pxi"

cdef extern from "theta12_12_12_c.h":
    cpdef char * theta_c_12_12_12(int, int, int, int, int, int, int)


def _theta12_12_12_cython_part(i_red_m):
    i_red, m = i_red_m

    if not (m in MatrixSpace(QQ, 3) and (2 * m in MatrixSpace(ZZ, 3)) and
            (m[a, a] in ZZ for a in range(3)) and m.transpose() == m):
        raise ValueError("m must be a half integral matrix of size 3.")
    l = [m[i, i] for i in range(3)] + [2*m[t] for t in [(1, 2), (0, 2), (0, 1)]]
    a, b, c, d, e, f = [int(x) for x in l]
    if max([a, b, c]) > 7:
        raise ValueError("Diagonal elements are too large.")
    sig_on()
    cdef char* c_str = theta_c_12_12_12(i_red, a, b, c, d, e, f)
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
def theta12_12_12_cython(m):
    p = Pool(processes=8)
    try:
        res = sum(p.map(_theta12_12_12_cython_part, zip(range(8), itertools.repeat(m, 8))))
    except KeyboardInterrupt:
        p.terminate()
        p.join()
    finally:
        p.close()
        p.join()
    return res

