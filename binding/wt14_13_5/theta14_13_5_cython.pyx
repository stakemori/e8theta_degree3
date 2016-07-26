from sage.rings.all import Integer, QQ, ZZ
from sage.misc.all import cached_function
from sage.modules.all import vector
from sage.matrix.all import MatrixSpace
from libc.stdlib cimport free
from e8theta_degree3.gl3_repn import GL3RepnElement
include "cysignals/signals.pxi"


cdef extern from "theta14_13_5_c.h":
    cpdef char * theta_c_14_13_5(int, int, int, int, int, int)


@cached_function
def theta(m):
    if not (m in MatrixSpace(QQ, 3) and (2 * m in MatrixSpace(ZZ, 3)) and
            (m[a, a] in ZZ for a in range(3)) and m.transpose() == m):
        raise ValueError("m must be a half integral matrix of size 3.")
    l = [m[i, i] for i in range(3)] + [2*m[t] for t in [(1, 2), (0, 2), (0, 1)]]
    a, b, c, d, e, f = [int(x) for x in l]
    if max([a, b, c]) > 7:
        raise ValueError("Diagonal elements are too large.")
    sig_on()
    cdef char* c_str = theta_c_14_13_5(a, b, c, d, e, f)
    cdef bytes py_str;
    try:
        py_str = c_str
    finally:
        free(c_str)
    sig_off()
    py_strs = py_str.split(",")
    res = [Integer(a) for a in py_strs]
    return vector(res)
