from sage.rings.all import Integer
from sage.misc.all import cached_function
from libc.stdlib cimport free
include "cysignals/signals.pxi"


cdef extern from "miyawaki_theta.h":
    cpdef char * miyawaki_theta_c(int, int, int, int, int, int)


@cached_function
def miyawaki_theta(m):
    try:
        l = [m[i, i] for i in range(3)] + [2*m[t] for t in [(1, 2), (0, 2), (0, 1)]]
        a, b, c, d, e, f = [int(x) for x in l]
    except:
        raise ValueError
    if max([a, b, c]) > 7:
        raise ValueError
    sig_on()
    cdef char* c_str = miyawaki_theta_c(a, b, c, d, e, f)
    cdef bytes py_str;
    try:
        py_str = c_str
    finally:
        free(c_str)
    res = Integer(py_str)
    sig_off()
    return res

# Local Variables:
# compile-command: "cd ..; make compile-cython"
# End:
