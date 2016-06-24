from sage.rings.all import Integer
from sage.misc.all import cached_function
from libc.stdlib cimport free
include "cysignals/signals.pxi"


cdef extern from "miyawaki_theta.h":
    cpdef char * miyawaki_theta_c(int, int, int, int, int, int)


@cached_function
def miyawaki_theta(a, b, c, d, e, f):
    try:
        l = [a, b, c, d, e, f]
        a, b, c, d, e, f = [int(x) for x in l]
    except:
        raise ValueError
    if max([a, b, c]) > 7:
        raise ValueError
    sig_on()
    s = miyawaki_theta_c(a, b, c, d, e, f)
    res = Integer(s)
    free(s)
    sig_off()
    return res

# Local Variables:
# compile-command: "cd ..; make compile-cython"
# End:
