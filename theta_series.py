try:
    # pylint: disable=no-name-in-module
    from e8theta_degree3.binding.e8theta import miyawaki_theta
except:
    raise
from e8theta_degree3.gl3_repn import element_constructor


def miyawaki_theta_dict(Ts, verbose=False):
    '''
    Ts: a list of instances of HalfIntMatElement.
    Return a dict of Fourier coefficients.
    Its key is an instance of HalfIntMatElement and its value is
    an instance of ReplSpaceElement.
    '''
    Vrho_const = element_constructor((12, 12, 12))
    res = {}
    for T in Ts:
        if verbose:
            print T
        res[T] = Vrho_const([miyawaki_theta(T.T)])
    return res
