try:
    # pylint: disable=no-name-in-module
    from e8theta_degree3.binding.e8theta import miyawaki_theta
except:
    raise
from e8theta_degree3.gl3_repn import GL3RepnElement


def miyawaki_theta_dict(Ts, verbose=False):
    '''
    Ts: a list of instances of HalfIntMatElement.
    Return a dict of Fourier coefficients.
    Its key is an instance of HalfIntMatElement and its value is
    an instance of ReplSpaceElement.
    '''
    res = {}
    for T in Ts:
        if verbose:
            print T
        res[T] = GL3RepnElement([miyawaki_theta(T.T)], (12, 12, 12))
    return res
