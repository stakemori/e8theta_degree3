from os.path import dirname, join
import e8theta_degree3
from sage.misc.all import cached_function


@cached_function
def data_dir():
    return join(dirname(e8theta_degree3.__file__), "results", "data")
