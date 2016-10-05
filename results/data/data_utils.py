from os.path import dirname, join
import e8theta_degree3
from sage.misc.all import cached_function


@cached_function
def data_dir():
    return join(dirname(e8theta_degree3.__file__), "results", "data")


def half_int_mat_to_list(t):
    return ([t.T[a, a] for a in range(3)] +
            [2 * t.T[i, j] for i, j in [(1, 2), (0, 2), (0, 1)]])


def sort_ts(ts):
    def key(t):
        T = t.T
        return (T.det(), T[0, 0], T[1, 1], T[2, 2], T[0, 1], T[0, 2], T[1, 2])
    return list(sorted(ts, key=key))
