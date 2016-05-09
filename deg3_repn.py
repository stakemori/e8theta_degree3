from sage.rings.all import PolynomialRing, QQ
from sage.matrix.all import matrix
from sage.misc.all import cached_function
from e8theta_degree3.young_tableau import YoungTableu


@cached_function
def matrix_var(base_field=QQ):
    R = PolynomialRing(base_field, names=[
                       'x%s%s' % (i, j) for i in range(3) for j in range(3)])
    return matrix(3, R.gens())


def _bideterminant(a, b):
    '''
    a, b: an instance of young_tableau.YoungTableu
    '''
    m = matrix_var()
    res = 1
    for l1, l2 in zip(a.col_numbers, b.col_numbers):
        res *= m.matrix_from_rows_and_columns(
            [i - 1 for i in l1], [j - 1 for j in l2]).det()
    return res


def _t_lambda(wt):
    return YoungTableu(n=3, row_numbers=[[i + 1 for _ in range(a)]
                                         for i, a in enumerate(wt)])
