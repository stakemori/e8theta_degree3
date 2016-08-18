from sage.functions.all import ceil, floor, sgn
from sage.matrix.all import identity_matrix, matrix
from sage.rings.all import QQ


def _nearest_integer(x):
    r = floor(x)
    if x - r > 0.5:
        return r + 1
    else:
        return r


def _gaussian_reduction(b1, b2, S):
    '''
    b1, b2: vectors of length 3
    S: symmetric matrix of size 3
    '''
    while True:
        nb1 = b1 * S * b1
        nb2 = b2 * S * b2
        if nb2 < nb1:
            b1, b2 = b2, b1
        x = (b2 * S * b1) / (b1 * S * b1)
        r = _nearest_integer(x)
        a = b2 - r * b1
        if a * S * a >= b2 * S * b2:
            return (b1, b2)
        else:
            b1, b2 = a, b1


def _minkowski_reduction(b1, b2, b3, S):

    def inner_prod(x, y):
        return x * S * y

    while True:
        b1, b2, b3 = sorted([b1, b2, b3], key=lambda b: b * S * b)

        b1, b2 = _gaussian_reduction(b1, b2, S)

        b11 = inner_prod(b1, b1)
        b12 = inner_prod(b1, b2)
        b13 = inner_prod(b1, b3)
        b22 = inner_prod(b2, b2)
        b23 = inner_prod(b2, b3)
        b33 = inner_prod(b3, b3)

        y1 = - (b13 / b11 - b12 * b23 / (b11 * b22)) / \
            (1 - b12 ** 2 / (b11 * b22))
        y2 = - (b23 / b22 - b12 * b13 / (b11 * b22)) / \
            (1 - b12 ** 2 / (b11 * b22))

        # Find integers x1, x2 so that norm(b3 + x2 * b2 + x1 * b1) is minimal.
        a_norms_alst = []

        for x1 in [floor(y1), ceil(y1)]:
            for x2 in [floor(y2), ceil(y2)]:
                a = b3 + x2 * b2 + x1 * b1
                a_norms_alst.append((x1, x2, a, inner_prod(a, a)))
        _inner_prod_a = min(x[-1] for x in a_norms_alst)
        x1, x2, a, _ = next(x for x in a_norms_alst if x[-1] == _inner_prod_a)

        if _inner_prod_a >= b33:
            # Change sings of b1, b2, b3 and terminate the alogrithm
            sngs = [sgn(b12), sgn(b13), sgn(b23)]
            bs = [b1, b2, b3]
            try:
                # If b12, b13 or b23 is zero, change sgns of b1, b2, b3 so that
                # b12, b13, b23 >= 0.
                zero_i = sngs.index(0)
                set_ls = [set([1, 2]), set([1, 3]), set([2, 3])]
                t = set_ls[zero_i]
                _other = [x for x in [1, 2, 3] if x not in t][0]
                for x in t:
                    i = set_ls.index(set([x, _other]))
                    if sngs[i] < 0:
                        bs[x - 1] *= -1
                b1, b2, b3 = bs
            except ValueError:
                # Else change sgns so that b12, b13 > 0
                if b12 < 0:
                    b2 = -b2
                if b13 < 0:
                    b3 = -b3
            return (b1, b2, b3)
        else:
            b3 = a


def _minkowski_reduction_transform_matrix(S):
    '''
    Return a unimodular matrix u such that u^t * S * u is reduced in Minkowski's sense.
    '''
    b1, b2, b3 = identity_matrix(QQ, 3).columns()
    c1, c2, c3 = _minkowski_reduction(b1, b2, b3, S)
    return matrix([c1, c2, c3]).transpose()
