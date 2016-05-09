from sage.modules.all import vector


# Same as degree2.utils.find_linearly_indep_indices
def find_linearly_indep_indices(vectors, r):
    '''
    Let vectors be a list of vectors or a list of list.
    Assume r be the rank of vectors.
    This function returns a list of indices I of length r
    such that the rank of [vectors[i] for i in I] is equal to r.
    '''
    acc = []
    ncls = len(vectors[0])
    if isinstance(vectors[0], list):
        vectors = [vector(a) for a in vectors]
    while True:
        if r == 0:
            return acc
        nrws = len(vectors)
        for a, i in zip(vectors, range(nrws)):
            if a != 0:
                first = a
                first_r_idx = i
                break

        for j in range(ncls):
            if not first[j] == 0:
                a = first[j]
                v = a ** (-1) * first
                nonzero_col_index = j
                break

        vectors1 = []
        for j in range(first_r_idx + 1, nrws):
            w = vectors[j]
            vectors1.append(w - w[nonzero_col_index] * v)
        vectors = vectors1
        r -= 1
        if acc == []:
            acc.append(first_r_idx)
        else:
            acc.append(first_r_idx + acc[-1] + 1)
