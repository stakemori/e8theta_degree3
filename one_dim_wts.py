# -*- coding: utf-8 -*-

'''
cf. Bergström, Faber, van der Geer, Siegel modualr forms of degree three and
the cohomology of local systems table 1
'''


def conjectual_one_dim_wts():
    cojectual_one_dim_wts = [(12, 4, 4), (11, 6, 3), (10, 8, 2), (14, 4, 4), (13, 6, 3),
                             (12, 8, 2), (11, 7, 4), (10, 8, 4), (13, 6, 1),
                             (16, 5, 1), (15, 6, 1), (14, 7, 1), (13, 9, 0), (13, 5, 4),
                             (12, 7, 3), (11, 10, 1), (14, 4, 2), (12, 6, 2), (10, 6, 4),
                             (8, 6, 6), (17, 4, 1), (16, 4, 2), (15, 5, 2), (14, 6, 2),
                             (13, 8, 1), (12, 6, 4), (11, 9, 2), (13, 4, 3), (10, 9, 1),
                             (15, 4, 3), (14, 5, 3), (13, 7, 2), (11, 8, 3)]
    return [(a + 4, b + 4, c + 4) for a, b, c in cojectual_one_dim_wts]

# [(t, poly_repn_dim(t)) for t in conjectual_one_dim_wts()]
# =>
# [((12, 10, 10), 6),
#  ((16, 8, 8), 45),
#  ((14, 12, 8), 60),
#  ((14, 10, 8), 60),
#  ((18, 8, 8), 66),
#  ((15, 11, 8), 90),
#  ((17, 9, 8), 99),
#  ((14, 13, 5), 99),
#  ((14, 12, 6), 105),
#  ((16, 10, 8), 105),
#  ((15, 10, 7), 120),
#  ((15, 14, 5), 120),
#  ((17, 8, 7), 120),
#  ((15, 12, 7), 120),
#  ((15, 13, 6), 132),
#  ((16, 11, 7), 165),
#  ((19, 8, 7), 168),
#  ((17, 10, 7), 192),
#  ((18, 9, 7), 195),
#  ((16, 12, 6), 210),
#  ((16, 10, 6), 210),
#  ((18, 8, 6), 231),
#  ((17, 11, 6), 273),
#  ((20, 8, 6), 312),
#  ((18, 10, 6), 315),
#  ((19, 9, 6), 330),
#  ((17, 10, 5), 336),
#  ((17, 12, 5), 336),
#  ((17, 13, 4), 375),
#  ((18, 11, 5), 420),
#  ((19, 10, 5), 480),
#  ((21, 8, 5), 504),
#  ((20, 9, 5), 510)]

# cojecture 7.7 (a)
# [(17, 13, 4), (15, 13, 6), (14, 13, 5)]
