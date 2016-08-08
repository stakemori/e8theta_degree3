#include <math.h>
#include "rank16_vectors.h"
#include <stdio.h>
#include <stdlib.h>
#include "vector_utils.h"
#include "memory.h"


/* Fourier coefficients of Eisenstein series of weight 8 */
int num_of_vectors_rk16[8] = {1, 480, 61920, 1050240, 7926240, 37500480, 135480960, 395301120};

Rk16VecInt cached_vectors_rk16[MAX_NORM_RK16 + 1][MAX_NM_OF_VECTORS_RK16][16];
static int cached_idx[MAX_NORM_RK16 + 1] = {0};

static void _cache_vectors(void)
{
  double m = sqrt(2 * MAX_NORM_RK16);
  int _s0_end = ceil(2 * m + 1);
  int _tmp;
  for (int s0 = floor(-2*m); s0 < _s0_end; s0++)
    {
      int nrm0 = s0 * s0;
      double ds0 = s0;
      int beg = floor(-m - ds0/2);
      int end = ceil(m - ds0/2 + 1);
      for (int s1 = beg; s1 < end; s1++)
        {
          _tmp = s0 + 2 * s1;
          int nrm1 = nrm0 + _tmp * _tmp;
          if (nrm1 < 8 * MAX_NORM_RK16 + 1)
            {
              for (int s2 = beg; s2 < end; s2++)
                {
                  _tmp = s0 + 2 * s2;
                  int nrm2 = nrm1 + _tmp * _tmp;
                  if (nrm2 < 8 * MAX_NORM_RK16 + 1)
                    {
                      for (int s3 = beg; s3 < end; s3++)
                        {
                          _tmp = s0 + 2 * s3;
                          int nrm3 = nrm2 + _tmp * _tmp;
                          if (nrm3 < 8 * MAX_NORM_RK16 + 1)
                            {
                              for (int s4 = beg; s4 < end; s4++)
                                {
                                  _tmp = s0 + 2 * s4;
                                  int nrm4 = nrm3 + _tmp * _tmp;
                                  if (nrm4 < 8 * MAX_NORM_RK16 + 1)
                                    {
                                      for (int s5 = beg; s5 < end; s5++)
                                        {
                                          _tmp = s0 + 2 * s5;
                                          int nrm5 = nrm4 + _tmp * _tmp;
                                          if (nrm5 < 8 * MAX_NORM_RK16 + 1)
                                            {
                                              for (int s6 = beg; s6 < end; s6++)
                                                {
                                                  _tmp = s0 + 2 * s6;
                                                  int nrm6 = nrm5 + _tmp * _tmp;
                                                  if (nrm6 < 8 * MAX_NORM_RK16 + 1)
                                                    {
                                                      for (int s7 = beg; s7 < end; s7++)
                                                        {
                                                          _tmp = s0 + 2 * s7;
                                                          int nrm7 = nrm6 + _tmp * _tmp;
                                                          if (nrm7 < 8 * MAX_NORM_RK16 + 1)
                                                            {
                                                              for (int s8 = beg; s8 < end; s8++)
                                                                {
                                                                  _tmp = s0 + 2 * s8;
                                                                  int nrm8 = nrm7 + _tmp * _tmp;
                                                                  if (nrm8 < 8 * MAX_NORM_RK16 + 1)
                                                                    {
                                                                      for (int s9 = beg; s9 < end; s9++)
                                                                        {
                                                                          _tmp = s0 + 2 * s9;
                                                                          int nrm9 = nrm8 + _tmp * _tmp;
                                                                          if (nrm9 < 8 * MAX_NORM_RK16 + 1)
                                                                            {
                                                                              for (int s10 = beg; s10 < end; s10++)
                                                                                {
                                                                                  _tmp = s0 + 2 * s10;
                                                                                  int nrm10 = nrm9 + _tmp * _tmp;
                                                                                  if (nrm10 < 8 * MAX_NORM_RK16 + 1)
                                                                                    {
                                                                                      for (int s11 = beg; s11 < end; s11++)
                                                                                        {
                                                                                          _tmp = s0 + 2 * s11;
                                                                                          int nrm11 = nrm10 + _tmp * _tmp;
                                                                                          if (nrm11 < 8 * MAX_NORM_RK16 + 1)
                                                                                            {
                                                                                              for (int s12 = beg; s12 < end; s12++)
                                                                                                {
                                                                                                  _tmp = s0 + 2 * s12;
                                                                                                  int nrm12 = nrm11 + _tmp * _tmp;
                                                                                                  if (nrm12 < 8 * MAX_NORM_RK16 + 1)
                                                                                                    {
                                                                                                      for (int s13 = beg; s13 < end; s13++)
                                                                                                        {
                                                                                                          _tmp = s0 + 2 * s13;
                                                                                                          int nrm13 = nrm12 + _tmp * _tmp;
                                                                                                          if (nrm13 < 8 * MAX_NORM_RK16 + 1)
                                                                                                            {
                                                                                                              for (int s14 = beg; s14 < end; s14++)
                                                                                                                {
                                                                                                                  _tmp = s0 + 2 * s14;
                                                                                                                  int nrm14 = nrm13 + _tmp * _tmp;
                                                                                                                  if (nrm14 < 8 * MAX_NORM_RK16 + 1)
                                                                                                                    {
                                                                                                                      double _centr = (2 * (s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14) + ds0)/4;
                                                                                                                      int _end = ceil(m/2 - _centr + 1);
                                                                                                                      for (int s15 = floor(-m/2 - _centr); s15 < _end; s15++)
                                                                                                                        {
                                                                                                                          _tmp = (s0 + 2 * (s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8 + s9 + s10 + s11 + s12 + s13 + s14) + 4 * s15);
                                                                                                                          int nrm15 = nrm14 + _tmp * _tmp;
                                                                                                                          if (nrm15 < 8 * MAX_NORM_RK16 + 1)
                                                                                                                            {
                                                                                                                              int v[16] = {s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15};
                                                                                                                              int _nm = nrm15 >> 3;
                                                                                                                              int idx = cached_idx[_nm]++;
                                                                                                                              for (int i = 0; i < 16; i++)
                                                                                                                                {
                                                                                                                                  cached_vectors_rk16[_nm][idx][i] = v[i];
                                                                                                                                }
                                                                                                                            }
                                                                                                                        }
                                                                                                                    }
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void cache_vectors_rk16(void)
{
  int b = 1;
  for (int i = 0; i < 16; i++)
    {
      b = (! cached_vectors_rk16[1][0][i]) & b;
    }
  if (b)
    {
      _cache_vectors();
    }
}


/* convert vec to an element of a sub lattice of (2Z)^16 */
/* (s0, ... , s15) to (s0, s1 + 2s0, ..., s0 + 2s14, s0 + 2(s1 + ... + s14) + 4s15) */
void _convert_to_euclid_vector_rk16(Rk16VecInt vec[16]){
  vec[15] = 4 * vec[15];
  for (int i = 1; i < 15; i++)
    {
      vec[15] += 2 * vec[i];
    }
  for (int i = 1; i < 15; i++)
    {
      vec[i] = 2 * vec[i] + vec[0];
    }
  vec[15] += vec[0];
}

/* inverse to _covert_to_euclid_vector */
void _convert_from_euclid_vector_rk16(Rk16VecInt vec[16]){
  for (int i = 1; i < 15; i++)
    {
      vec[i] = (vec[i] - vec[0]) >> 1;
    }
  vec[15] -= vec[0];
  for (int i = 1; i < 15; i++)
    {
      vec[15] -= 2*vec[i];
    }
  vec[15] = vec[15] >> 2;
}

void normalize_vec_rk16(Rk16VecInt vec[16])
{
  int mul = 1;
  int vec1[9];
  for (int i = 0; i < 9; i++)
    {
      vec1[i] = vec[i + 7];
    }
  for (int i = 0; i < 9; i++)
    {
      if (vec1[i] == 0)
        {
          mul = 0;
          break;
        }
      else if (vec1[i] < 0)
        {
          mul *= -1;
        }
    }
  for (int i = 0; i < 9; i++)
    {
      vec1[i] = abs(vec1[i]);
    }
  sort_int_vec(vec1, 9);
  if (mul == -1)
    {
      vec1[8] *= -1;
    }
  for (int i = 0; i < 9; i++)
    {
      vec[i + 7] = vec1[i];
    }
}

int repr_modulo_autom_rk16(int n, int repr[MAX_NM_REPRS_RK16][17])
/* Let L be the even unimodular lattice of rank 16 in the euclid space
   G a subgroup of autom group generated by permutations of entries and even sign changes that
   do not change the first 7 entries.
   repr: an array
   Let _repr the set of representatives cached_vectors_rk16[n]/G.
   Set repr {(s0, ..., s15, a) for (s0, ..., s15) in _repr}, where a is the number
   of vectors in L equivalent to the representative (s0, ..., s15).
   Finally, return the number of representatives.
 */
{
  int num = 0;
  Rk16VecInt vec[16] = {0};
  for (int i = 0; i < num_of_vectors_rk16[n]; i++)
    {
      memcpy(vec, cached_vectors_rk16[n][i], sizeof(Rk16VecInt) * 16);
      _convert_to_euclid_vector_rk16(vec);
      normalize_vec_rk16(vec);
      _convert_from_euclid_vector_rk16(vec);
      int found = 0;
      for (int j = 0; j < num; j++)
        {
          int bl = 1;
          for (int k = 0; k < 16; k++)
            {
              bl = bl & (repr[j][k] == vec[k]);
            }
          if (bl)
            {
              found = 1;
              repr[j][16]++;
              break;
            }
        }
      if (! found)
        {
          for (int j = 0; j < 16; j++)
            {
              repr[num][j] = vec[j];
            }
          repr[num][16] = 1;
          num++;
        }
    }
  return num;
}
