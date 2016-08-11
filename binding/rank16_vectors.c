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

void _normalize_vec_w_sign(Rk16VecInt vec[16], size_t len)
/* Let G be an subgroup of the auto group of the lattice generated by
 even sign changes and permutations which change only first len element.
This function converts vec to the representative by the action of G.*/
{
  int mul = 1;
  int len_int = (int)len;
  for (int i = 0; i < len_int; i++)
    {
      if (vec[i] == 0)
        {
          mul = 0;
          break;
        }
      else if (vec[i] < 0)
        {
          mul *= -1;
        }
    }
  for (int i = 0; i < len_int; i++)
    {
      vec[i] = abs(vec[i]);
    }
  sort_int_vec(vec, len);
  if (mul == -1)
    {
      vec[0] *= -1;
    }
}

void _normalize_vec_wo_sign(Rk16VecInt vec[16], size_t len)
/* Similar to _normalize_vec_w_sign but only for permutations. */
{
  sort_int_vec(vec, len);
}

void normalize_vec_rk16_w_indices(Rk16VecInt vec[16],
                                  int w_sign_indices[16],
                                  int wo_sign_indices_array[8][16])
/* Let v be an int array of length 16 s.t.
   v = {i_1, ..., i_a, 0, ..., 0} where 0 < i_k < 16.
   Denote by G(v) a subgroup of the autom group of the lattice
   generated by even sign changes and permutations which change only
   {i_1, ... i_a}th entries.
   Denote by H(v) a similar group for only permutations.

   Let G be a group generated by G(w_sign_indices) and
   H(v) for non-zero array v in wo_sign_indices_array.
   This function converts vec to the representative by the action of G.*/
{
  Rk16VecInt vec1[16];
  int len = 0;
  for (int i = 0; w_sign_indices[i]; i++)
    {
      vec1[i] = vec[w_sign_indices[i]];
      len++;
    }
  _normalize_vec_w_sign(vec1, len);
  for (int i = 0; w_sign_indices[i]; i++)
    {
      vec[w_sign_indices[i]] = vec1[i];
    }

  for (int j = 0; wo_sign_indices_array[j][0]; j++)
    {
      len = 0;
      for (int i = 0; wo_sign_indices_array[j][i]; i++)
        {
          vec1[i] = vec[wo_sign_indices_array[j][i]];
          len++;
        }
      _normalize_vec_wo_sign(vec1, len);
      for (int i = 0; wo_sign_indices_array[j][i]; i++)
        {
          vec[wo_sign_indices_array[j][i]] = vec1[i];
        }
    }
}

void normalize_vec_rk16_last9(Rk16VecInt vec[16])
{
  /* Action of -1. */
  int first_nonzero_idx;
  for (first_nonzero_idx = 0; first_nonzero_idx < 16; first_nonzero_idx++)
    {
      if (vec[first_nonzero_idx] != 0)
        {
          break;
        }
    }

  if (vec[first_nonzero_idx] < 0)
    {
      for (int i = 0; i < 16; i++)
        {
          vec[i] = - vec[i];
        }
    }

  /* Action of the last 9 elements */
  Rk16VecInt vec1[9];
  for (int i = 0; i < 9; i++)
    {
      vec1[i] = vec[i + 7];
    }
  _normalize_vec_w_sign(vec1, 9);
  for (int i = 0; i < 9; i++)
    {
      vec[i + 7] = vec1[i];
    }
}

static inline int vec_equal(Rk16VecInt const vec1[16], Rk16VecInt const vec2[16])
{
  return !(memcmp(vec1, vec2, 16 * sizeof(Rk16VecInt)));
}

int repr_modulo_autom_rk16(int n, int reprs[MAX_NM_REPRS_RK16][16], int num_of_classes[MAX_NM_REPRS_RK16])
/* Let L be the even unimodular lattice of rank 16 in the euclid space
   G a subgroup of autom group generated by permutations of entries and even sign changes that
   do not change the first 7 entries.
   reprs: an array
   Set reprs the set of representatives cached_vectors_rk16[n]/G.
   Set num_of_repres so that ith element of num_of_repres is the number of
   vectors in L L equivalent to the ith element of reprs.
   Finally, return the number of representatives.
 */
{
  int num = 0;
  int found;
  Rk16VecInt vec[16] = {0};
  for (int i = 0; i < num_of_vectors_rk16[n]; i++)
    {
      memcpy(vec, cached_vectors_rk16[n][i], sizeof(Rk16VecInt) * 16);
      _convert_to_euclid_vector_rk16(vec);
      normalize_vec_rk16_last9(vec);
      _convert_from_euclid_vector_rk16(vec);
      found = 0;
      for (int j = 0; j < num; j++)
        {
          found = vec_equal(reprs[j], vec);
          if (found)
            {
              num_of_classes[j]++;
              break;
            }
        }
      if (! found)
        {
          memcpy(reprs[num], vec, 16 * sizeof(Rk16VecInt));
          num_of_classes[num] = 1;
          num++;
        }
    }
  return num;
}

int repr_modulo_autom_rk16_w_indices(Rk16VecInt vecs[MAX_NM_OF_VECTORS_RK16][16],
                                     int num_of_vecs,
                                     int reprs[MAX_NM_REPRS_RK16][16],
                                     int num_of_classes[MAX_NM_REPRS_RK16],
                                     int w_sign_indices[16],
                                     int wo_sign_indices_array[8][16])
/* Similar to repr_modulo_autom_rk16 but the representative is computed by
normalize_vec_rk16_w_indices and cached_vectors_rk16[n] is replaced by vecs.
num_of_vecs is the actual length of vecs.*/
{
  int num = 0;
  int found;
  Rk16VecInt vec[16] = {0};
  for (int i = 0; i < num_of_vecs; i++)
    {
      memcpy(vec, vecs[i], sizeof(Rk16VecInt) * 16);
      _convert_to_euclid_vector_rk16(vec);
      normalize_vec_rk16_w_indices(vec, w_sign_indices, wo_sign_indices_array);
      _convert_from_euclid_vector_rk16(vec);
      found = 0;
      for (int j = 0; j < num; j++)
        {
          found = vec_equal(reprs[j], vec);
          if (found)
            {
              num_of_classes[j]++;
              break;
            }
        }
      if (! found)
        {
          memcpy(reprs[num], vec, 16 * sizeof(Rk16VecInt));
          num_of_classes[num] = 1;
          num++;
        }
    }
  return num;
}
