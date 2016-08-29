#include <math.h>
#include "rank16_vectors.h"
#include <stdio.h>
#include <stdlib.h>
#include "vector_utils.h"
#include "memory.h"


/* Fourier coefficients of Eisenstein series of weight 8 */
int num_of_vectors_rk16[8] = {1, 480, 61920, 1050240, 7926240, 37500480, 135480960, 395301120};

static int cached_vectors0[1][16] = {0};
static int cached_vectors1[480][16] = {0};
static int cached_vectors2[61920][16] = {0};
static int cached_vectors3[1050240][16] = {0};
static int cached_vectors4[7926240][16] = {0};

int *cached_vectors_rk16_ptr[] = {cached_vectors0[0], cached_vectors1[0],
                                  cached_vectors2[0], cached_vectors3[0],
                                  cached_vectors4[0]};

static int cached_idx[MAX_NORM_RK16 + 1] = {0};

/* convert vec to an element of a sub lattice of (2Z)^16 */
/* (s0, ... , s15) to (s0, s1 + 2s0, ..., s0 + 2s14, s0 + 2(s1 + ... + s14) + 4s15) */
static void _convert_to_euclid_vector_rk16(int vec[16]){
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
                                                                                                                              _convert_to_euclid_vector_rk16(v);
                                                                                                                              int _nm = nrm15 >> 3;
                                                                                                                              int idx = cached_idx[_nm]++;
                                                                                                                              memcpy(cached_vectors_rk16_ptr[_nm] + 16 * idx, v, 16 * sizeof(int));
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
      b = (! (cached_vectors_rk16_ptr[1] + 16 * 0)[i]) && b;
    }
  if (b)
    {
      _cache_vectors();
    }
}


inline int vec_equal(const int vec1[16], const int vec2[16])
{
  return !(memcmp(vec1, vec2, 16 * sizeof(int)));
}

int repr_modulo_autom_rk16(int n, int reprs[MAX_NM_REPRS_RK16][16],
                           unsigned int num_of_classes[MAX_NM_REPRS_RK16])
/* Let L be the even unimodular lattice of rank 16 in the euclid space
   G a subgroup of autom group generated by permutations of entries, even sign changes that
   do not change the first 7 entries and -1.
   reprs: an array
   Set reprs the set of representatives cached_vectors_rk16[n]/G.
   Set num_of_repres so that ith element of num_of_repres is the number of
   vectors in L L equivalent to the ith element of reprs.
   Finally, return the number of representatives.
*/
{
  int num = 0;
  int found;
  int vec[16] = {0};
  int * ptr = cached_vectors_rk16_ptr[n];
  for (int i = 0; i < num_of_vectors_rk16[n]; i++, ptr += 16)
    {
      memcpy(vec, ptr, sizeof(int) * 16);
      normalize_vec_last_len_elts(vec, 16, 9);
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
          memcpy(reprs[num], vec, 16 * sizeof(int));
          num_of_classes[num] = 1;
          num++;
        }
    }
  return num;
}


inline int repr_modulo_autom_rk16_w_indices(int vecs[MAX_NM_OF_VECTORS_RK16][16],
                                            int num_of_vecs,
                                            int reprs[MAX_NM_REPRS_RK16][16],
                                            unsigned int num_of_classes[MAX_NM_REPRS_RK16],
                                            int w_sign_indices[16],
                                            int wo_sign_indices_array[8][16])
/* Similar to repr_modulo_autom_rk16 but the representative is computed by
   normalize_vec_w_indices and cached_vectors_rk16[n] is replaced by vecs.
   num_of_vecs is the actual length of vecs.*/
{
  int num = 0;
  int found;
  int vec[16] = {0};
  for (int i = 0; i < num_of_vecs; i++)
    {
      memcpy(vec, vecs[i], sizeof(int) * 16);
      normalize_vec_w_indices(vec, w_sign_indices, wo_sign_indices_array);
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
          memcpy(reprs[num], vec, 16 * sizeof(int));
          num_of_classes[num] = 1;
          num++;
        }
    }
  return num;
}

void set_w_sign_indices_rk16(int indices[16], const int vec[16])
{
  int idx = 0;
  /* 7 is hard-coded */
  for (int k = 7; k < 16; k++)
    {
      if (vec[k] == 0)
        {
          indices[idx++] = k;
        }
    }
}

void set_w_sign_indices_rk16_2(int indices[16], const int vec1[16],
                               const int vec2[16])
{
  int idx = 0;
  /* 7 is hard-coded */
  for (int k = 7; k < 16; k++)
    {
      if ((vec1[k] == 0) && (vec2[k] == 0))
        {
          indices[idx++] = k;
        }
    }
}

void set_wo_sign_indices_array(int indices_array[8][16], const int vec[16])
{
  int max = vec[0];
  int min = vec[0];
  /* Set max and min of vec. */
  for (int i = 0; i < 16; i++)
    {
      if (vec[i] > max)
        {
          max = vec[i];
        }
      if (vec[i] < min)
        {
          min = vec[i];
        }
    }

  int idx = 0;
  for (int k = min; k < max + 1; k++)
    {
      if (k)
        {
          int count = 0;
          int idx_vec[16] = {0};
          for (int l = 7; l < 16; l++)
            {
              if (k == vec[l])
                {
                  idx_vec[count++] = l;
                }
            }
          if (count > 1)
            {
              memcpy(indices_array[idx++], idx_vec, sizeof(int) * count);
            }
        }
    }
}

void set_wo_sign_indices_array2(int indices_array[8][16], const int vec1[16],
                                const int vec2[16])
{
  int indices_array2[8][16] = {0};
  set_wo_sign_indices_array(indices_array2, vec2);
  int max = vec1[0];
  int min = vec1[0];
  /* Set max and min of vec1. */
  for (int i = 7; i < 16; i++)
    {
      if (vec1[i] > max)
        {
          max = vec1[i];
        }
      if (vec1[i] < min)
        {
          min = vec1[i];
        }
    }

  int idx = 0;

  /* non-zero entries of vec2 and any entries of vec1. */
  for (int idx2 = 0; indices_array2[idx2][0]; idx2++)
    {
      for (int i = min; i < max + 1; i++)
        {
          int count = 0;
          int idx_vec[16] = {0};
          for (int j = 0; indices_array2[idx2][j]; j++)
            {
              if (i == vec1[indices_array2[idx2][j]])
                {
                  idx_vec[count++] = indices_array2[idx2][j];
                }
            }
          if (count > 1)
            {
              memcpy(indices_array[idx++], idx_vec, sizeof(int) * count);
            }
        }
    }

  /* zero entries of vec2 and any enties of vec1. */
  int zero_idcs[16] = {0};
  set_w_sign_indices_rk16(zero_idcs, vec2);
  for (int i = min; i < max + 1; i++)
    {
      int count = 0;
      int idx_vec[16] = {0};
      for (int j = 0; zero_idcs[j]; j++)
        {
          if (i == vec1[zero_idcs[j]])
            {
              idx_vec[count++] = zero_idcs[j];
            }
        }
      if (count > 1)
        {
          memcpy(indices_array[idx++], idx_vec, sizeof(int) * count);
        }
    }
}
/* Local Variables: */
/* compile-command: "cd ..; make compile-theta_vectors" */
/* End: */
