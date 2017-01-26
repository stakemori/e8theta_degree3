#define _GNU_SOURCE             /* for asprintf */
#include <stdio.h>
#include "e8vectors.h"
#include <stdlib.h>
#include <mpir.h>
#include "memory.h"
#include "vector_utils.h"

inline int inner_prod(int s[8], int t[8])
{
  return (s[0]*t[0] + s[1]*t[1] + s[2]*t[2] + s[3]*t[3] +
          s[4]*t[4] + s[5]*t[5] + s[6]*t[6] + s[7]*t[7]) >> 2;
}




char * theta_c_12_12_12(int i_red, int a, int b, int c, int d, int e, int f)
{
  /* mat: [1, 0, 0, i, 0, 0, 0, 0, 0, 1, 0, 0, i, 0, 0, 0, 0, 0, 1, 0, 0, i, 0, 0], quad_field: x^2 + 1, real_part: True */
  /* young tableaux of the basis: [[[1, 1, 1, 1, 1, 1, 1, 1], [2, 2, 2, 2, 2, 2, 2, 2], [3, 3, 3, 3, 3, 3, 3, 3]]] */

  cache_vectors();

  mpz_t a0;
  mpz_t a1;
  mpz_t a2;
  mpz_t a3;
  mpz_t a4;
  mpz_t a5;
  mpz_t a6;
  mpz_t a7;
  mpz_t a8;
  mpz_t tmp;
  mpz_t res0;
  mpz_t real_part0;
  mpz_t imag_part0;
  mpz_t s0;
  mpz_t s1;
  mpz_t s2;
  mpz_t s3;
  mpz_t s4;
  mpz_t s5;
  mpz_t s6;
  mpz_t s7;
  mpz_t t0;
  mpz_t t1;
  mpz_t t2;
  mpz_t t3;
  mpz_t t4;
  mpz_t t5;
  mpz_t t6;
  mpz_t t7;
  mpz_t u0;
  mpz_t u1;
  mpz_t u2;
  mpz_t u3;
  mpz_t u4;
  mpz_t u5;
  mpz_t u6;
  mpz_t u7;
  mpz_t normalizig_num;

  mpz_init(a0);
  mpz_init(a1);
  mpz_init(a2);
  mpz_init(a3);
  mpz_init(a4);
  mpz_init(a5);
  mpz_init(a6);
  mpz_init(a7);
  mpz_init(a8);
  mpz_init(tmp);
  mpz_init(res0);
  mpz_init(real_part0);
  mpz_init(imag_part0);
  mpz_init(s0);
  mpz_init(s1);
  mpz_init(s2);
  mpz_init(s3);
  mpz_init(s4);
  mpz_init(s5);
  mpz_init(s6);
  mpz_init(s7);
  mpz_init(t0);
  mpz_init(t1);
  mpz_init(t2);
  mpz_init(t3);
  mpz_init(t4);
  mpz_init(t5);
  mpz_init(t6);
  mpz_init(t7);
  mpz_init(u0);
  mpz_init(u1);
  mpz_init(u2);
  mpz_init(u3);
  mpz_init(u4);
  mpz_init(u5);
  mpz_init(u6);
  mpz_init(u7);
  mpz_init(normalizig_num);

  char * res_str;
  mpz_set_str(normalizig_num, "1", 10);
  mpz_mul_2exp(normalizig_num, normalizig_num, 24);


  static int reprs_k[MAX_NM_REPRS][8];
  static unsigned int num_of_classes_k[MAX_NM_REPRS];

  static int reprs_j[MAX_NM_REPRS][8];
  static unsigned int num_of_classes_j[MAX_NM_REPRS];
  static int vecs_j[MAX_NM_OF_VECTORS][8];
  int num_of_vecs_j;

  static int reprs_i[MAX_NM_REPRS][8];
  static unsigned int num_of_classes_i[MAX_NM_REPRS];
  static int vecs_i[MAX_NM_OF_VECTORS][8];

  int num_of_reprs_k = repr_modulo_autom(c, reprs_k, num_of_classes_k);
  int num_of_reprs_j;


  for (int k = i_red; k < num_of_reprs_k; k += 8)
    {

      int wo_sign_indices_array[8][16] = {{0}};
      int w_sign_indices[8] = {0};
      set_w_sign_indices(w_sign_indices, reprs_k[k], 8, 2);
      set_wo_sign_indices_array(wo_sign_indices_array, reprs_k[k], 8, 2);

      num_of_vecs_j = 0;
      int * cached_vec_b = cached_vectors_ptr[b];
      for (int j = 0; j < num_of_vectors[b]; j++, cached_vec_b += 8)
        {
          if (inner_prod(cached_vec_b, reprs_k[k]) == d)
            {
              memcpy(vecs_j[num_of_vecs_j++], cached_vec_b, sizeof(int) * 8);
            }
        }
      num_of_reprs_j = repr_modulo_autom_w_indices(vecs_j, num_of_vecs_j, reprs_j,
                                                   num_of_classes_j,
                                                   w_sign_indices, wo_sign_indices_array);

      for (int j = 0; j < num_of_reprs_j; j++)
        {

          int wo_sign_indices_array[8][16] = {{0}};
          int w_sign_indices[8] = {0};
          set_wo_sign_indices_array2(wo_sign_indices_array, reprs_j[j], reprs_k[k], 8, 2);
          set_w_sign_indices_2(w_sign_indices, reprs_j[j], reprs_k[k], 8, 2);

          int num_of_vecs_i = 0;
          int * cached_vec_a = cached_vectors_ptr[a];
          for (int l = 0; l < num_of_vectors[a]; l++, cached_vec_a += 8)
            {
              if ((inner_prod(cached_vec_a, reprs_k[k]) == e) &&
                  (inner_prod(cached_vec_a, reprs_j[j]) == f))
                {
                  memcpy(vecs_i[num_of_vecs_i++], cached_vec_a, sizeof(int) * 8);
                }
            }
          int num_of_reprs_i = repr_modulo_autom_w_indices(vecs_i, num_of_vecs_i, reprs_i,
                                                           num_of_classes_i,
                                                           w_sign_indices,
                                                           wo_sign_indices_array);

          for (int i = 0; i < num_of_reprs_i; i++)
            {

              mpz_set_si(s0, reprs_i[i][0]);
              mpz_set_si(s1, reprs_i[i][1]);
              mpz_set_si(s2, reprs_i[i][2]);
              mpz_set_si(s3, reprs_i[i][3]);
              mpz_set_si(s4, reprs_i[i][4]);
              mpz_set_si(s5, reprs_i[i][5]);
              mpz_set_si(s6, reprs_i[i][6]);
              mpz_set_si(s7, reprs_i[i][7]);

              mpz_set_si(t0, reprs_j[j][0]);
              mpz_set_si(t1, reprs_j[j][1]);
              mpz_set_si(t2, reprs_j[j][2]);
              mpz_set_si(t3, reprs_j[j][3]);
              mpz_set_si(t4, reprs_j[j][4]);
              mpz_set_si(t5, reprs_j[j][5]);
              mpz_set_si(t6, reprs_j[j][6]);
              mpz_set_si(t7, reprs_j[j][7]);

              mpz_set_si(u0, reprs_k[k][0]);
              mpz_set_si(u1, reprs_k[k][1]);
              mpz_set_si(u2, reprs_k[k][2]);
              mpz_set_si(u3, reprs_k[k][3]);
              mpz_set_si(u4, reprs_k[k][4]);
              mpz_set_si(u5, reprs_k[k][5]);
              mpz_set_si(u6, reprs_k[k][6]);
              mpz_set_si(u7, reprs_k[k][7]);

              /* Computation of real_part0 = real part of 1 * (x02*x11*x20 - x01*x12*x20 - x02*x10*x21 + x00*x12*x21 + x01*x10*x22 - x00*x11*x22) */
              mpz_set(a0, s0);
              mpz_neg(a1, s1);
              mpz_set(a2, s3);
              mpz_neg(a3, s4);
              mpz_mul(a3, a3, t0);
              mpz_addmul(a3, a2, t1);
              mpz_addmul(a3, a1, t3);
              mpz_addmul(a3, a0, t4);
              mpz_neg(a0, s0);
              mpz_set(a1, s2);
              mpz_neg(a2, s3);
              mpz_set(a4, s5);
              mpz_mul(a4, a4, t0);
              mpz_addmul(a4, a2, t2);
              mpz_addmul(a4, a1, t3);
              mpz_addmul(a4, a0, t5);
              mpz_set(a0, s1);
              mpz_neg(a1, s2);
              mpz_set(a2, s4);
              mpz_neg(a5, s5);
              mpz_mul(a5, a5, t1);
              mpz_addmul(a5, a2, t2);
              mpz_addmul(a5, a1, t4);
              mpz_addmul(a5, a0, t5);
              mpz_set(a0, s3);
              mpz_neg(a1, s4);
              mpz_neg(a2, s0);
              mpz_set(a6, s1);
              mpz_mul(a6, a6, t0);
              mpz_addmul(a6, a2, t1);
              mpz_addmul(a6, a1, t3);
              mpz_addmul(a6, a0, t4);
              mpz_neg(a0, s3);
              mpz_set(a1, s5);
              mpz_set(a2, s0);
              mpz_neg(a7, s2);
              mpz_mul(a7, a7, t0);
              mpz_addmul(a7, a2, t2);
              mpz_addmul(a7, a1, t3);
              mpz_addmul(a7, a0, t5);
              mpz_set(a0, s4);
              mpz_neg(a1, s5);
              mpz_neg(a2, s1);
              mpz_set(a8, s2);
              mpz_mul(a8, a8, t1);
              mpz_addmul(a8, a2, t2);
              mpz_addmul(a8, a1, t4);
              mpz_addmul(a8, a0, t5);
              mpz_mul(a8, a8, u0);
              mpz_addmul(a8, a7, u1);
              mpz_addmul(a8, a6, u2);
              mpz_addmul(a8, a5, u3);
              mpz_addmul(a8, a4, u4);
              mpz_addmul(a8, a3, u5);
              mpz_set(real_part0, a8);


              /* Computation of imag_part0 = imaginary part of 1 * (x02*x11*x20 - x01*x12*x20 - x02*x10*x21 + x00*x12*x21 + x01*x10*x22 - x00*x11*x22) */
              mpz_set(a0, s3);
              mpz_neg(a1, s4);
              mpz_neg(a2, s0);
              mpz_set(a3, s1);
              mpz_mul(a3, a3, t0);
              mpz_addmul(a3, a2, t1);
              mpz_addmul(a3, a1, t3);
              mpz_addmul(a3, a0, t4);
              mpz_neg(a0, s3);
              mpz_set(a1, s5);
              mpz_set(a2, s0);
              mpz_neg(a4, s2);
              mpz_mul(a4, a4, t0);
              mpz_addmul(a4, a2, t2);
              mpz_addmul(a4, a1, t3);
              mpz_addmul(a4, a0, t5);
              mpz_set(a0, s4);
              mpz_neg(a1, s5);
              mpz_neg(a2, s1);
              mpz_set(a5, s2);
              mpz_mul(a5, a5, t1);
              mpz_addmul(a5, a2, t2);
              mpz_addmul(a5, a1, t4);
              mpz_addmul(a5, a0, t5);
              mpz_neg(a0, s0);
              mpz_set(a1, s1);
              mpz_neg(a2, s3);
              mpz_set(a6, s4);
              mpz_mul(a6, a6, t0);
              mpz_addmul(a6, a2, t1);
              mpz_addmul(a6, a1, t3);
              mpz_addmul(a6, a0, t4);
              mpz_set(a0, s0);
              mpz_neg(a1, s2);
              mpz_set(a2, s3);
              mpz_neg(a7, s5);
              mpz_mul(a7, a7, t0);
              mpz_addmul(a7, a2, t2);
              mpz_addmul(a7, a1, t3);
              mpz_addmul(a7, a0, t5);
              mpz_neg(a0, s1);
              mpz_set(a1, s2);
              mpz_neg(a2, s4);
              mpz_set(a8, s5);
              mpz_mul(a8, a8, t1);
              mpz_addmul(a8, a2, t2);
              mpz_addmul(a8, a1, t4);
              mpz_addmul(a8, a0, t5);
              mpz_mul(a8, a8, u0);
              mpz_addmul(a8, a7, u1);
              mpz_addmul(a8, a6, u2);
              mpz_addmul(a8, a5, u3);
              mpz_addmul(a8, a4, u4);
              mpz_addmul(a8, a3, u5);
              mpz_set(imag_part0, a8);

              /* Computation of res0 = real_part0^8 - 28*real_part0^6*imag_part0^2 + 70*real_part0^4*imag_part0^4 - 28*real_part0^2*imag_part0^6 + imag_part0^8 */
              mpz_set(a0, imag_part0);
              mpz_mul_ui(a1, real_part0, 28); mpz_neg(a1, a1);
              mpz_mul(a1, a1, real_part0);
              mpz_addmul(a1, a0, imag_part0);
              mpz_mul(a1, a1, imag_part0);
              mpz_mul_ui(a0, real_part0, 70);
              mpz_mul(a0, a0, real_part0);
              mpz_mul(a0, a0, real_part0);
              mpz_mul(a0, a0, real_part0);
              mpz_addmul(a0, a1, imag_part0);
              mpz_mul(a0, a0, imag_part0);
              mpz_mul_ui(a1, real_part0, 28); mpz_neg(a1, a1);
              mpz_mul(a1, a1, real_part0);
              mpz_mul(a1, a1, real_part0);
              mpz_mul(a1, a1, real_part0);
              mpz_mul(a1, a1, real_part0);
              mpz_mul(a1, a1, real_part0);
              mpz_addmul(a1, a0, imag_part0);
              mpz_mul(a1, a1, imag_part0);
              mpz_set(a0, real_part0);
              mpz_mul(a0, a0, real_part0);
              mpz_mul(a0, a0, real_part0);
              mpz_mul(a0, a0, real_part0);
              mpz_mul(a0, a0, real_part0);
              mpz_mul(a0, a0, real_part0);
              mpz_mul(a0, a0, real_part0);
              mpz_mul(a0, a0, real_part0);
              mpz_addmul(a0, a1, imag_part0);
              mpz_set(tmp, a0);
              mpz_mul_ui(tmp, tmp, num_of_classes_i[i]);
              mpz_mul_ui(tmp, tmp, num_of_classes_j[j]);
              mpz_mul_ui(tmp, tmp, num_of_classes_k[k]);
              mpz_add(res0, res0, tmp);
            }
        }
    }

  int buf_size = asprintf(&res_str, "%s,%s", mpz_get_str(NULL, 10, normalizig_num), mpz_get_str(NULL, 10, res0));

  if (buf_size == -1) {
    exit(1);
  }

  mpz_clear(a0);
  mpz_clear(a1);
  mpz_clear(a2);
  mpz_clear(a3);
  mpz_clear(a4);
  mpz_clear(a5);
  mpz_clear(a6);
  mpz_clear(a7);
  mpz_clear(a8);
  mpz_clear(tmp);
  mpz_clear(res0);
  mpz_clear(real_part0);
  mpz_clear(imag_part0);
  mpz_clear(s0);
  mpz_clear(s1);
  mpz_clear(s2);
  mpz_clear(s3);
  mpz_clear(s4);
  mpz_clear(s5);
  mpz_clear(s6);
  mpz_clear(s7);
  mpz_clear(t0);
  mpz_clear(t1);
  mpz_clear(t2);
  mpz_clear(t3);
  mpz_clear(t4);
  mpz_clear(t5);
  mpz_clear(t6);
  mpz_clear(t7);
  mpz_clear(u0);
  mpz_clear(u1);
  mpz_clear(u2);
  mpz_clear(u3);
  mpz_clear(u4);
  mpz_clear(u5);
  mpz_clear(u6);
  mpz_clear(u7);
  mpz_clear(normalizig_num);

  /* The first element of res_str is a number for normalization and
    res_str must be freed later. */
  return res_str;
}
