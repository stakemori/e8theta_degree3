#include "e8vectors.h"
#include "rank16_vectors.h"
#include <stdio.h>
#include "memory.h"
#include <math.h>
#include "vector_utils.h"
#include "mpir.h"

Rk16VecInt inner_prod_rk16(Rk16VecInt s[16], Rk16VecInt t[16])
{
  int sum = 0;
  for (int i = 0; i < 16; i++)
    {
      sum += s[i] * t[i];
    }
  return sum >> 2;
}

static void print_vec(int * vec, int a)
{
  for (int i = 0; i < a; i++)
    {
      printf("%d, ", vec[i]);
    }
  printf("\n");
}

int test_norm_vec_rk16(void)
{
  cache_vectors_rk16();
  Rk16VecInt vec1[16];
  Rk16VecInt vec2[16];
  Rk16VecInt vec3[16];
  Rk16VecInt vec4[16];
  int bl = 1;
  for (int i = 0; i < num_of_vectors_rk16[2]; i++)
    {
      memcpy(vec1, cached_vectors_rk16[2][i], sizeof(int) * 16);
      memcpy(vec2, vec1, sizeof(int) * 16);
      normalize_vec_rk16_last9(vec1);
      memcpy(vec3, vec1, sizeof(int) * 16);
      memcpy(vec4, vec2, sizeof(int) * 16);

      int mul1 = 1;
      int mul2 = 1;
      for (int j = 0; j < 16; j++)
        {
          mul1 *= vec1[j];
        }
      for (int j = 0; j < 16; j++)
        {
          mul2 *= vec2[j];
        }

      /* Test first if 7 entries are only sign change. */
      for (int j = 0; j < 7; j++)
        {
          bl = bl && (vec1[j] == vec2[j]);
        }
      if (! bl)
        {
          bl = 1;
          for (int j = 0; j < 7; j++)
            {
              bl = bl && (vec1[j] == -vec2[j]);
            }
          if (! bl) {
              printf("False: unchagne.\n");
              print_vec(vec1, 16);
              print_vec(vec2, 16);
              return 0;
          }
        }

      for (int j = 0; j < 16; j++)
        {
          vec1[j] = abs(vec1[j]);
          vec2[j] = abs(vec2[j]);
        }

      sort_int_vec(vec1, 16);
      sort_int_vec(vec2, 16);

      for (int j = 0; j < 16; j++)
        {
          if (vec1[j] != vec2[j])
            {
              printf("False: sorting\n");
              print_vec(vec1, 16);
              print_vec(vec2, 16);
              return 0;
            }
        }
      if (mul1 != mul2)
        {
          print_vec(vec3, 16);
          print_vec(vec4, 16);
          printf("False: sign\n");
          return 0;
        }
    }
  return 1;
}

int test_repr_rk16(void)
{
  static Rk16VecInt _reprs_rk16[MAX_NM_REPRS_RK16][16];
  static unsigned int num_of_reprs_rk16[MAX_NM_REPRS_RK16];
  int n = 2;
  int num;
  Rk16VecInt vec[16];
  for (int i = 0; i < n + 1; i++)
    {
      num = repr_modulo_autom_rk16(i, _reprs_rk16, num_of_reprs_rk16);
      printf("%d\n", num);
    }
  int s = 0;
  for (int i = 0; i < num; i++)
    {
      s += num_of_reprs_rk16[i];
      memcpy(vec, _reprs_rk16[i], sizeof(Rk16VecInt) * 16);
      if (inner_prod_rk16(vec, vec) != 2 * n)
        {
          printf("False: norm\n");
          return 0;
        }
      print_vec(vec, 16);
    }
  if (s == num_of_vectors_rk16[n])
    {
      return 1;
    }
  else
    {
      return 0;
    }
  /* cache_vectors_rk16(); */
  /* Rk16VecInt vec[16]; */
  /* for (int i = 0; i < num_of_vectors_rk16[2]; i++) */
  /*   { */
  /*     for (int j = 0; j < 16; j++) */
  /*       { */
  /*         vec[j] = cached_vectors_rk16[2][i][j]; */
  /*       } */
  /*     print_vec(vec, 16); */
  /*     normalize_vec_rk16_last9(vec); */
  /*     print_vec(vec, 16); */
  /*   } */

}


void number_of_loops(int a, int b, int c, int d, int e, int f)
{
  cache_vectors_rk16();
  static Rk16VecInt _reprs[MAX_NM_REPRS_RK16][16];
  static unsigned int num_of_classes[MAX_NM_REPRS_RK16];
  int num_of_reprs = repr_modulo_autom_rk16(c, _reprs, num_of_classes);

  mpz_t res, one;
  mpz_init(res);
  mpz_init(one);
  mpz_set_si(one, 1);
  for (int i = 0; i < num_of_vectors_rk16[a]; i++)
    {
      for (int j = 0; j < num_of_vectors_rk16[b]; j++)
        {
          for (int k = 0; k < num_of_reprs; k++)
            {
              if (inner_prod_rk16(cached_vectors_rk16[a][i], cached_vectors_rk16[b][j]) == f)
                {
                  if (inner_prod_rk16(cached_vectors_rk16[a][i], _reprs[k]) == e)
                    {
                      if (inner_prod_rk16(cached_vectors_rk16[b][j], _reprs[k]) == d)
                        {
                          mpz_add(res, res, one);
                        }
                    }
                }
            }
        }
    }
  printf("%s\n", mpz_get_str(NULL, 10, res));
  mpz_clear(res);
}

int vec_equal(Rk16VecInt * vec1, Rk16VecInt * vec2, int start, int end)
{
  int is_equal = 1;
  for (int i = start; i < end; i++)
    {
      is_equal = is_equal && (vec1[i] == vec2[i]);
    }
  return is_equal;
}

void abs_vec(Rk16VecInt * vec, int len)
{
  for (int i = 0; i < len; i++)
    {
      vec[i] = abs(vec[i]);
    }
}

int equal_as_multiset(Rk16VecInt vec1[16], Rk16VecInt vec2[16], int start, int end)
{
  Rk16VecInt _vec1[16];
  Rk16VecInt _vec2[16];
  for (int i = start; i < end; i++)
    {
      _vec1[i - start] = vec1[i];
      _vec2[i - start] = vec2[i];
    }
  sort_int_vec(_vec1, end - start);
  sort_int_vec(_vec2, end - start);
  return vec_equal(_vec1, _vec2, 0, end - start);
}

int test_normalize_vec_rk16_w_indices(void)
{
  cache_vectors_rk16();
  int zero_idcs[16] = {0};
  for (int i = 7; i < 12; i++)
    {
      zero_idcs[i-7] = i;
    }
  int non_zero_idcss[8][16] = {0};
  for (int i = 12; i < 16; i++)
    {
      non_zero_idcss[0][i - 12] = i;
    }
  Rk16VecInt vec[16];
  Rk16VecInt vec1[16];
  int num = 2;
  for (int i = 0; i < num_of_vectors_rk16[num]; i++)
    {
      memcpy(vec, cached_vectors_rk16[num][i], 16 * sizeof(Rk16VecInt));
      memcpy(vec1, vec, 16 * sizeof(Rk16VecInt));
      normalize_vec_rk16_w_indices(vec, zero_idcs, non_zero_idcss);
      /* printf("raw:\n"); */
      /* print_vec(vec1, 16); */
      /* print_vec(vec, 16); */
      if (! vec_equal(vec, vec1, 0, 7))
        {
          printf("False: unchanged.\n");
          print_vec(vec1, 16);
          print_vec(vec, 16);
          return 0;
        }
      if (! equal_as_multiset(vec, vec1, 12, 16))
        {
          printf("False: multiset.\n");
          print_vec(vec1, 16);
          print_vec(vec, 16);
          return 0;
        }
      /* Destructive operator */
      abs_vec(vec, 16);
      abs_vec(vec1, 16);
      if (! equal_as_multiset(vec, vec1, 7, 12))
        {
          printf("False: multiset sign.\n");
          print_vec(vec1, 16);
          print_vec(vec, 16);
          return 0;
        }
    }
  return 1;
}

void test_normalize_vec_rk16_w_indices_1(void)
{
  int vec[16] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2};
  int vec1[16] = {0, 0, 0, 0, 0, 0, 0, 3, 4, 5, 1, 3, -100, 100, 30, -9};
  int wo_sign_indices_array[8][16] = {0};
  int idx_array[16] = {0};
  set_w_sign_indices_rk16(idx_array, vec);
  set_wo_sign_indices_array(wo_sign_indices_array, vec);
  normalize_vec_rk16_w_indices(vec1, idx_array, wo_sign_indices_array);
  print_vec(vec1, 16);
  for (int j = 0; wo_sign_indices_array[j][0]; j++)
    {
      print_vec(wo_sign_indices_array[j], 16);
    }
}

void test_repr_rk16_w_idices(void)
{
  cache_vectors_rk16();
  static Rk16VecInt reprs1[MAX_NM_REPRS_RK16][16];
  static unsigned int num_of_classes1[MAX_NM_REPRS_RK16];
  static Rk16VecInt _reprs2[1050240][16];
  static unsigned int num_of_classes2[1050240];
  static Rk16VecInt vecs[MAX_NM_OF_VECTORS_RK16][16];
  int num_of_vecs;
  int num = 2;
  int num_of_reprs1 = repr_modulo_autom_rk16(num, reprs1, num_of_classes1);
  mpz_t res;
  mpz_t one;
  mpz_init(res);
  mpz_init(one);
  mpz_set_si(one, 1);
  Rk16VecInt vec[16];
  printf("%d\n", num_of_reprs1);
  int num_of_reprs2;
  for (int i = 0; i < num_of_reprs1; i++)
    {
      int wo_sign_indices_array[8][16] = {0};
      int w_sign_indices[16] = {0};
      memcpy(vec, reprs1[i], 16 * sizeof(Rk16VecInt));
      printf("i: %d\n", i);
      print_vec(vec, 16);
      set_w_sign_indices_rk16(w_sign_indices, vec);
      set_wo_sign_indices_array(wo_sign_indices_array, vec);
      print_vec(w_sign_indices, 16);
      for (int j = 0; wo_sign_indices_array[j][0]; j++)
        {
          print_vec(wo_sign_indices_array[j], 16);
        }
      num_of_vecs = 0;
      for (int k = 0; k < num_of_vectors_rk16[num]; k++)
        {
          if (inner_prod_rk16(cached_vectors_rk16[num][k], reprs1[i]) == 2)
            {
              memcpy(vecs[num_of_vecs++], cached_vectors_rk16[num][k], sizeof(Rk16VecInt) * 16);
            }
        }
      num_of_reprs2 = repr_modulo_autom_rk16_w_indices(vecs, num_of_vecs, _reprs2,
                                                       num_of_classes2,
                                                       w_sign_indices, wo_sign_indices_array);
      printf("%d\n", num_of_reprs2);
      for (int j = 0; j < num_of_reprs2; j++)
        {
          mpz_add(res, res, one);
        }
    }
  printf("%s\n", mpz_get_str(NULL, 10, res));
  /* num = 3 => res = 3222705
     num = 2 => res = 98149 */
  mpz_clear(res);
}

void test_set_idices_2(void)
{
  {
    Rk16VecInt vec1[16] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    Rk16VecInt vec2[16] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0};
    int indices_array[8][16] = {0};
    set_wo_sign_indices_array2(indices_array, vec1, vec2);
    for (int i = 0; indices_array[i][0]; i++)
      {
        print_vec(indices_array[i], 16);
      }
  }
  {
    Rk16VecInt vec1[16] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    Rk16VecInt vec2[16] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 0, 2, 2, 0};
    int indices_array[8][16] = {0};
    set_wo_sign_indices_array2(indices_array, vec1, vec2);
    for (int i = 0; indices_array[i][0]; i++)
      {
        print_vec(indices_array[i], 16);
      }
  }
}

int main()
{
  /* printf("test_norm_vec_rk16\n"); */
  /* if (test_norm_vec_rk16()) */
  /*   { */
  /*     printf("OK\n"); */
  /*   } */
  /* printf("test_repr_rk16\n"); */
  /* if (test_repr_rk16()) */
  /*   { */
  /*     printf("OK\n"); */
  /*   } */
  /* printf("test_normalize_vec_rk16_w_indices\n"); */
  /* if (test_normalize_vec_rk16_w_indices()) */
  /*   { */
  /*     printf("Ok\n"); */
  /*   } */

  test_repr_rk16_w_idices();

  return 0;
}


/* Local Variables: */
/* compile-command: "gcc test_representatives.c -o test_representatives.o -O3 -le8vectors -lrank16_vectors -lvector_utils -lm -lmpir -L../lib -std=c11 -Wall -Wextra" */
/* End: */
