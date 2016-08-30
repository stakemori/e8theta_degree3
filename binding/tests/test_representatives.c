#include "e8vectors.h"
#include "rank16_vectors.h"
#include <stdio.h>
#include "memory.h"
#include <math.h>
#include "vector_utils.h"
#include "mpir.h"

int inner_prod(int s[8], int t[8])
{
  int sum = 0;
  for (int i = 0; i < 8; i++)
    {
      sum += s[i] * t[i];
    }
  return sum >> 2;
}


int inner_prod_rk16(int s[16], int t[16])
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
  int vec1[16];
  int vec2[16];
  int vec3[16];
  int vec4[16];
  int bl = 1;
  int * cached_vec = cached_vectors_rk16_ptr[2];
  for (int i = 0; i < num_of_vectors_rk16[2]; i++, cached_vec += 16)
    {
      memcpy(vec1, cached_vec, sizeof(int) * 16);
      memcpy(vec2, vec1, sizeof(int) * 16);
      normalize_vec_last_len_elts(vec1, 16, 9);
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
  static int _reprs_rk16[MAX_NM_REPRS_RK16][16];
  static unsigned int num_of_reprs_rk16[MAX_NM_REPRS_RK16];
  int n = 2;
  int num;
  int vec[16];
  for (int i = 0; i < n + 1; i++)
    {
      num = repr_modulo_autom_rk16(i, _reprs_rk16, num_of_reprs_rk16);
      printf("%d\n", num);
    }
  int s = 0;
  for (int i = 0; i < num; i++)
    {
      s += num_of_reprs_rk16[i];
      memcpy(vec, _reprs_rk16[i], sizeof(int) * 16);
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
  /* int vec[16]; */
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


int vec_equal(int * vec1, int * vec2, int start, int end)
{
  int is_equal = 1;
  for (int i = start; i < end; i++)
    {
      is_equal = is_equal && (vec1[i] == vec2[i]);
    }
  return is_equal;
}

void abs_vec(int * vec, int len)
{
  for (int i = 0; i < len; i++)
    {
      vec[i] = abs(vec[i]);
    }
}

int equal_as_multiset(int vec1[16], int vec2[16], int start, int end)
{
  int _vec1[16];
  int _vec2[16];
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
  int vec[16];
  int vec1[16];
  int num = 2;
  int * cached_vec = cached_vectors_rk16_ptr[num];
  for (int i = 0; i < num_of_vectors_rk16[num]; i++, cached_vec += 16)
    {
      memcpy(vec, cached_vec, 16 * sizeof(int));
      memcpy(vec1, vec, 16 * sizeof(int));
      normalize_vec_w_indices(vec, zero_idcs, non_zero_idcss);
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
  set_w_sign_indices(idx_array, vec, 16, 9);
  set_wo_sign_indices_array(wo_sign_indices_array, vec, 16, 9);
  normalize_vec_w_indices(vec1, idx_array, wo_sign_indices_array);
  print_vec(vec1, 16);
  for (int j = 0; wo_sign_indices_array[j][0]; j++)
    {
      print_vec(wo_sign_indices_array[j], 16);
    }
}

void test_repr_rk16_w_idices(int a, int b, int c, int d, int e, int f)
{
  cache_vectors_rk16();
  static int reprs_k[MAX_NM_REPRS_RK16][16];
  static unsigned int num_of_classes_k[MAX_NM_REPRS_RK16];

  static int reprs_j[MAX_NM_REPRS_RK16][16];
  static unsigned int num_of_classes_j[MAX_NM_REPRS_RK16];
  static int vecs_j[MAX_NM_OF_VECTORS_RK16][16];
  int num_of_vecs_j;

  static int reprs_i[MAX_NM_REPRS_RK16][16];
  static unsigned int num_of_classes_i[MAX_NM_REPRS_RK16];
  static int vecs_i[MAX_NM_OF_VECTORS_RK16][16];

  int num_of_reprs_k = repr_modulo_autom_rk16(c, reprs_k, num_of_classes_k);

  mpz_t res;
  mpz_t one;
  mpz_init(res);
  mpz_init(one);
  mpz_set_si(one, 1);
  printf("%d\n", num_of_reprs_k);
  int num_of_reprs_j;
  int i_reprs_max = 100;
  for (int k = 0; k < num_of_reprs_k; k += 6)
    {
      int wo_sign_indices_array[8][16] = {0};
      int w_sign_indices[16] = {0};
      printf("k: %d\n", k);
      print_vec(reprs_k[k], 16);
      set_w_sign_indices(w_sign_indices, reprs_k[k], 16, 9);
      set_wo_sign_indices_array(wo_sign_indices_array, reprs_k[k], 16, 9);
      print_vec(w_sign_indices, 16);

      num_of_vecs_j = 0;
      int * cached_vec_b = cached_vectors_rk16_ptr[b];
      for (int l = 0; l < num_of_vectors_rk16[b]; l++, cached_vec_b += 16)
        {
          if (inner_prod_rk16(cached_vec_b, reprs_k[k]) == d)
            {
              memcpy(vecs_j[num_of_vecs_j++], cached_vec_b, sizeof(int) * 16);
            }
        }
      num_of_reprs_j = repr_modulo_autom_rk16_w_indices(vecs_j, num_of_vecs_j, reprs_j,
                                                        num_of_classes_j,
                                                        w_sign_indices, wo_sign_indices_array);
      printf("%d\n", num_of_reprs_j);
      for (int j = 0; j < num_of_reprs_j; j++)
        {
          int wo_sign_indices_array[8][16] = {0};
          int w_sign_indices[16] = {0};
          set_wo_sign_indices_array2(wo_sign_indices_array, reprs_j[j], reprs_k[k], 16, 9);
          set_w_sign_indices_2(w_sign_indices, reprs_j[j], reprs_k[k], 16, 9);

          int num_of_vecs_i = 0;
          int * cached_vec_a = cached_vectors_rk16_ptr[a];
          for (int l = 0; l < num_of_vectors_rk16[a]; l++, cached_vec_a += 16)
            {
              if ((inner_prod_rk16(cached_vec_a, reprs_k[k]) == e) &&
                  (inner_prod_rk16(cached_vec_a, reprs_j[j]) == f))
                {
                  memcpy(vecs_i[num_of_vecs_i++], cached_vec_a, sizeof(int) * 16);
                }
            }
          int num_of_reprs_i = repr_modulo_autom_rk16_w_indices(vecs_i, num_of_vecs_i, reprs_i,
                                                                num_of_classes_i,
                                                                w_sign_indices,
                                                                wo_sign_indices_array);
          if (num_of_reprs_i > i_reprs_max)
            {
              i_reprs_max = num_of_reprs_i;
            }

          for (int i = 0; i < num_of_reprs_i; i++)
            {
              mpz_add(res, res, one);
            }
        }
    }
  printf("ireprs_max:%d\n", i_reprs_max);
  printf("%s\n", mpz_get_str(NULL, 10, res));
  /* num = 3 => res = 3222705 */
  /* num = 2 => res = 98149 */
  mpz_clear(res);
}

void test_set_idices_2(void)
{
  {
    int vec1[16] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int vec2[16] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0};
    int indices_array[8][16] = {0};
    set_wo_sign_indices_array2(indices_array, vec1, vec2, 16, 9);
    for (int i = 0; indices_array[i][0]; i++)
      {
        print_vec(indices_array[i], 16);
      }
  }
  {
    int vec1[16] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int vec2[16] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 0, 2, 2, 0};
    int indices_array[8][16] = {0};
    set_wo_sign_indices_array2(indices_array, vec1, vec2, 16, 9);
    for (int i = 0; indices_array[i][0]; i++)
      {
        print_vec(indices_array[i], 16);
      }
  }
  {
    int vec1[8] = {0, 0, 0, 0, 0, 0, 1, 1};
    int vec2[16] = {0, 0, 0, 0, 1, 1, 1, 1};
    int indices_array[8][16] = {0};
    set_wo_sign_indices_array2(indices_array, vec1, vec2, 8, 2);
    for (int i = 0; indices_array[i][0]; i++)
      {
        print_vec(indices_array[i], 16);
      }

  }
}

void test_repr_e8_w_indices(int a, int b, int c, int d, int e, int f)
{
  cache_vectors();
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

  mpz_t res;
  mpz_t one;
  mpz_init(res);
  mpz_init(one);
  mpz_set_si(one, 1);
  printf("%d\n", num_of_reprs_k);
  int num_of_reprs_j;
  for (int k = 0; k < num_of_reprs_k; k += 6)
    {
      int wo_sign_indices_array[8][16] = {0};
      int w_sign_indices[8] = {0};
      printf("k: %d\n", k);
      print_vec(reprs_k[k], 8);
      set_w_sign_indices(w_sign_indices, reprs_k[k], 8, 2);
      set_wo_sign_indices_array(wo_sign_indices_array, reprs_k[k], 8, 2);
      print_vec(w_sign_indices, 8);

      num_of_vecs_j = 0;
      int * cached_vec_b = cached_vectors_ptr[b];
      for (int l = 0; l < num_of_vectors[b]; l++, cached_vec_b += 8)
        {
          if (inner_prod(cached_vec_b, reprs_k[k]) == d)
            {
              memcpy(vecs_j[num_of_vecs_j++], cached_vec_b, sizeof(int) * 8);
            }
        }
      num_of_reprs_j = repr_modulo_autom_w_indices(vecs_j, num_of_vecs_j, reprs_j,
                                                   num_of_classes_j,
                                                   w_sign_indices, wo_sign_indices_array);
      printf("%d\n", num_of_reprs_j);
      for (int j = 0; j < num_of_reprs_j; j++)
        {
          int wo_sign_indices_array[8][16] = {0};
          int w_sign_indices[8] = {0};
          set_wo_sign_indices_array2(wo_sign_indices_array, reprs_j[j], reprs_k[k], 8, 2);
          set_w_sign_indices_2(w_sign_indices, reprs_j[j], reprs_k[k], 8, 2);

          int num_of_vecs_i = 0;
          int * cached_vec_a = cached_vectors_ptr[a];
          for (int l = 0; l < num_of_vectors_rk16[a]; l++, cached_vec_a += 8)
            {
              if ((inner_prod_rk16(cached_vec_a, reprs_k[k]) == e) &&
                  (inner_prod_rk16(cached_vec_a, reprs_j[j]) == f))
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
              mpz_add(res, res, one);
            }
        }
    }
}

int main()
{
  printf("test_norm_vec_rk16\n");
  if (test_norm_vec_rk16())
    {
      printf("OK\n");
    }
  printf("test_repr_rk16\n");
  if (test_repr_rk16())
    {
      printf("OK\n");
    }
  printf("test_normalize_vec_rk16_w_indices\n");
  if (test_normalize_vec_rk16_w_indices())
    {
      printf("Ok\n");
    }

  /* test_repr_rk16_w_idices(1, 1, 2, 0, 0, 0); */
  /* test_repr_rk16_w_idices(3, 3, 3, -2, 2, 2); */
  /* test_repr_rk16_w_idices(1, 4, 4, 0, 0, 0); */
  /* test_repr_rk16_w_idices(2, 2, 4, 0, 0, 0); */
  /* test_repr_rk16_w_idices(1, 3, 3, -2, 2, 2); */
  /* test_repr_rk16_w_idices(1, 1, 2, 0, 0, 0); */

  /* test_set_idices_2(); */
  /* test_repr_e8_w_indices(1, 3, 3, 1, 1, 1); */
  return 0;
}


/* Local Variables: */
/* compile-command: "gcc test_representatives.c -o test_representatives.o -O3 -le8vectors -lrank16_vectors -lvector_utils -lm -lmpir -L../lib -std=c11 -Wall -Wextra" */
/* End: */
