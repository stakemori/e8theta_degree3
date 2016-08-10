#include "e8vectors.h"
#include "rank16_vectors.h"
#include <stdio.h>
#include "memory.h"
#include <math.h>
#include "vector_utils.h"
#include "mpir.h"

Rk16VecInt inner_prod_rk16(Rk16VecInt s[16], Rk16VecInt t[16])
{
  int a = (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] +
           s[8] + s[9] + s[10] + s[11] + s[12] + s[13] + s[14] + 2 * s[15]);
  return ((3 * s[0] + a - s[15]) * t[0] + (a + s[1]) * t[1] + (a + s[2]) * t[2] +(a + s[3]) * t[3] +
          (a + s[4]) * t[4] + (a + s[5]) * t[5] + (a + s[6]) * t[6] + (a + s[7]) * t[7] +
          (a + s[8]) * t[8] + (a + s[9]) * t[9] + (a + s[10]) * t[10] + (a + s[11]) * t[11] +
          (a + s[12]) * t[12] + (a + s[13]) * t[13] + (a + s[14]) * t[14] +
          (-s[0] + 2 * a) * t[15]);
}

static void print_vec(int * vec, int a)
{
  for (int i = 0; i < a; i++)
    {
      printf("%d, ", vec[i]);
    }
  printf("\n");
}

static int test_norm_vec_rk16(void)
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
      _convert_to_euclid_vector_rk16(vec1);
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
          bl = bl & (vec1[j] == vec2[j]);
        }
      if (! bl)
        {
          bl = 1;
          for (int j = 0; j < 7; j++)
            {
              bl = bl & (vec1[j] == -vec2[j]);
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

static int test_repr_rk16(void)
{
  static Rk16VecInt _reprs_rk16[MAX_NM_REPRS_RK16][16];
  static int num_of_reprs_rk16[MAX_NM_REPRS_RK16];
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
      _convert_to_euclid_vector_rk16(vec);
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
  /*     _convert_to_euclid_vector_rk16(vec); */
  /*     print_vec(vec); */
  /*     normalize_vec_rk16_last9(vec); */
  /*     print_vec(vec); */
  /*   } */

}


void number_of_loops(int a, int b, int c, int d, int e, int f)
{
  cache_vectors_rk16();
  static Rk16VecInt _reprs[MAX_NM_REPRS_RK16][16];
  static int num_of_classes[MAX_NM_REPRS_RK16];
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
  /* number_of_loops(1, 1, 1, 0, 0, 0); */
  /* number_of_loops(1, 1, 3, 1, 1, 1); */
  /* number_of_loops(2, 2, 2, 1, 1, 1); */
  /* number_of_loops(2, 2, 2, 0, 0, 0); */
  /* number_of_loops(1, 3, 3, 1, 0, 0); */
  return 0;
}


/* Local Variables: */
/* compile-command: "gcc test_representatives.c -o test_representatives.o -le8vectors -lrank16_vectors -lvector_utils -lm -lmpir -L../lib -std=c11" */
/* End: */
