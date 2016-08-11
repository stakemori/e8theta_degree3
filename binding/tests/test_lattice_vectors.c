#include "e8vectors.h"
#include "rank16_vectors.h"
#include <stdio.h>

static inline Rk16VecInt norm_rk16(Rk16VecInt s[16])
{
  int sum = 0;
  for (int i = 0; i < 16; i++)
    {
      sum += s[i] * s[i];
    }
  return sum >> 2;
}


/* static void print_vec(int * vec, int a) */
/* { */
/*   for (int i = 0; i < a; i++) */
/*     { */
/*       printf("%d, ", vec[i]); */
/*     } */
/*   printf("\n"); */
/* } */


static inline int norm(int s[8])
{
  int sum = 0;
  for (int i = 0; i < 8; i++)
    {
      sum += s[i] * s[i];
    }
  return sum >> 2;
}


int test_e8vectors(void)
{
  for (int i = 1; i < MAX_NORM + 1; i++)
    {
      int count = 0;
      cache_vectors();
      for (int j = 0; j < MAX_NM_OF_VECTORS; j++)
        {
          if (norm(cached_vectors[i][j]) == i * 2)
            {
              count++;
            }
        }
      if (!(count == num_of_vectors[i]))
        {
          printf("i: %d, count %d\n", i, count);
          return 0;
        }
    }
  return 1;
}

int test_rank16_vectors(void)
{
  cache_vectors_rk16();
  for (int i = 1; i < MAX_NORM_RK16 + 1; i++)
    {
      int count = 0;
      for (int j = 0; j < MAX_NM_OF_VECTORS_RK16; j++)
        {
          if (norm_rk16(cached_vectors_rk16[i][j]) >> 1 == i)
            {
              count++;
            }
        }
      if (!(count == num_of_vectors_rk16[i]))
        {
          return 0;
        }
    }
  return 1;
}

int main(void)
{
  if (test_e8vectors() & test_rank16_vectors())
    {
      printf("OK\n");
    }
  return 0;
}

/* Local Variables: */
/* compile-command: "gcc test_lattice_vectors.c -o test_lattice_vectors.o -le8vectors -lrank16_vectors -lm -L../lib -std=c11" */
/* End: */
