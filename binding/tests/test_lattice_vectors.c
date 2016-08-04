#include "e8vectors.h"
#include "rank16_vectors.h"
#include <stdio.h>

static inline short norm_rk16(short s[16])
{
  return ((4*s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9] + s[10] + s[11] + s[12] + s[13] + s[14] + s[15])*s[0] +
          (s[0] + 2*s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9] + s[10] + s[11] + s[12] + s[13] + s[14] + 2*s[15])*s[1] +
          (s[0] + s[1] + 2*s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9] + s[10] + s[11] + s[12] + s[13] + s[14] + 2*s[15])*s[2] +
          (s[0] + s[1] + s[2] + 2*s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9] + s[10] + s[11] + s[12] + s[13] + s[14] + 2*s[15])*s[3] +
          (s[0] + s[1] + s[2] + s[3] + 2*s[4] + s[5] + s[6] + s[7] + s[8] + s[9] + s[10] + s[11] + s[12] + s[13] + s[14] + 2*s[15])*s[4] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + 2*s[5] + s[6] + s[7] + s[8] + s[9] + s[10] + s[11] + s[12] + s[13] + s[14] + 2*s[15])*s[5] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + 2*s[6] + s[7] + s[8] + s[9] + s[10] + s[11] + s[12] + s[13] + s[14] + 2*s[15])*s[6] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + 2*s[7] + s[8] + s[9] + s[10] + s[11] + s[12] + s[13] + s[14] + 2*s[15])*s[7] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + 2*s[8] + s[9] + s[10] + s[11] + s[12] + s[13] + s[14] + 2*s[15])*s[8] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + 2*s[9] + s[10] + s[11] + s[12] + s[13] + s[14] + 2*s[15])*s[9] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9] + 2*s[10] + s[11] + s[12] + s[13] + s[14] + 2*s[15])*s[10] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9] + s[10] + 2*s[11] + s[12] + s[13] + s[14] + 2*s[15])*s[11] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9] + s[10] + s[11] + 2*s[12] + s[13] + s[14] + 2*s[15])*s[12] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9] + s[10] + s[11] + s[12] + 2*s[13] + s[14] + 2*s[15])*s[13] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9] + s[10] + s[11] + s[12] + s[13] + 2*s[14] + 2*s[15])*s[14] +
          (s[0] + 2*s[1] + 2*s[2] + 2*s[3] + 2*s[4] + 2*s[5] + 2*s[6] + 2*s[7] + 2*s[8] + 2*s[9] + 2*s[10] + 2*s[11] + 2*s[12] + 2*s[13] + 2*s[14] + 4*s[15])*s[15]);
}


static inline int norm(int s[8])
{
  return ((2*s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7])*s[0] +
          (s[0] + 2*s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + 2*s[7])*s[1] +
          (s[0] + s[1] + 2*s[2] + s[3] + s[4] + s[5] + s[6] + 2*s[7])*s[2] +
          (s[0] + s[1] + s[2] + 2*s[3] + s[4] + s[5] + s[6] + 2*s[7])*s[3] +
          (s[0] + s[1] + s[2] + s[3] + 2*s[4] + s[5] + s[6] + 2*s[7])*s[4] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + 2*s[5] + s[6] + 2*s[7])*s[5] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + 2*s[6] + 2*s[7])*s[6] +
          (s[0] + 2*s[1] + 2*s[2] + 2*s[3] + 2*s[4] + 2*s[5] + 2*s[6] + 4*s[7])*s[7]);
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
/* compile-command: "gcc test_lattice_vectors.c -o test_lattice_vectors -le8vectors -lrank16_vectors -lm -L../lib -std=c11" */
/* End: */
