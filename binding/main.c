#include <fmpz.h>
#include <math.h>
#include "e8vectors.h"

static inline int inner_prod(int s[8], int t[8])
{
  return ((2*s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7]) * t[0] +
          (s[0] + 2*s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + 2*s[7]) * t[1] +
          (s[0] + s[1] + 2*s[2] + s[3] + s[4] + s[5] + s[6] + 2*s[7]) * t[2] +
          (s[0] + s[1] + s[2] + 2*s[3] + s[4] + s[5] + s[6] + 2*s[7]) * t[3] +
          (s[0] + s[1] + s[2] + s[3] + 2*s[4] + s[5] + s[6] + 2*s[7]) * t[4] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + 2*s[5] + s[6] + 2*s[7]) * t[5] +
          (s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + 2*s[6] + 2*s[7]) * t[6] +
          (s[0] + 2*s[1] + 2*s[2] + 2*s[3] + 2*s[4] + 2*s[5] + 2*s[6] + 4*s[7]) * t[7]);
}


void miyawaki_theta(int a, int b, int c, int d, int e, int f)
{
  int vs1[MAX_NM_OF_VECTORS][8];
  int vs2[MAX_NM_OF_VECTORS][8];
  int vs3[MAX_NM_OF_VECTORS][8];

  _set_vs3(vs1, vs2, vs3, a, b, c);

  fmpz_t res, rl_pt, im_pt;
  fmpz_t a0; fmpz_t a1; fmpz_t a2; fmpz_t a3; fmpz_t a4; fmpz_t a5; fmpz_t a6; fmpz_t a7; fmpz_t a8;
  fmpz_t tmp;

  fmpz_init(res); fmpz_init(rl_pt); fmpz_init(im_pt);
  fmpz_init(a0); fmpz_init(a1); fmpz_init(a2); fmpz_init(a3); fmpz_init(a4); fmpz_init(a5); fmpz_init(a6); fmpz_init(a7); fmpz_init(a8);
  fmpz_init(tmp);

  fmpz_zero(res);
  int s0, s1, s2, s3, s4, s5, s6, s7, t0, t1, t2, t3, t4, t5, t6, t7, u0, u1, u2, u3, u4, u5, u6, u7;
  int i, j, k;
  for (i = 0; i < num_of_vectors[a]; i++)
    {
      for (j = 0; j < num_of_vectors[b]; j++)
        {
          for (k = 0; k < num_of_vectors[c]; k++)
            {
              if (inner_prod(vs1[i], vs2[j]) == f)
                {
                  if (inner_prod(vs1[i], vs3[k]) == e)
                    {
                      if (inner_prod(vs2[j], vs3[k]) == d)
                        {
                          s0 = vs1[i][0];
                          s1 = vs1[i][1];
                          s2 = vs1[i][2];
                          s3 = vs1[i][3];
                          s4 = vs1[i][4];
                          s5 = vs1[i][5];
                          s6 = vs1[i][6];
                          s7 = vs1[i][7];

                          t0 = vs2[j][0];
                          t1 = vs2[j][1];
                          t2 = vs2[j][2];
                          t3 = vs2[j][3];
                          t4 = vs2[j][4];
                          t5 = vs2[j][5];
                          t6 = vs2[j][6];
                          t7 = vs2[j][7];

                          u0 = vs3[k][0];
                          u1 = vs3[k][1];
                          u2 = vs3[k][2];
                          u3 = vs3[k][3];
                          u4 = vs3[k][4];
                          u5 = vs3[k][5];
                          u6 = vs3[k][6];
                          u7 = vs3[k][7];

                          fmpz_set_si(rl_pt, (-s2 * t1 - s3 * t1 + s5 * t1 + s1 * t2 + s3 * t2 - s4 * t2 + s1 * t3 - s2 * t3 + s4 * t3 - s5 * t3 + s2 * t4 - s3 * t4 + s5 * t4 - s1 * t5 + s3 * t5 - s4 * t5) * u0 + (s2 * t0 + s3 * t0 - s5 * t0 - s0 * t2 - s0 * t3 - 2 * s5 * t3 + s0 * t5 + 2 * s3 * t5) * u1 + (-s1 * t0 - s3 * t0 + s4 * t0 + s0 * t1 + s0 * t3 + 2 * s4 * t3 - s0 * t4 - 2 * s3 * t4) * u2 + (-s1 * t0 + s2 * t0 - s4 * t0 + s5 * t0 + s0 * t1 + 2 * s5 * t1 - s0 * t2 - 2 * s4 * t2 + s0 * t4 + 2 * s2 * t4 - s0 * t5 - 2 * s1 * t5) * u3 + (-s2 * t0 + s3 * t0 - s5 * t0 + s0 * t2 + 2 * s3 * t2 - s0 * t3 - 2 * s2 * t3 + s0 * t5) * u4 + (s1 * t0 - s3 * t0 + s4 * t0 - s0 * t1 - 2 * s3 * t1 + s0 * t3 + 2 * s1 * t3 - s0 * t4) * u5);
                          fmpz_set_si(im_pt, (-s2 * t1 + s3 * t1 - s5 * t1 + s1 * t2 - s3 * t2 + s4 * t2 - s1 * t3 + s2 * t3 + s4 * t3 - s5 * t3 - s2 * t4 - s3 * t4 + s5 * t4 + s1 * t5 + s3 * t5 - s4 * t5) * u0 + (s2 * t0 - s3 * t0 + s5 * t0 - s0 * t2 - 2 * s3 * t2 + s0 * t3 + 2 * s2 * t3 - s0 * t5) * u1 + (-s1 * t0 + s3 * t0 - s4 * t0 + s0 * t1 + 2 * s3 * t1 - s0 * t3 - 2 * s1 * t3 + s0 * t4) * u2 + (s1 * t0 - s2 * t0 - s4 * t0 + s5 * t0 - s0 * t1 - 2 * s2 * t1 + s0 * t2 + 2 * s1 * t2 + s0 * t4 + 2 * s5 * t4 - s0 * t5 - 2 * s4 * t5) * u3 + (s2 * t0 + s3 * t0 - s5 * t0 - s0 * t2 - s0 * t3 - 2 * s5 * t3 + s0 * t5 + 2 * s3 * t5) * u4 + (-s1 * t0 - s3 * t0 + s4 * t0 + s0 * t1 + s0 * t3 + 2 * s4 * t3 - s0 * t4 - 2 * s3 * t4) * u5);

                          /* Computation of a0 = (rl_pt ** 8 - 28 * rl_pt ** 6 * im_pt ** 2 + 70 * rl_pt ** 4 * im_pt ** 4 - 28 * rl_pt ** 2
* im_pt ** 6 + im_pt ** 8)  */

                          fmpz_zero(tmp);
                          fmpz_pow_ui(a0, im_pt, 6);
                          fmpz_pow_ui(a1, rl_pt, 2);
                          fmpz_mul(a0, a0, a1);
                          fmpz_submul_ui(tmp, a0, 28);
                          fmpz_pow_ui(a0, rl_pt, 8);
                          fmpz_add(tmp, tmp, a0);
                          fmpz_pow_ui(a0, im_pt, 2);
                          fmpz_pow_ui(a1, rl_pt, 6);
                          fmpz_mul(a0, a0, a1);
                          fmpz_submul_ui(tmp, a0, 28);
                          fmpz_pow_ui(a0, im_pt, 4);
                          fmpz_pow_ui(a1, rl_pt, 4);
                          fmpz_mul(a0, a0, a1);
                          fmpz_addmul_ui(tmp, a0, 70);
                          fmpz_pow_ui(a0, im_pt, 8);
                          fmpz_add(tmp, tmp, a0);

                          fmpz_add(res, res, tmp);
                        }
                    }
                }
            }
        }
    }

  fmpz_print(res);

  fmpz_clear(a0); fmpz_clear(a1); fmpz_clear(a2); fmpz_clear(a3); fmpz_clear(a4); fmpz_clear(a5); fmpz_clear(a6); fmpz_clear(a7); fmpz_clear(a8);
  fmpz_clear(rl_pt); fmpz_clear(im_pt); fmpz_clear(res); fmpz_clear(tmp);

}

int main(void)
{
  _cache_vectors();

  printf("Computation of cached vectors done.\n");

  /* miyawaki_theta(vs1, vs1, vs1, 1, 1, 1, 0, 0, 0); */
  /* printf("\n"); */
  miyawaki_theta(1, 1, 1, 1, 1, 1);
  printf("\n");
  miyawaki_theta(2, 1, 1, 1, 1, 1);
  printf("\n");
  miyawaki_theta(2, 2, 1, 1, 1, 1);

  printf("\n");
  return 0;
}

/* Local Variables: */
/* compile-command: "cd ..; make compile" */
/* End: */
