#include "e8vectors.h"

inline int inner_prod(int s[8], int t[8])
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

char * miyawaki_theta_c(int a, int b, int c, int d, int e, int f)
{
  cache_vectors();
  /* Use static to avoid segmentation fault */
  static int vs1[MAX_NM_OF_VECTORS][8];
  static int vs2[MAX_NM_OF_VECTORS][8];
  static int vs3[MAX_NM_OF_VECTORS][8];

  _set_vs3(vs1, vs2, vs3, a, b, c);

  fmpz_t res, rl_pt, im_pt;
  fmpz_t a0; fmpz_t a1; fmpz_t a2; fmpz_t a3; fmpz_t a4; fmpz_t a5; fmpz_t a6; fmpz_t a7; fmpz_t a8;
  fmpz_t tmp;

  fmpz_t s0, s1, s2, s3, s4, s5, s6, s7;
  fmpz_t t0, t1, t2, t3, t4, t5, t6, t7;
  fmpz_t u0, u1, u2, u3, u4, u5, u6, u7;

  fmpz_init(res); fmpz_init(rl_pt); fmpz_init(im_pt);
  fmpz_init(a0); fmpz_init(a1); fmpz_init(a2); fmpz_init(a3); fmpz_init(a4); fmpz_init(a5); fmpz_init(a6); fmpz_init(a7); fmpz_init(a8);
  fmpz_init(tmp);

  fmpz_init(s0); fmpz_init(s1); fmpz_init(s2); fmpz_init(s3); fmpz_init(s4); fmpz_init(s5); fmpz_init(s6); fmpz_init(s7);
  fmpz_init(t0); fmpz_init(t1); fmpz_init(t2); fmpz_init(t3); fmpz_init(t4); fmpz_init(t5); fmpz_init(t6); fmpz_init(t7);
  fmpz_init(u0); fmpz_init(u1); fmpz_init(u2); fmpz_init(u3); fmpz_init(u4); fmpz_init(u5); fmpz_init(u6); fmpz_init(u7);


  fmpz_zero(res);
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
                          fmpz_set_si(s0, vs1[i][0]);
                          fmpz_set_si(s1, vs1[i][1]);
                          fmpz_set_si(s2, vs1[i][2]);
                          fmpz_set_si(s3, vs1[i][3]);
                          fmpz_set_si(s4, vs1[i][4]);
                          fmpz_set_si(s5, vs1[i][5]);
                          fmpz_set_si(s6, vs1[i][6]);
                          fmpz_set_si(s7, vs1[i][7]);

                          fmpz_set_si(t0, vs2[j][0]);
                          fmpz_set_si(t1, vs2[j][1]);
                          fmpz_set_si(t2, vs2[j][2]);
                          fmpz_set_si(t3, vs2[j][3]);
                          fmpz_set_si(t4, vs2[j][4]);
                          fmpz_set_si(t5, vs2[j][5]);
                          fmpz_set_si(t6, vs2[j][6]);
                          fmpz_set_si(t7, vs2[j][7]);

                          fmpz_set_si(u0, vs3[k][0]);
                          fmpz_set_si(u1, vs3[k][1]);
                          fmpz_set_si(u2, vs3[k][2]);
                          fmpz_set_si(u3, vs3[k][3]);
                          fmpz_set_si(u4, vs3[k][4]);
                          fmpz_set_si(u5, vs3[k][5]);
                          fmpz_set_si(u6, vs3[k][6]);
                          fmpz_set_si(u7, vs3[k][7]);

                          fmpz_set(a0, s4);
                          fmpz_add(a0, a0, s1);
                          fmpz_sub(a0, a0, s2);
                          fmpz_sub(a0, a0, s5);
                          fmpz_set(a1, t5);
                          fmpz_add(a1, a1, t2);
                          fmpz_sub(a1, a1, t4);
                          fmpz_sub(a1, a1, t1);
                          fmpz_neg(a2, s4);
                          fmpz_sub(a2, a2, s1);
                          fmpz_set(a3, s5);
                          fmpz_add(a3, a3, s2);
                          fmpz_set(a4, s1);
                          fmpz_sub(a4, a4, s4);
                          fmpz_neg(a5, s2);
                          fmpz_add(a5, a5, s5);
                          fmpz_mul(a5, a5, t1);
                          fmpz_addmul(a5, a4, t2);
                          fmpz_addmul(a5, a3, t4);
                          fmpz_addmul(a5, a2, t5);
                          fmpz_addmul(a5, a1, s3);
                          fmpz_addmul(a5, a0, t3);
                          fmpz_neg(a0, s4);
                          fmpz_add(a0, a0, s2);
                          fmpz_add(a0, a0, s5);
                          fmpz_sub(a0, a0, s1);
                          fmpz_set(a1, u4);
                          fmpz_sub(a1, a1, u5);
                          fmpz_add(a1, a1, u1);
                          fmpz_sub(a1, a1, u2);
                          fmpz_set(a2, s4);
                          fmpz_add(a2, a2, s1);
                          fmpz_neg(a3, s5);
                          fmpz_sub(a3, a3, s2);
                          fmpz_neg(a4, s1);
                          fmpz_add(a4, a4, s4);
                          fmpz_set(a6, s2);
                          fmpz_sub(a6, a6, s5);
                          fmpz_mul(a6, a6, u1);
                          fmpz_addmul(a6, a4, u2);
                          fmpz_addmul(a6, a3, u4);
                          fmpz_addmul(a6, a2, u5);
                          fmpz_addmul(a6, a1, s3);
                          fmpz_addmul(a6, a0, u3);
                          fmpz_neg(a0, t5);
                          fmpz_add(a0, a0, t4);
                          fmpz_add(a0, a0, t1);
                          fmpz_sub(a0, a0, t2);
                          fmpz_neg(a1, u4);
                          fmpz_sub(a1, a1, u1);
                          fmpz_add(a1, a1, u5);
                          fmpz_add(a1, a1, u2);
                          fmpz_neg(a2, t4);
                          fmpz_sub(a2, a2, t1);
                          fmpz_set(a3, t2);
                          fmpz_add(a3, a3, t5);
                          fmpz_set(a4, t1);
                          fmpz_sub(a4, a4, t4);
                          fmpz_neg(a7, t2);
                          fmpz_add(a7, a7, t5);
                          fmpz_mul(a7, a7, u1);
                          fmpz_addmul(a7, a4, u2);
                          fmpz_addmul(a7, a3, u4);
                          fmpz_addmul(a7, a2, u5);
                          fmpz_addmul(a7, a1, t3);
                          fmpz_addmul(a7, a0, u3);
                          fmpz_mul_2exp(a0, s1, 1); fmpz_neg(a0, a0);
                          fmpz_mul_2exp(a1, s2, 1);
                          fmpz_mul_2exp(a2, s4, 1); fmpz_neg(a2, a2);
                          fmpz_mul_2exp(a3, s5, 1);
                          fmpz_mul(a3, a3, t1);
                          fmpz_addmul(a3, a2, t2);
                          fmpz_addmul(a3, a1, t4);
                          fmpz_addmul(a3, a0, t5);
                          fmpz_mul_2exp(a0, s1, 1);
                          fmpz_mul_2exp(a1, s2, 1); fmpz_neg(a1, a1);
                          fmpz_mul_2exp(a2, s4, 1);
                          fmpz_mul_2exp(a4, s5, 1); fmpz_neg(a4, a4);
                          fmpz_mul(a4, a4, u1);
                          fmpz_addmul(a4, a2, u2);
                          fmpz_addmul(a4, a1, u4);
                          fmpz_addmul(a4, a0, u5);
                          fmpz_mul_2exp(a0, t1, 1); fmpz_neg(a0, a0);
                          fmpz_mul_2exp(a1, t2, 1);
                          fmpz_mul_2exp(a2, t4, 1); fmpz_neg(a2, a2);
                          fmpz_mul_2exp(a8, t5, 1);
                          fmpz_mul(a8, a8, u1);
                          fmpz_addmul(a8, a2, u2);
                          fmpz_addmul(a8, a1, u4);
                          fmpz_addmul(a8, a0, u5);
                          fmpz_mul(a8, a8, s3);
                          fmpz_addmul(a8, a4, t3);
                          fmpz_addmul(a8, a3, u3);
                          fmpz_addmul(a8, a7, s0);
                          fmpz_addmul(a8, a6, t0);
                          fmpz_addmul(a8, a5, u0);
                          fmpz_set(rl_pt, a8);

                          fmpz_set(a0, s4);
                          fmpz_sub(a0, a0, s1);
                          fmpz_add(a0, a0, s2);
                          fmpz_sub(a0, a0, s5);
                          fmpz_set(a1, t5);
                          fmpz_sub(a1, a1, t2);
                          fmpz_sub(a1, a1, t4);
                          fmpz_add(a1, a1, t1);
                          fmpz_neg(a2, s4);
                          fmpz_add(a2, a2, s1);
                          fmpz_set(a3, s5);
                          fmpz_sub(a3, a3, s2);
                          fmpz_set(a4, s1);
                          fmpz_add(a4, a4, s4);
                          fmpz_neg(a5, s2);
                          fmpz_sub(a5, a5, s5);
                          fmpz_mul(a5, a5, t1);
                          fmpz_addmul(a5, a4, t2);
                          fmpz_addmul(a5, a3, t4);
                          fmpz_addmul(a5, a2, t5);
                          fmpz_addmul(a5, a1, s3);
                          fmpz_addmul(a5, a0, t3);
                          fmpz_neg(a0, s4);
                          fmpz_sub(a0, a0, s2);
                          fmpz_add(a0, a0, s5);
                          fmpz_add(a0, a0, s1);
                          fmpz_set(a1, u4);
                          fmpz_sub(a1, a1, u5);
                          fmpz_sub(a1, a1, u1);
                          fmpz_add(a1, a1, u2);
                          fmpz_set(a2, s4);
                          fmpz_sub(a2, a2, s1);
                          fmpz_neg(a3, s5);
                          fmpz_add(a3, a3, s2);
                          fmpz_neg(a4, s1);
                          fmpz_sub(a4, a4, s4);
                          fmpz_set(a6, s2);
                          fmpz_add(a6, a6, s5);
                          fmpz_mul(a6, a6, u1);
                          fmpz_addmul(a6, a4, u2);
                          fmpz_addmul(a6, a3, u4);
                          fmpz_addmul(a6, a2, u5);
                          fmpz_addmul(a6, a1, s3);
                          fmpz_addmul(a6, a0, u3);
                          fmpz_neg(a0, t5);
                          fmpz_add(a0, a0, t4);
                          fmpz_sub(a0, a0, t1);
                          fmpz_add(a0, a0, t2);
                          fmpz_neg(a1, u4);
                          fmpz_add(a1, a1, u1);
                          fmpz_add(a1, a1, u5);
                          fmpz_sub(a1, a1, u2);
                          fmpz_neg(a2, t4);
                          fmpz_add(a2, a2, t1);
                          fmpz_neg(a3, t2);
                          fmpz_add(a3, a3, t5);
                          fmpz_set(a4, t1);
                          fmpz_add(a4, a4, t4);
                          fmpz_neg(a7, t2);
                          fmpz_sub(a7, a7, t5);
                          fmpz_mul(a7, a7, u1);
                          fmpz_addmul(a7, a4, u2);
                          fmpz_addmul(a7, a3, u4);
                          fmpz_addmul(a7, a2, u5);
                          fmpz_addmul(a7, a1, t3);
                          fmpz_addmul(a7, a0, u3);
                          fmpz_mul_2exp(a0, s4, 1); fmpz_neg(a0, a0);
                          fmpz_mul_2exp(a1, s5, 1);
                          fmpz_mul_2exp(a2, s1, 1);
                          fmpz_mul_2exp(a3, s2, 1); fmpz_neg(a3, a3);
                          fmpz_mul(a3, a3, t1);
                          fmpz_addmul(a3, a2, t2);
                          fmpz_addmul(a3, a1, t4);
                          fmpz_addmul(a3, a0, t5);
                          fmpz_mul_2exp(a0, s4, 1);
                          fmpz_mul_2exp(a1, s5, 1); fmpz_neg(a1, a1);
                          fmpz_mul_2exp(a2, s1, 1); fmpz_neg(a2, a2);
                          fmpz_mul_2exp(a4, s2, 1);
                          fmpz_mul(a4, a4, u1);
                          fmpz_addmul(a4, a2, u2);
                          fmpz_addmul(a4, a1, u4);
                          fmpz_addmul(a4, a0, u5);
                          fmpz_mul_2exp(a0, t4, 1); fmpz_neg(a0, a0);
                          fmpz_mul_2exp(a1, t5, 1);
                          fmpz_mul_2exp(a2, t1, 1);
                          fmpz_mul_2exp(a8, t2, 1); fmpz_neg(a8, a8);
                          fmpz_mul(a8, a8, u1);
                          fmpz_addmul(a8, a2, u2);
                          fmpz_addmul(a8, a1, u4);
                          fmpz_addmul(a8, a0, u5);
                          fmpz_mul(a8, a8, s3);
                          fmpz_addmul(a8, a4, t3);
                          fmpz_addmul(a8, a3, u3);
                          fmpz_addmul(a8, a7, s0);
                          fmpz_addmul(a8, a6, t0);
                          fmpz_addmul(a8, a5, u0);
                          fmpz_set(im_pt, a8);

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

  char *res_str = _store_fmpz_str_using_malloc(res);

  fmpz_clear(a0); fmpz_clear(a1); fmpz_clear(a2); fmpz_clear(a3); fmpz_clear(a4); fmpz_clear(a5); fmpz_clear(a6); fmpz_clear(a7); fmpz_clear(a8);
  fmpz_clear(rl_pt); fmpz_clear(im_pt); fmpz_clear(res); fmpz_clear(tmp);

  fmpz_clear(s0); fmpz_clear(s1); fmpz_clear(s2); fmpz_clear(s3); fmpz_clear(s4); fmpz_clear(s5); fmpz_clear(s6); fmpz_clear(s7); fmpz_clear(t0); fmpz_clear(t1); fmpz_clear(t2); fmpz_clear(t3); fmpz_clear(t4); fmpz_clear(t5); fmpz_clear(t6); fmpz_clear(t7); fmpz_clear(u0); fmpz_clear(u1); fmpz_clear(u2); fmpz_clear(u3); fmpz_clear(u4); fmpz_clear(u5); fmpz_clear(u6); fmpz_clear(u7);
  return res_str;
}

/* int main(void) */
/* { */
/*   char *s; */
/*   s = miyawaki_theta_c(1, 1, 1, 1, 1, 1); */
/*   printf("%s\n", s); */
/*   free(s); */
/*   s = NULL; */

/*   s = miyawaki_theta_c(2, 1, 1, 1, 1, 1); */
/*   printf("%s\n", s); */
/*   free(s); */
/*   s = NULL; */

/*   return 0; */
/* } */

/* Local Variables: */
/* compile-command: "cd ..; make compile-miyawaki-theta" */
/* End: */
