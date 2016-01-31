#include <fmpz.h>
#include <math.h>
#include "e8theta.h"
int inner_prod(int[8], int[8]);
int norm(int[8]);
void _cache_vectors(void);

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



inline int norm(int s[8])
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

/* Fourier coefficients of Eisenstein series of weight 4 */
int num_of_vectors[100] =
  {1, 240, 2160, 6720, 17520, 30240, 60480,
   82560, 140400, 181680, 272160, 319680, 490560,
   527520, 743040, 846720, 1123440, 1179360, 1635120,
   1646400, 2207520, 2311680, 2877120, 2920320, 3931200,
   3780240, 4747680, 4905600, 6026880, 5853600, 7620480,
   7150080, 8987760, 8951040, 10614240, 10402560,
   13262640, 12156960, 14817600, 14770560, 17690400,
   16541280, 20805120, 19081920, 23336640, 22891680,
   26282880, 24917760, 31456320, 28318320, 34022160,
   33022080, 38508960, 35730720, 44150400, 40279680,
   48297600, 46099200, 52682400, 49291200, 61810560,
   54475680, 64350720, 62497920, 71902320, 66467520,
   80559360, 72183360, 86093280, 81768960, 93623040,
   85898880, 106282800, 93364320, 109412640, 105846720,
   120187200, 109969920, 132935040, 118329600, 141553440,
   132451440, 148871520, 137229120, 168752640, 148599360,
   171737280, 163900800, 187012800, 169192800, 206025120,
   181466880, 213183360, 200202240, 224259840, 207446400,
   251657280, 219041760, 254864880, 241997760};

int cached_vectors[MAX_NORM + 1][MAX_NM_OF_VECTORS][8];
int cached_idx[MAX_NORM + 1] = {0};

void _cache_vectors(void)
{

  double m = sqrt(2 * MAX_NORM);
  int _s0_end = ceil(2 * m + 1);
  for (int s0 = floor(-2*m); s0 < _s0_end; s0++)
    {
      double ds0 = s0;
      int beg = floor(-m - ds0/2);
      int end = ceil(m - ds0/2 + 1);
      for (int s1 = beg; s1 < end; s1++)
        {
          for (int s2 = beg; s2 < end; s2++)
            {
              for (int s3 = beg; s3 < end; s3++)
                {
                  for (int s4 = beg; s4 < end; s4++)
                    {
                      for (int s5 = beg; s5 < end; s5++)
                        {
                          for (int s6 = beg; s6 < end; s6++)
                            {
                              double _centr = (2 * (s1 + s2 + s3 + s4 + s5 + s6) + ds0)/4;
                              int _end = ceil(m/2 - _centr + 1);
                              for (int s7 = floor(-m/2 - _centr); s7 < _end; s7++)
                                {
                                  int v[8] = {s0, s1, s2, s3, s4, s5, s6, s7};
                                  int _nm = norm(v)/2;
                                  if (_nm < MAX_NORM + 1)
                                    {
                                      int idx = ++cached_idx[_nm];
                                      for (int i = 0; i < 8; i++)
                                        {
                                          cached_vectors[_nm][idx][i] = v[i];
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


void miyawaki_theta(int vs1[MAX_NM_OF_VECTORS][8], int vs2[MAX_NM_OF_VECTORS][8], int vs3[MAX_NM_OF_VECTORS][8],
                    int a, int b, int c, int d, int e, int f)
{
  fmpz_t s0; fmpz_t s1; fmpz_t s2; fmpz_t s3; fmpz_t s4; fmpz_t s5; fmpz_t s6; fmpz_t s7;
  fmpz_t t0; fmpz_t t1; fmpz_t t2; fmpz_t t3; fmpz_t t4; fmpz_t t5; fmpz_t t6; fmpz_t t7;
  fmpz_t u0; fmpz_t u1; fmpz_t u2; fmpz_t u3; fmpz_t u4; fmpz_t u5; fmpz_t u6; fmpz_t u7;
  fmpz_t rl_pt; fmpz_t im_pt; fmpz_t res;
  fmpz_t a0; fmpz_t a1; fmpz_t a2; fmpz_t a3; fmpz_t a4; fmpz_t a5; fmpz_t a6; fmpz_t a7; fmpz_t a8;
  fmpz_t tmp;

  fmpz_init(s0); fmpz_init(s1); fmpz_init(s2); fmpz_init(s3); fmpz_init(s4); fmpz_init(s5); fmpz_init(s6); fmpz_init(s7);
  fmpz_init(t0); fmpz_init(t1); fmpz_init(t2); fmpz_init(t3); fmpz_init(t4); fmpz_init(t5); fmpz_init(t6); fmpz_init(t7);
  fmpz_init(u0); fmpz_init(u1); fmpz_init(u2); fmpz_init(u3); fmpz_init(u4); fmpz_init(u5); fmpz_init(u6); fmpz_init(u7);
  fmpz_init(rl_pt); fmpz_init(im_pt); fmpz_init(res);
  fmpz_init(a0); fmpz_init(a1); fmpz_init(a2); fmpz_init(a3); fmpz_init(a4); fmpz_init(a5); fmpz_init(a6); fmpz_init(a7); fmpz_init(a8);
  fmpz_init(tmp);

  fmpz_zero(res);

  for (int i = 0; i < num_of_vectors[a]; i++)
    {
      for (int j = 0; j < num_of_vectors[b]; j++)
        {
          for (int k = 0; k < num_of_vectors[c]; k++)
            {
              if (inner_prod(vs1[i], vs2[j]) == d)
                {
                  if (inner_prod(vs1[i], vs3[k]) == e)
                    {
                      if (inner_prod(vs2[j], vs3[k]) == f)
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

                          /* Compute rl_pt */
                          fmpz_zero(a0);
                          fmpz_sub(a0, a0, s0);
                          fmpz_zero(a1);
                          fmpz_add(a1, a1, s1);
                          fmpz_zero(a2);
                          fmpz_sub(a2, a2, s3);
                          fmpz_zero(a3);
                          fmpz_add(a3, a3, s4);
                          fmpz_mul(a3, a3, t0);
                          fmpz_addmul(a3, a2, t1);
                          fmpz_addmul(a3, a1, t3);
                          fmpz_addmul(a3, a0, t4);
                          fmpz_zero(a0);
                          fmpz_add(a0, a0, s0);
                          fmpz_zero(a1);
                          fmpz_sub(a1, a1, s2);
                          fmpz_zero(a2);
                          fmpz_add(a2, a2, s3);
                          fmpz_zero(a4);
                          fmpz_sub(a4, a4, s5);
                          fmpz_mul(a4, a4, t0);
                          fmpz_addmul(a4, a2, t2);
                          fmpz_addmul(a4, a1, t3);
                          fmpz_addmul(a4, a0, t5);
                          fmpz_zero(a0);
                          fmpz_sub(a0, a0, s1);
                          fmpz_zero(a1);
                          fmpz_add(a1, a1, s2);
                          fmpz_zero(a2);
                          fmpz_sub(a2, a2, s4);
                          fmpz_zero(a5);
                          fmpz_add(a5, a5, s5);
                          fmpz_mul(a5, a5, t1);
                          fmpz_addmul(a5, a2, t2);
                          fmpz_addmul(a5, a1, t4);
                          fmpz_addmul(a5, a0, t5);
                          fmpz_zero(a0);
                          fmpz_sub(a0, a0, s3);
                          fmpz_zero(a1);
                          fmpz_add(a1, a1, s4);
                          fmpz_zero(a2);
                          fmpz_add(a2, a2, s0);
                          fmpz_zero(a6);
                          fmpz_sub(a6, a6, s1);
                          fmpz_mul(a6, a6, t0);
                          fmpz_addmul(a6, a2, t1);
                          fmpz_addmul(a6, a1, t3);
                          fmpz_addmul(a6, a0, t4);
                          fmpz_zero(a0);
                          fmpz_add(a0, a0, s3);
                          fmpz_zero(a1);
                          fmpz_sub(a1, a1, s5);
                          fmpz_zero(a2);
                          fmpz_sub(a2, a2, s0);
                          fmpz_zero(a7);
                          fmpz_add(a7, a7, s2);
                          fmpz_mul(a7, a7, t0);
                          fmpz_addmul(a7, a2, t2);
                          fmpz_addmul(a7, a1, t3);
                          fmpz_addmul(a7, a0, t5);
                          fmpz_zero(a0);
                          fmpz_sub(a0, a0, s4);
                          fmpz_zero(a1);
                          fmpz_add(a1, a1, s5);
                          fmpz_zero(a2);
                          fmpz_add(a2, a2, s1);
                          fmpz_zero(a8);
                          fmpz_sub(a8, a8, s2);
                          fmpz_mul(a8, a8, t1);
                          fmpz_addmul(a8, a2, t2);
                          fmpz_addmul(a8, a1, t4);
                          fmpz_addmul(a8, a0, t5);
                          fmpz_mul(a8, a8, u0);
                          fmpz_addmul(a8, a7, u1);
                          fmpz_addmul(a8, a6, u2);
                          fmpz_addmul(a8, a5, u3);
                          fmpz_addmul(a8, a4, u4);
                          fmpz_addmul(a8, a3, u5);
                          fmpz_set(rl_pt, a8);

                          /* Compute im_pt */
                          fmpz_zero(a0);
                          fmpz_sub(a0, a0, s3);
                          fmpz_zero(a1);
                          fmpz_add(a1, a1, s4);
                          fmpz_zero(a2);
                          fmpz_add(a2, a2, s0);
                          fmpz_zero(a3);
                          fmpz_sub(a3, a3, s1);
                          fmpz_mul(a3, a3, t0);
                          fmpz_addmul(a3, a2, t1);
                          fmpz_addmul(a3, a1, t3);
                          fmpz_addmul(a3, a0, t4);
                          fmpz_zero(a0);
                          fmpz_add(a0, a0, s3);
                          fmpz_zero(a1);
                          fmpz_sub(a1, a1, s5);
                          fmpz_zero(a2);
                          fmpz_sub(a2, a2, s0);
                          fmpz_zero(a4);
                          fmpz_add(a4, a4, s2);
                          fmpz_mul(a4, a4, t0);
                          fmpz_addmul(a4, a2, t2);
                          fmpz_addmul(a4, a1, t3);
                          fmpz_addmul(a4, a0, t5);
                          fmpz_zero(a0);
                          fmpz_sub(a0, a0, s4);
                          fmpz_zero(a1);
                          fmpz_add(a1, a1, s5);
                          fmpz_zero(a2);
                          fmpz_add(a2, a2, s1);
                          fmpz_zero(a5);
                          fmpz_sub(a5, a5, s2);
                          fmpz_mul(a5, a5, t1);
                          fmpz_addmul(a5, a2, t2);
                          fmpz_addmul(a5, a1, t4);
                          fmpz_addmul(a5, a0, t5);
                          fmpz_zero(a0);
                          fmpz_add(a0, a0, s0);
                          fmpz_zero(a1);
                          fmpz_sub(a1, a1, s1);
                          fmpz_zero(a2);
                          fmpz_add(a2, a2, s3);
                          fmpz_zero(a6);
                          fmpz_sub(a6, a6, s4);
                          fmpz_mul(a6, a6, t0);
                          fmpz_addmul(a6, a2, t1);
                          fmpz_addmul(a6, a1, t3);
                          fmpz_addmul(a6, a0, t4);
                          fmpz_zero(a0);
                          fmpz_sub(a0, a0, s0);
                          fmpz_zero(a1);
                          fmpz_add(a1, a1, s2);
                          fmpz_zero(a2);
                          fmpz_sub(a2, a2, s3);
                          fmpz_zero(a7);
                          fmpz_add(a7, a7, s5);
                          fmpz_mul(a7, a7, t0);
                          fmpz_addmul(a7, a2, t2);
                          fmpz_addmul(a7, a1, t3);
                          fmpz_addmul(a7, a0, t5);
                          fmpz_zero(a0);
                          fmpz_add(a0, a0, s1);
                          fmpz_zero(a1);
                          fmpz_sub(a1, a1, s2);
                          fmpz_zero(a2);
                          fmpz_add(a2, a2, s4);
                          fmpz_zero(a8);
                          fmpz_sub(a8, a8, s5);
                          fmpz_mul(a8, a8, t1);
                          fmpz_addmul(a8, a2, t2);
                          fmpz_addmul(a8, a1, t4);
                          fmpz_addmul(a8, a0, t5);
                          fmpz_mul(a8, a8, u0);
                          fmpz_addmul(a8, a7, u1);
                          fmpz_addmul(a8, a6, u2);
                          fmpz_addmul(a8, a5, u3);
                          fmpz_addmul(a8, a4, u4);
                          fmpz_addmul(a8, a3, u5);
                          fmpz_set(im_pt, a8);

                          /* Computation of tmp = (rl_pt ** 8 - 28 * rl_pt ** 6 * im_pt ** 2 + 70 * rl_pt ** 4 * im_pt ** 4 - 28 * rl_pt ** 2 * im_pt ** 6 + im_pt ** 8)  */
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
  fmpz_clear(s0); fmpz_clear(s1); fmpz_clear(s2); fmpz_clear(s3); fmpz_clear(s4); fmpz_clear(s5); fmpz_clear(s6); fmpz_clear(s7);
  fmpz_clear(t0); fmpz_clear(t1); fmpz_clear(t2); fmpz_clear(t3); fmpz_clear(t4); fmpz_clear(t5); fmpz_clear(t6); fmpz_clear(t7);
  fmpz_clear(u0); fmpz_clear(u1); fmpz_clear(u2); fmpz_clear(u3); fmpz_clear(u4); fmpz_clear(u5); fmpz_clear(u6); fmpz_clear(u7);
  fmpz_clear(rl_pt); fmpz_clear(im_pt); fmpz_clear(res); fmpz_clear(tmp);

}

static void _set_vs(int vs[MAX_NM_OF_VECTORS][8], int a)
{
  int n = num_of_vectors[a];
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < 8; j++)
        {
          vs[i][j] = cached_vectors[a][i][j];
        }
    }
}

int main(void)
{
  _cache_vectors();

  printf("computation for cached vectors done.\n");

  int vs1[MAX_NM_OF_VECTORS][8];

  _set_vs(vs1, num_of_vectors[1]);

  miyawaki_theta(vs1, vs1, vs1, 1, 1, 1, 0, 0, 0);
  return 0;
}

/* Local Variables: */
/* compile-command: "make compile" */
/* End: */
