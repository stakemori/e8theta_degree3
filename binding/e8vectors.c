#include <math.h>
#include "e8vectors.h"

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
                                      int idx = cached_idx[_nm]++;
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


void _set_vs(int vs[MAX_NM_OF_VECTORS][8], int a)
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

void _set_vs3(int vs1[MAX_NM_OF_VECTORS][8],
              int vs2[MAX_NM_OF_VECTORS][8],
              int vs3[MAX_NM_OF_VECTORS][8],
              int a, int b, int c)
{
  _set_vs(vs1, a);
  _set_vs(vs2, b);
  _set_vs(vs3, c);
}


/* Local Variables: */
/* compile-command: "gcc e8vectors.c -I. -O2 -shared -fPIC -o libe8vectors.so" */
/* End: */
