#include <math.h>
#include "e8vectors.h"
#include <stdlib.h>
#include <memory.h>

static inline int norm(int s[8])
{
  return (s[0]*s[0] + s[1]*s[1] + s[2]*s[2] + s[3]*s[3] + s[4]*s[4] + s[5]*s[5] + s[6]*s[6] +
          s[7]*s[7]) >> 2;
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

/* Convert (s0, ..., s7) to eulidian vector */
static void _convert_to_euclid_vector_e8(int vec[8])
{
  vec[7] = 4 * vec[7];
  for (int i = 1; i < 7; i++)
    {
      vec[7] += 2 * vec[i];
    }
  vec[7] += vec[0];
  for (int i = 1; i < 7; i++)
    {
      vec[i] = 2 * vec[i] + vec[0];
    }
}

static int cached_vectors0[1][8] = {0};
static int cached_vectors1[240][8] = {0};
static int cached_vectors2[2160][8] = {0};
static int cached_vectors3[6720][8] = {0};
static int cached_vectors4[17520][8] = {0};
static int cached_vectors5[30240][8] = {0};
static int cached_vectors6[60480][8] = {0};

int * cached_vectors_ptr[] = {cached_vectors0[0],
                              cached_vectors1[0],
                              cached_vectors2[0],
                              cached_vectors3[0],
                              cached_vectors4[0],
                              cached_vectors5[0],
                              cached_vectors6[0]};
static int cached_idx[MAX_NORM + 1] = {0};

static void _cache_vectors(void)
{

  double m = sqrt(2 * MAX_NORM);
  int _s0_end = ceil(2 * m + 1);
  int _tmp;
  for (int s0 = floor(-2*m); s0 < _s0_end; s0++)
    {
      int nrm0 = s0 * s0;
      double ds0 = s0;
      int beg = floor(-m - ds0/2);
      int end = ceil(m - ds0/2 + 1);
      for (int s1 = beg; s1 < end; s1++)
        {
          _tmp = s0 + 2 * s1;
          int nrm1 = nrm0 + _tmp * _tmp;
          if (nrm1 < 8 * MAX_NORM + 1)
            {
              for (int s2 = beg; s2 < end; s2++)
                {
                  _tmp = s0 + 2 * s2;
                  int nrm2 = nrm1 + _tmp * _tmp;
                  if (nrm2 < 8 * MAX_NORM + 1)
                    {
                      for (int s3 = beg; s3 < end; s3++)
                        {
                          _tmp = s0 + 2 * s3;
                          int nrm3 = nrm2 + _tmp * _tmp;
                          if (nrm3 < 8 * MAX_NORM + 1)
                            {
                              for (int s4 = beg; s4 < end; s4++)
                                {
                                  _tmp = s0 + 2 * s4;
                                  int nrm4 = nrm3 + _tmp * _tmp;
                                  if (nrm4 < 8 * MAX_NORM + 1)
                                    {
                                      for (int s5 = beg; s5 < end; s5++)
                                        {
                                          _tmp = s0 + 2 * s5;
                                          int nrm5 = nrm4 + _tmp * _tmp;
                                          if (nrm5 < 8 * MAX_NORM + 1)
                                            {
                                              for (int s6 = beg; s6 < end; s6++)
                                                {
                                                  _tmp = s0 + 2 * s6;
                                                  int nrm6 = nrm5 + _tmp * _tmp;
                                                  if (nrm6 < 8 * MAX_NORM + 1)
                                                    {
                                                      double _centr = (2 * (s1 + s2 + s3 + s4 + s5 + s6) + ds0)/4;
                                                      int _end = ceil(m/2 - _centr + 1);
                                                      for (int s7 = floor(-m/2 - _centr); s7 < _end; s7++)
                                                        {
                                                          _tmp = (s0 + 2 * (s1 + s2 + s3 + s4 + s5 + s6) + 4 * s7);
                                                          int nrm7 = nrm6 + _tmp * _tmp;
                                                          if (nrm7 < 8 * MAX_NORM + 1)
                                                            {
                                                              int v[8] = {s0, s1, s2, s3, s4, s5, s6, s7};
                                                              _convert_to_euclid_vector_e8(v);
                                                              int _nm = nrm7 >> 3;
                                                              int idx = cached_idx[_nm]++;
                                                              memcpy(cached_vectors_ptr[_nm] + 8 * idx, v, 8 * sizeof(int));
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
                }
            }
        }
    }
}

void cache_vectors(void)
{
  int b = 1;
  for (int i = 0; i < 8; i++)
    {
      b = (! (cached_vectors_ptr[1] + 8 * 0)[i]) && b;
    }
  if (b)
    {
      _cache_vectors();
    }
}

/* Local Variables: */
/* compile-command: "cd ..; make compile-theta_vectors" */
/* End: */
