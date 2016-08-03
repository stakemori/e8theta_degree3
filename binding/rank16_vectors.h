#ifndef _RANK16_VECTORS_H_
#define _RANK16_VECTORS_H_

#define MAX_NORM_RK16 3
/* MAX_NORM_RK16 th Fourier coefficient of the Eisenstein series of weight 8 */
#define MAX_NM_OF_VECTORS_RK16 1050240

extern int num_of_vectors_rk16[8];
short cached_vectors_rk16[MAX_NORM_RK16 + 1][MAX_NM_OF_VECTORS_RK16][16];

void cache_vectors_rk16(void);

#endif /* _RANK16_VECTORS_H_ */
