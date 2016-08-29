#ifndef _RANK16_VECTORS_H_
#define _RANK16_VECTORS_H_

#define MAX_NORM_RK16 4
/* MAX_NORM_RK16 th Fourier coefficient of the Eisenstein series of weight 8 */
#define MAX_NM_OF_VECTORS_RK16 7926240

#define MAX_NM_REPRS_RK16 116240

extern int num_of_vectors_rk16[8];
extern int * cached_vectors_rk16_ptr[];

void cache_vectors_rk16(void);
int repr_modulo_autom_rk16(int n, int reprs[MAX_NM_REPRS_RK16][16], unsigned int num_of_classes[MAX_NM_REPRS_RK16]);
int repr_modulo_autom_rk16_w_indices(int vecs[MAX_NM_OF_VECTORS_RK16][16],
                                     int num_of_vecs,
                                     int reprs[MAX_NM_REPRS_RK16][16],
                                     unsigned int num_of_classes[MAX_NM_REPRS_RK16],
                                     int w_sign_indices[16],
                                     int wo_sign_indices_array[8][16]);
void set_w_sign_indices_rk16(int indices[16], const int vec[16]);
void set_wo_sign_indices_array(int indices_array[8][16], const int vec[16]);
void set_w_sign_indices_rk16_2(int indices[16], const int vec1[16],
                               const int vec2[16]);
void set_wo_sign_indices_array2(int indices_array[8][16], const int vec1[16],
                                const int vec2[16]);

#endif /* _RANK16_VECTORS_H_ */
