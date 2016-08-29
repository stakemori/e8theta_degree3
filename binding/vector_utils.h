#ifndef _VECTOR_UTILS_H_
#define _VECTOR_UTILS_H_
#include <stdlib.h>


void sort_int_vec(int * base, size_t len);
void _normalize_vec_w_sign(int * vec, size_t len);
void normalize_vec_w_indices(int * vec,
                             int w_sign_indices[16],
                             int wo_sign_indices_array[8][16]);

void normalize_vec_last_len_elts(int * vec, int vec_len, int len);

void set_w_sign_indices(int * indices, const int * vec, int vec_len, int zero_len);

void set_w_sign_indices_2(int * indices, const int * vec1, const int * vec2,
                          int vec_len, int zero_len);

void set_wo_sign_indices_array(int indices_array[8][16], const int * vec, int vec_len,
                               int zero_len);


void set_wo_sign_indices_array2(int indices_array[8][16], const int * vec1,
                                const int * vec2, int vec_len, int zero_len);
#endif /* _VECTOR_UTILS_H_ */
