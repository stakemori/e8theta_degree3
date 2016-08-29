#ifndef _VECTOR_UTILS_H_
#define _VECTOR_UTILS_H_
#include <stdlib.h>


void sort_int_vec(int * base, size_t len);
void _normalize_vec_w_sign(int * vec, size_t len);
void normalize_vec_w_indices(int * vec,
                             int w_sign_indices[16],
                             int wo_sign_indices_array[8][16]);

void normalize_vec_last_len_elts(int * vec, int vec_len, int len);


#endif /* _VECTOR_UTILS_H_ */
