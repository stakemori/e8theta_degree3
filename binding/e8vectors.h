#ifndef _E8VECTORS_H_
#define _E8VECTORS_H_

#include <mpir.h>

#define MAX_NORM 7
#define MAX_NM_OF_VECTORS 82560
#define MAX_NM_REPRS 82560

extern int * cached_vectors_ptr[];

extern int num_of_vectors[100];

void cache_vectors(void);
int repr_modulo_autom(int n, int reprs[MAX_NM_REPRS][8],
                      unsigned int num_of_classes[MAX_NM_REPRS]);
int repr_modulo_autom_w_indices(int vecs[MAX_NM_OF_VECTORS][8],
                                int num_of_vecs,
                                int reprs[MAX_NM_REPRS][8],
                                unsigned int num_of_classes[MAX_NM_REPRS],
                                int w_sign_indices[8],
                                int wo_sign_indices_array[8][16]);

#endif /* _E8VECTORS_H_ */
