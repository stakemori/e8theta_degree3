#ifndef _E8VECTORS_H_
#define _E8VECTORS_H_

#include <mpir.h>

#define MAX_NORM 6
#define MAX_NM_OF_VECTORS 60480

extern int num_of_vectors[100];

void _set_vs3(int vs1[MAX_NM_OF_VECTORS][8],
              int vs2[MAX_NM_OF_VECTORS][8],
              int vs3[MAX_NM_OF_VECTORS][8],
              int a, int b, int c);

void cache_vectors(void);

void _set_vs(int vs[MAX_NM_OF_VECTORS][8], int a);

#endif /* _E8VECTORS_H_ */
