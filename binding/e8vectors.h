#ifndef _E8VECTORS_H_
#define _E8VECTORS_H_

#include <fmpz.h>

#define MAX_NORM 6
#define MAX_NM_OF_VECTORS 60480

extern int num_of_vectors[100];

void _set_vs3(int vs1[MAX_NM_OF_VECTORS][8],
              int vs2[MAX_NM_OF_VECTORS][8],
              int vs3[MAX_NM_OF_VECTORS][8],
              int a, int b, int c);

void _cache_vectors(void);

char * _store_fmpz_str_using_malloc(fmpz_t x);

#endif /* _E8VECTORS_H_ */
