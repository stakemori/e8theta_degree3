#ifndef _E8VECTORS_H_
#define _E8VECTORS_H_

#include <mpir.h>

#define MAX_NORM 6
#define MAX_NM_OF_VECTORS 60480

extern int * cached_vectors_ptr[];

extern int num_of_vectors[100];

void cache_vectors(void);

#endif /* _E8VECTORS_H_ */
