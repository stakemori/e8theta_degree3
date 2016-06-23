#ifndef _E8VECTORS_H_
#define _E8VECTORS_H_

#define MAX_NORM 5
#define MAX_NM_OF_VECTORS 30240

extern int num_of_vectors[100];

void _set_vs3(int vs1[MAX_NM_OF_VECTORS][8],
              int vs2[MAX_NM_OF_VECTORS][8],
              int vs3[MAX_NM_OF_VECTORS][8],
              int a, int b, int c);

void _cache_vectors(void);

#endif /* _E8VECTORS_H_ */
