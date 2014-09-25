#ifndef V_DEFS
#define V_DEFS

#define SQR(x)			((x) * (x))
#define CUBE(x)			((x) * (x) * (x))
#define MIN(x1, x2)		(((x1) < (x2)) ? (x1) : (x2))
#define MAX(x1, x2)		(((x1) > (x2)) ? (x1) : (x2))

typedef struct {double x, y, z;} VecR;
typedef struct {int x, y, z;} VecI;

#define MAKE_FILENAME(full, file)						\
	{strcpy (full, foldername);							\
	strcat (full,file);}

/* memory allocation */
#define ALLOC_MEM(a, n, t) a = (t *) malloc ((n) * sizeof (t))

#define ALLOC_MEM2(a, n1, n2, t)						\
	ALLOC_MEM (a, n1, t *);								\
	ALLOC_MEM (a[0], (n1) * (n2), t);					\
	for (k__0 = 1; k__0 < (n1); k__0 ++) a[k__0] = a[k__0 - 1] + (n2);

/*#define ALLOC_MEM3(a, n1, n2, n3, t) 					\
		t (* a)[n2][n3] = (t (*)[n2][n3]) malloc ((n1) * (n2) * (n3) * sizeof(t));
*/
#define ALLOC_MEM3(var, x, y, z, type)							\
		do { 													\
			(var) = malloc(sizeof(type **) * (x));				\
			if ((var) == NULL) {								\
				perror("Error allocating memory.\n");			\
			}													\
			int i, j;											\
			for (i = 0; i < x; i++) {							\
				(var)[i] = malloc(sizeof(type*) * (y));			\
				if ((var)[i] == NULL) {							\
					perror("Error allocating memory.\n");		\
				}												\
				for (j = 0; j < y; j++) {						\
					(var)[i][j] = malloc(sizeof(type) * (z));	\
					if ((var)[i][j] == NULL) {					\
						perror("Error allocating memory.\n");	\
					}											\
				}												\
			}													\
		} while(0)

/*#define ALLOC_MEM3(a, n1, n2, n3, t)									\
		ALLOC_MEM2 (a, n1, n2, t *);									\
		ALLOC_MEM2 (a[0][0], (n1), (n2) * (n3), t *);						\
		for (k__1 = 1; k__1 < n1; k__1 ++) {									\
			for (k__2 = 1; k__2 < n2; k__2 ++) a[k__1][k__2] = a[k__1-1][k__2-1] + n3;}
*/
#include "in_proto.h"

#define SET_TO_ZERO(array_element, i, j, k)						\
	array_element[i][j][k] = 0.;


#define V_SET(v, sx, sy, sz)							\
	(v).x = sx;											\
	(v).y = sy;											\
	(v).z = sz;

#define V_SET_ALL(v, s)									\
	V_SET (v, s, s, s);

#define V_ZERO(v)										\
	V_SET_ALL (v, 0);

#endif
