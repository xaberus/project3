#ifndef _CARRAY_H_
#define _CARRAY_H_

#include <complex.h>

typedef struct carray {
  int            alloc;
  int            length;
  complex double data[];
  // field elements are allocated after the array
} carray_t;

carray_t * carray_new_sized(int length, int alloc);
carray_t * carray_append(carray_t * z, complex double x);

#endif /* _CARRAY_H_ */