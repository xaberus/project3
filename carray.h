#ifndef _CARRAY_H_
#define _CARRAY_H_

#include <complex.h>

/*! A rezizeable vector of complex doubles */
typedef struct carray {
  int            alloc;    /**< allocated length */
  int            length;   /**< used length */
  /*! c99 flexible array tail with \a alloc members
      (field elements are allocated after the header) */
  complex double data[];
} carray_t;

carray_t * carray_new(int length);
carray_t * carray_new_sized(int length, int alloc);
carray_t * carray_append(carray_t * z, complex double x);

#endif /* _CARRAY_H_ */