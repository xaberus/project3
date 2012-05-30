#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "carray.h"

carray_t * carray_new_sized(int length, int alloc) {
  assert(alloc > length);
  carray_t * r = malloc(sizeof(carray_t) + sizeof(complex double) * alloc);
  if (!r) { abort(); }
  r->length = length;
  r->alloc = length;
  return r;
}

carray_t * carray_append(carray_t * z, complex double x) {
  if (z->length < z->alloc) {
    z->data[z->length++] = x;
    return z;
  } else {
    int alloc = (z->length+1) * 2;
    carray_t * t = realloc(z, sizeof(carray_t) + sizeof(complex double) * alloc); assert(t);
    t->data[t->length++] = x;
    t->alloc = alloc;
    return t;
  }
}