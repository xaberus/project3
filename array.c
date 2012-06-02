#include <stdlib.h>
#include <assert.h>

#include "array.h"

/** \memberof array
  creates a new array of length \a length
  while preallocating \a alloc members */
array_t * array_new_sized(int length, int alloc)
{
  assert(alloc > length);
  array_t * r = malloc(sizeof(array_t) + sizeof(double) * alloc); assert(r);
  if (!r) { abort(); }
  r->length = length;
  r->alloc = alloc;
  return r;
}

/** \memberof array
  creates a new array of length \a length */
array_t * array_new(int length)
{
  array_t * r = malloc(sizeof(array_t) + sizeof(double) * length); assert(r);
  if (!r) { abort(); }
  r->length = length;
  r->alloc = length;
  return r;
}

/*! \memberof array
 creates a array with an equipatition of the
 intervall [start;end] with \a length-1 subintevalls */
array_t * array_equipart(double start, double end, int length)
{
  int k;
  double range = end - start;
  array_t * r = array_new(length);

  for (k = 0; k < length; k++) {
    r->data[k] = start + range / (length - 1) * k;
  }

  return r;
}

/*! \memberof array
 appends \a x to array \a z reallocating \a z if needed

 the typical usage is:
 \code
     z = array_append(z, x);
 \endcode
 */
array_t * array_append(array_t * z, double x)
{
  if (z->length < z->alloc) {
    z->data[z->length++] = x;
    return z;
  } else {
    int alloc = (z->length+1) * 2;
    array_t * t = realloc(z, sizeof(array_t) + sizeof(double) * alloc); assert(t);
    t->data[t->length++] = x;
    t->alloc = alloc;
    return t;
  }
}
