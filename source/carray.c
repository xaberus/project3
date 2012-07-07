/******************************************************************************
* Copyright (C) 2012 Pavel Sterin
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
******************************************************************************/

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "carray.h"

#include <stdio.h>

/** \memberof carray
  creates a new array of length \a length */
carray_t * carray_new(int length)
{
  carray_t * r = malloc(sizeof(carray_t) + sizeof(complex double) * length);
  assert(r);
  r->length = length;
  r->alloc = length;
  return r;
}

/** \memberof carray
  creates a new array of length \a length
  while preallocating \a alloc members */
carray_t * carray_new_sized(int length, int alloc)
{
  assert(alloc > length);
  carray_t * r = malloc(sizeof(carray_t) + sizeof(complex double) * alloc);
  assert(r);
  r->length = length;
  r->alloc = alloc;
  return r;
}

/*! \memberof carray
 appends \a x to array \a z reallocating \a z if needed

 the typical usage is:
 \code
     z = carray_append(z, x);
 \endcode
 */
carray_t * carray_append(carray_t * z, complex double x)
{
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

/** \memberof carray
  returns an array of absolute values of \a z using prermuation \a index */
array_t * carray_abs(carray_t * z, int * index)
{
  int length = z->length;
  array_t * r = array_new(length);
  complex double * s = z->data;
  double * d = r->data;
  if (index) {
    for (int k = 0; k < length; k++) {
      *(d++) = cabs(s[index[k]]);
    }
  } else {
    for (int k = 0; k < length; k++) {
      *(d++) = cabs(*(s++));
    }
  }
  return r;
}

carray_t * carray_pcopy(int length, complex double raw[length])
{
  carray_t * r = carray_new(length);
  for (int k = 0; k < length; k++) {
    r->data[k] = raw[k];
  }
  return r;
}
