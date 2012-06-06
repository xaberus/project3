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

/*! \memberof array
 3-point-rule:
 calculates first numerical derrivative of \a f using
 \a h as delta and \a gm \a gp as boundary points */
array_t * array_first_diff_3(double h, array_t * f, double gm, double gp)
{
  int k, N = f->length;
  array_t * df = array_new(N);
  for (k = 0; k < N; k++) {
    // benutze Randpunkte am Rand und nur am Rand ^_^
    double fm, fp;
    if (k == 0) {
      fm = gm;
      fp = f->data[k+1];
    } else if (k == N - 1) {
      fm = f->data[k-1];
      fp = gp;
    } else {
      fm = f->data[k-1];
      fp = f->data[k+1];
    }
    df->data[k] = (fp - fm) / (2 * h);
  }
  return df;
}

/*! \memberof array
 searches for sign changes in f samples on h grid using 3-point-rule
 and assuming a boundary condition of f[-1] = gm, f[n+1] = gp
 node: returns indices cast to double */
array_t * search_der_sign_change_3(double h, array_t * f, double gm, double gp)
{
  array_t * df = array_first_diff_3(h, f, gm, gp);
  array_t * sc = array_new(0);

  double si = df->data[0];
  for (int l = 1; l < df->length - 1; l++) {
    double d = df->data[l];
    //printf("%.17g %.17g\n", s->data[l], d);
    if ((si < 0 && d > 0) || (si > 0 && d < 0)) {
      sc = array_append(sc, l);
    }
    si = d;
  }

  free(df);

  return sc;
}