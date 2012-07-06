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
#include <stdio.h>
#include <stdarg.h>
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
array_t * search_der_sign_change_3(double h, array_t * f, double gm, double gp, double tol)
{
  array_t * df = array_first_diff_3(h, f, gm, gp);
  array_t * sc = array_new(0);

  double si = df->data[0];
  for (int l = 1; l < df->length - 1; l++) {
    double d = df->data[l];
    //printf("%.17g %.17g\n", s->data[l], d);
    if ((si < 0 && d > 0) || (si > 0 && d < 0)) {
      if (fabs((f->data[l-1] - f->data[l+1]) / h) > tol) {
        sc = array_append(sc, l);
      }
    }
    si = d;
  }

  free(df);

  return sc;
}

/*! \memberof array
 returns an array with r[i]=fn(z[i]) */
array_t * array_map(array_t * z, double (*fn)(double))
{
  int length = z->length;
  array_t * r = array_new(z->length);
  double * d = r->data, * s = z->data;
  for (int k = 0; k < length; k++) {
    *(d++) = fn(*(s++));
  }
  return r;
}

/*! \memberof array
 returns an array with r[i]=fn(z[i], arg) */
array_t * array_mapv(array_t * z, double (*fn)(double, void*), void * arg)
{
  int length = z->length;
  array_t * r = array_new(z->length);
  double * d = r->data, * s = z->data;
  for (int k = 0; k < length; k++) {
    *(d++) = fn(*(s++), arg);
  }
  return r;
}


/*! \memberof array
 dumps \a argc arrays to file \a name */
void array_dump_to_file(const char name[], const char sep[], int argc, ...)
{
  array_t *argv[argc];
  va_list ap;

  va_start(ap, argc);
  for (int i = 0; i < argc; i++) {
    argv[i] = va_arg(ap, array_t *);
  }
  va_end(ap);

  int length = 0;

  for (int i = 0; i < argc; i++) {
    if (argv[i]) { // NULL means counter
      int l = argv[i]->length;
      length = length < l ? l : length;
    }
  }

  FILE * fp = fopen(name, "w");

  if (fp) {
    for (int n = 0; n < length; n++) {
      for (int i = 0; i < argc; i++) {
        if (argv[i]) {
          double p = argv[i]->data[n];
          if (i > 0) {
            fprintf(fp, "%s%.17e", sep, p);
          } else {
            fprintf(fp, "%.17e", p);
          }
        } else {
          if (i > 0) {
            fprintf(fp, "%s%d", sep, n);
          } else {
            fprintf(fp, "%d", n);
          }
        }
      }
      fputs("\n", fp);
    }
    fclose(fp);
  }
}

/*! \memberof array
 prepare cubic spline interpolation by calculating the a[i]'s */
array_t * array_cspline_prepare(array_t * f, double h)
{
  int k;
  int length = f->length;
  // because the matrix is symmetic and the subdiagonals are never touched by the
  // calculation only two new vectors are needed: diag and b

  // create and initialize matrix diagonal
  array_t * diag = array_new(length - 2);
  for (k = 0; k < diag->length; k++) { diag->data[k] = 4 * h; }

  // create and initialize value vector
  array_t * b = array_new(diag->length);
  for (k = 0; k < b->length; k++) {
    double f1 = f->data[k];
    double f2 = f->data[k + 1];
    double f3 = f->data[k + 2];
    b->data[k] = 6 * ((f3 - f2) / h - (f2 - f1) / h);
  }

  // calculate the new diagonal of the twodiagonal matrix
  for (k = 1 /* ! */; k < b->length; k++) {
    double dd = diag->data[k - 1];
    diag->data[k] += - h * h / dd;
    double tt = b->data[k - 1];
    b->data[k] += - h * tt / dd;
  }

  // this will be the coefficients vector
  array_t * a = array_new(length);

  // natural cspline boundary condition
  a->data[0] = 0;
  a->data[a->length - 1] = 0;

  // calculate the rest of a[i]
  for (k = b->length - 1; k >= 0; k--) {
    a->data[k + 1] = (b->data[k] - h *  a->data[k + 2]) / diag->data[k];
  }

  // dispose of garbage
  free(diag);
  free(b);

  return a;
}

/*! \memberof array
 takes a partition and returns i so, that z is in [v[i],v[i+1]]
 or -1 if no such i exsists */
int array_getmaxindex(array_t * v, double z)
{
  int k;
  if (z >= v->data[0]) {
    for (k = 0; k < v->length; k++) {
      if (v->data[k + 1] >= z) {
        return k;
      }
    }
  }
  return -1;
}

/*! \memberof array
 returns a vector of interpolated values of f for points in x */
array_t * array_cspline_interpolate(array_t * x, array_t * s, array_t * f, array_t * a, double h)
{
  int k;
  array_t * p = array_new(x->length);
  for (k = 0; k < x->length; k++) {
    double z = x->data[k];
    int idx = array_getmaxindex(s, z);
    if (idx >= 0) {
      double dt1 = (z - s->data[idx]);
      double dt2 = (z - s->data[idx + 1]);
      double a1 = a->data[idx];
      double a2 = a->data[idx + 1];
      p->data[k] =
        (dt1 * dt1 * dt1 * a2 - dt2 * dt2 * dt2 * a1) / ( 6 * h)
        + dt1 * (f->data[idx + 1] / h - h * a2 / 6)
        + dt2 * (h * a1 / 6 - f->data[idx] / h);
    } else {
      p->data[k] = 0;
    }
  }

  return p;
}

/*! \memberof array
 returns a interpolated value of f for point \a z */
double array_cspline_ddinterpolate1(double z, array_t * s, array_t * a, double h)
{
  int idx = array_getmaxindex(s, z);
  if (idx >= 0) {
    double dt1 = (z - s->data[idx]);
    double dt2 = (z - s->data[idx + 1]);
    double a1 = a->data[idx];
    double a2 = a->data[idx + 1];
    return (a2 * dt1 - a1 * dt2) / h;
  }
  return 0;
}

/*! \memberof array
 returns a vector of interpolated values of f for points in x */
array_t * array_cspline_dinterpolate(array_t * x, array_t * s, array_t * f, array_t * a, double h)
{
  int k;
  array_t * p = array_new(x->length);
  for (k = 0; k < x->length; k++) {
    double z = x->data[k];
    int idx = array_getmaxindex(s, z);
    if (idx >= 0) {
      double s1 = s->data[idx];
      double s2 = s->data[idx + 1];
      double f1 = f->data[idx];
      double f2 = f->data[idx + 1];
      double a1 = a->data[idx];
      double a2 = a->data[idx + 1];
      double d1 = z - s1;
      double d2 = z - s2;
      p->data[k] = (-(a2*(pow(h,2) - 3*pow(d1,2))) + a1*(pow(h,2) - 3*pow(d2,2)) - 6*f1 + 6*f2)/(6.*h);
    } else {
      p->data[k] = 0;
    }
  }

  return p;
}

/*! \memberof array
 returns cspline extrema with (second derrivative * \a c) < 0 */
array_t * array_cspline_zroots(array_t * s, array_t * f, array_t * a, double h, double c)
{
  printf("array_cspline_zroots(%g) : %d control points\n", c, a->length);

  array_t * p = array_new_sized(0, 100);

  for (int k = 0; k < s->length; k++) {
    double z = s->data[k];
    int idx = k > 0 ? k - 1 : 0;
    double s1 = s->data[idx];
    double s2 = s->data[idx + 1];
    double f1 = f->data[idx];
    double f2 = f->data[idx + 1];
    double a1 = a->data[idx];
    double a2 = a->data[idx + 1];
    double r =
      pow(6*a2*s1 - 6*a1*s2,2)
      + 12*(a1 - a2)*(-6*f1 + 6*f2 - a2*(pow(h,2) - 3*pow(s1,2)) + a1*(pow(h,2) - 3*pow(s2,2)));
    if (r >= 0) {
      double xp = -(sqrt(r) + 6*a2*s1 - 6*a1*s2)/(6.*(a1 - a2));
      double xm = (sqrt(r) - 6*a2*s1 + 6*a1*s2)/(6.*(a1 - a2));
      //printf("-->>-- {%g, (%g, %g), %g} \n", s1, xp, xm, s2);
      if (xm >= s1 && xm <= s2 && array_cspline_ddinterpolate1(xm, s, a, h) * c < 0) {
        printf("--- %g :: %g\n", z, xm);
        p = array_append(p, xm);
      } else if (xp >= s1 && xp <= s2 && array_cspline_ddinterpolate1(xp, s, a, h) * c < 0) {
        printf("+++ %g :: %g\n", z, xp);
        p = array_append(p, xp);
      }
    }
  }

  return p;
}


array_t * array_pcopy(int length, double raw[length])
{
  array_t * r = array_new(length);
  for (int k = 0; k < length; k++) {
    r->data[k] = raw[k];
  }
  return r;
}