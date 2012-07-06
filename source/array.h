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

#ifndef _ARRAY_H_
#define _ARRAY_H_

/*! A rezizeable vector of doubles */
typedef struct array {
  int    alloc;   /**< allocated length */
  int    length;  /**< used length */
  /*! c99 flexible array tail with \a alloc members
      (field elements are allocated after the header) */
  double data[];
} array_t;

array_t * array_new(int length);
array_t * array_new_sized(int length, int alloc);
array_t * array_equipart(double start, double end, int points);
array_t * array_append(array_t * z, double x);
array_t * array_map(array_t * z, double (*fn)(double));
array_t * array_mapv(array_t * z, double (*fn)(double, void*), void * arg);
void array_dump_to_file(const char name[], const char sep[], int argc, ...);

array_t * array_first_diff_3(double h, array_t * f, double gm, double gp);
array_t * search_der_sign_change_3(double h, array_t * f, double gm, double gp, double tol);

array_t * array_cspline_prepare(array_t * f, double h);
int array_getmaxindex(array_t * v, double z);
array_t * array_cspline_interpolate(array_t * x, array_t * s, array_t * f, array_t * a, double h);
array_t * array_cspline_zroots(array_t * s, array_t * f, array_t * a, double h, double c);
array_t * array_cspline_dinterpolate(array_t * x, array_t * s, array_t * f, array_t * a, double h);

array_t * array_pcopy(int length, double raw[length]);

#endif /* _ARRAY_H_ */
