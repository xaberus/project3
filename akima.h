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

#ifndef _AKIMA_H_
#define _AKIMA_H_

#include "array.h"

/*! struct holding akima interpolation data */
typedef struct akima {
  int      n;  /**< number of points */
  double * s;
  double * a;
  double * b;
  double * c;
  double * d;
} akima_t;

akima_t * akima_new(array_t * x, array_t * f);
double akima_interpolate1(akima_t * a, double x);
array_t * akima_interpolate(akima_t * a, array_t * x);
array_t * akima_dinterpolate(akima_t * a, array_t * x);
void akima_free(akima_t * a);
array_t * akima_zroots(akima_t * a, double s);


#endif /* _AKIMA_H_ */
