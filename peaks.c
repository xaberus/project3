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
#include <math.h>

#include "peaks.h"

/*!
 calculates mean value of \a a */
double vec_mean(int length, double * a)
{
  double s = 0;
  for (int k = 0; k < length; k++) {
    s += *(a++) / length;
  }
  return s;
}

/*!
 calculates deviation of \a a from \a ref

 i.e \a ref = mean gives standard deviation */
double vec_deviation(int length, double * a, double ref)
{
  double var = 0;
  for (int k = 0; k < length; k++) {
    // better sum over small numbers at expense of more calculation
    double f = *(a++) - ref;
    var += f * f / length;
  }
  return sqrt(var);
}

/*! this function searches for peaks in dataset \a data using
 the 3-point-rule for picking and standard deviation for selecting
 the peaks */
array_t * peaks_find(array_t * data, int swindow)
{
  array_t * cand = search_der_sign_change_3(1, data, 0, 0);

  int length = data->length;
  int window = swindow * swindow;

  array_t * out = array_new_sized(0, window * 2 + 1);
  array_t * ret = array_new(0);

  for (int k = 0; k < cand->length; k++) {
    int pos = cand->data[k], is, ie;
    if (pos - window < 0) { is = 0; } else { is = pos - window; }
    if (pos + window >= length) {ie = length - 1; } else { ie = pos + window; }

    (void) ie; (void) is;

    int span = ie - is + 1;

    double m = vec_mean(span, data->data + is);
    double s = vec_deviation(span, data->data + is, m);
    double h = 1;

    if (pos - swindow < 0) { is = 0; } else { is = pos - swindow; }
    if (pos + swindow >= length) {ie = length - 1; } else { ie = pos + swindow; }

    out->length = 0;

    for (int i = is; i <= ie; i++) {
      double ak = data->data[i];
      if (ak > 0 && (ak - m) > (h * s)) {
        out = array_append(out, i);
      }
    }

    for (int k = 0, l = 1; l < out->length; l++) {
      int i = out->data[k];
      int j = out->data[l];
      if (abs(j - i) <= window) {
        double xi = data->data[i];
        double xj = data->data[j];
        if (xi < xj) {
          out->data[k] = -1;
          k = l;
        } else {
          out->data[l] = -1;
        }
      } else {
        k = l;
      }
    }

    for (int k = 0; k < out->length; k++) {
      if (out->data[k] >= 0) {
        ret = array_append(ret, out->data[k]);
      }
    }
  }

  free(out);
  free(cand);

  return ret;
}
