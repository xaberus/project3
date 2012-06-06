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

inline static
void peek_data(int l, double d[l], int i, int k, int w, double * N1, double * N2)
{
  int n = i - k;
  int m = i + k;

  int is, ie;

  if (n < 0) { is = 0; } else { is = n; }
  if (m >= l) { ie = l - 1; } else { ie = m; }

  for (int p = n, s = 0, ss = 0; p <= m + w; p++) {
    if (p < is) {
      if (p != i) {
        N1[s++] = d[is]/fabs(p-i);
      }
      N2[ss++] = d[is]/fabs(p-i);
    }
    else if (p > ie) {
      if (p != i) {
        N1[s++] = d[ie]/fabs(p-i);
      }
      N2[ss++] = d[ie]/fabs(p-i);
    } else {
      if (p != i) {
        N1[s++] = d[p];
      }
      N2[ss++] = d[p];
    }
  }
}

double kernel(double x)
{
  return exp(-x*x/2)/sqrt(2 * M_PI);
}

double entropy_cmpare(int k, int w, double N1[2*k+w+1], double N2[2*k+w+2])
{
  int m1 = 2 * k;
  double h1 = 0;
  for (int i = 0; i < m1; i++) {
    double q = fabs(N1[i]-N1[i+w]), s = 0;
    for (int j = 0; j < m1; j++) {
      s += kernel((N1[i]-N1[j])/q);
    }
    double p = 1.0/(m1*q) * s;
    h1 += - p * log(p);
  }
  int m2 = 2 * k + 1;
  double h2 = 0;
  for (int i = 0; i < m2; i++) {
    double q = fabs(N2[i]-N2[i+w]), s = 0;
    for (int j = 0; j < m2; j++) {
      s += kernel((N2[i]-N2[j])/q);
    }
    double p = 1.0/((m2)*q) * s;
    h1 += - p * log(p);
  }
  return h1 - h2;
}

double max_cmpare(array_t * d, int i, int k)
{
  int n = i - k;
  int m = i + k;

  int is, ie;

  if (n < 0) { is = 0; } else { is = n; }
  if (m >= d->length) { ie = d->length - 1; } else { ie = m; }

  double xi = d->data[i];
  double m1 = -1.0/0.0, m2 = -1.0/0.0;
  for (int p = is; p < i; p++) {
    double r = xi - d->data[p];
    if (r > m1) {
      m1 = r;
    }
  }
  for (int p = i + 1; p <= ie; p++) {
    double r = xi - d->data[p];
    if (r > m2) {
      m2 = r;
    }
  }

  return (m1 + m2)/2;
}

void peak_erase(array_t * a, int i, int k)
{
  int n = i - k;
  int m = i + k;

  int is, ie;

  if (n < 0) { is = 0; } else { is = n; }
  if (m >= a->length) { ie = a->length - 1; } else { ie = m; }

  for (int p = is; p < ie; p++) {
    double x = fabs(p - i);
    a->data[p] *= exp(-k/x);
  }
}

int cmp_double(const void * a, const void * b)
{
  double A = *((double *)a);
  double B = *((double *)b);
  if (A < B) {
    return -1;
  } else if (A > B) {
    return 1;
  }
  return 0;
}

/*! this function searches for peaks in dataset \a data using
 the 3-point-rule for picking and standard deviation for selecting
 the peaks */
array_t * peaks_find(array_t * data, int swindow)
{
  double h = 1.9;

  int length = data->length;

  array_t * out = array_new(0);
  array_t * ret = array_new(0);
  array_t * a = array_new(length);

  for (int k = 0; k < length; k++) {
    /*double N1[swindow * 2 + w + 1], N2[swindow * 2 + w + 2];
    peek_data(length, data->data, k, swindow, w, N1, N2);
    double r = entropy_cmpare(swindow, w, N1, N2);*/
    double r = max_cmpare(data, k, swindow);
    a->data[k] = r > 0 ? r : 0;
  }

  int w = swindow;
  do {
    w *= 1.5;
    h *= h;
    out->length = 0;
    double mean = vec_mean(length, a->data);
    double sdev = vec_deviation(length, a->data, mean);

    for (int k = 0; k < length; k++) {
      double ak = a->data[k];
      if (ak > 0 && (ak - mean) > (h * sdev)) {
        out = array_append(out, k);
      }
    }

    for (int k = 0, l = 1; l < out->length; l++) {
      int i = out->data[k];
      int j = out->data[l];
      if (abs(j - i) <= swindow) {
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
        peak_erase(a, out->data[k], w);
      }
    }
  } while(out->length);

  free(a);
  free(out);

  qsort(ret->data, sizeof(double), ret->length, cmp_double);

  return ret;
}
