#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "akima.h"

akima_t * akima_new(array_t * x, array_t * f)
{
  assert(x->length == f->length);
  akima_t * a = malloc(sizeof(akima_t)); assert(a);
  a->n = x->length - 1;
  assert(a->n >= 2);
  a->s = malloc(sizeof(double) * (a->n + 1)); assert(a->s);
  a->a = malloc(sizeof(double) * (a->n + 1)); assert(a->a);
  a->c = malloc(sizeof(double) * a->n); assert(a->c);
  a->d = malloc(sizeof(double) * a->n); assert(a->d);

  double * h = malloc(sizeof(double) * a->n); assert(h);
  double * m = malloc(sizeof(double) * (a->n + 4)); assert(h);
  double * tl = malloc(sizeof(double) * (a->n + 1)); assert(h);
  double * tr = malloc(sizeof(double) * (a->n + 1)); assert(h);

  int n = a->n;

  for (int i = 0; i <= n; i++) {
    a->s[i] = x->data[i];
    a->a[i] = f->data[i];
  }

  for (int i = 0; i < n; i++) {
    h[i] = a->s[i+1] - a->s[i];
    m[2 + i] = (f->data[i+1] - f->data[i]) / h[i];
  }

  m[2 - 2] = 3 * m[2 + 0] - 2 * m[2 + 1];
  m[2 - 1] = 2 * m[2 + 0] - m[2 + 1];
  m[2 + n] = 2 * m[2 + n - 1] - m[2 + n - 2];
  m[2 + n + 1] = 3 * m[2 + n - 1] - 2 * m[2 + n - 2];

  for (int i = 0; i <= n; i++) {
    double l = fabs(m[2 + i - 2] - m[2 + i - 1]);
    double ne = l + fabs(m[2 + i] - m[2 + i + 1]);
    if (ne > 0) {
      double alpha = l / ne;
      tr[i] = tl[i] = (1 - alpha) * m[2 + i - 1] + alpha * m[2 + i];
    } else {
      tl[i] = m[2 + i - 1];
      tr[i] = m[2 + i];
    }
  }

  a->b = tr;

  for (int i = 0; i < n; i++) {
    a->c[i] = (3 * m[2 + i] - 2 * tr[i] - tl[i + 1]) / h[i];
    a->d[i] = (tr[i] + tl[i + 1] - 2 * m[2 + i]) / h[i] / h[i];
  }

  free(tl);
  free(h);
  free(m);

  return a;
}

int getmaxindex(int length, double v[length], double z)
{
  if (length > 0) {
    int k;
    if (z >= v[0]) {
      for (k = 0; k < length; k++) {
        if (v[k + 1] >= z) {
          return k;
        }
      }
    }
  }
  return -1;
}

double akima_interpolate1(akima_t * a, double x)
{
  int i = getmaxindex(a->n + 1, a->s, x);
  double d = (x - a->s[i]);
  if (i >= 0) {
    return
      a->a[i]
        + a->b[i] * d
        + a->c[i] * d * d
        + a->d[i] * d * d * d;
  }
  return 0;
}

double akima_dinterpolate1(akima_t * a, double x)
{
  int i = getmaxindex(a->n + 1, a->s, x);
  double d = (x - a->s[i]);
  if (i >= 0) {
    return a->b[i] + 2 * a->c[i] * d + 3 * a->d[i] * d * d;
  }
  return 0;
}

double akima_ddinterpolate1(akima_t * a, double x)
{
  int i = getmaxindex(a->n + 1, a->s, x);
  double d = (x - a->s[i]);
  if (i >= 0) {
    return 2 * a->c[i] + 6 * a->d[i] * d;
  }
  return 0;
}

array_t * akima_dinterpolate(akima_t * a, array_t * x)
{
  array_t * r = array_new(x->length);
  for (int k = 0; k < x->length; k++) {
    r->data[k] = akima_dinterpolate1(a, x->data[k]);
  }
  return r;
}

array_t * akima_interpolate(akima_t * a, array_t * x)
{
  array_t * r = array_new(x->length);
  for (int k = 0; k < x->length; k++) {
    r->data[k] = akima_interpolate1(a, x->data[k]);
  }
  return r;
}

void akima_free(akima_t * a)
{
  free(a->a);
  free(a->b);
  free(a->c);
  free(a->d);
  free(a->s);
  free(a);
}

array_t * akima_zroots(akima_t * a, double s)
{
  printf("akima_zroots(%g) : %d control points\n", s, a->n+1);

  array_t * p = array_new_sized(0, 100);

  for (int k = 0; k < a->n; k++) {
    double s1 = a->s[k];
    double s2 = a->s[k + 1];
    double b = a->b[k];
    double c = a->c[k];
    double d = a->d[k];
    double r = c * c - 3 * b * d;
    if (r >= 0) {
      double xm = (-c - sqrt(r)) / (3 * d) + s1;
      double xp = (-c + sqrt(r)) / (3 * d) + s1;

      //printf("-->>-- {%g, (%g, %g), %g} \n", s1, xp, xm, s2);
      if (xm >= s1 && xm <= s2 && akima_ddinterpolate1(a, xm) * s < 0) {
        //printf("--- %g :: %g\n", s1, xm);
        p = array_append(p, xm);
      } else if (xp >= s1 && xp <= s2 && akima_ddinterpolate1(a, xp) * s < 0) {
        //printf("+++ %g :: %g\n", s1, xp);
        p = array_append(p, xp);
      }
    }
  }

  return p;
}
