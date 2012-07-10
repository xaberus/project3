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

#include "simulation.h"
#include "numerov.h"

int maximal_iterations = 5000;

typedef struct {
  array_t * V;
  double h;
  double psi1;
  double psi2;
} numerov_t;

inline static
long double numerov_step(double T1, double T2, double T3, long double f1, long double f2)
{
  return (f1 * (T1 - 1.) + 2. * f2 * (1. + 5. * T2)) / (1. - T3);
}

long double numerov_integrate(long double E, void * arg)
{
  numerov_t * num = arg;
  int length = num->V->length;
  double * V = num->V->data;
  long double psi1 = num->psi1, psi2 = num->psi2;
  double hh = num->h * num->h;
  for (int k = 2; k < length; k++) {
    double T1 = - hh * (E - V[k - 2]) / 12.;
    double T2 = - hh * (E - V[k - 1]) / 12.;
    double T3 = - hh * (E - V[k]) / 12.;
    double t = numerov_step(T1, T2, T3, psi1, psi2);
    psi1 = psi2;
    psi2 = t;
  }
  return psi2;
}

/* newton: secant method
 * seek for zero in the intervall, with wt most one change of sign of first derrivative of fn */
double search_zero(long double (*fn)(long double, void*), void * arg, double min, double step, double max, double fmax) {
  //printf("searching in: [%.17g:%.17g] : %.17g ~ %.17g\n", min, max, step, fmax);
  long double xm, x, xp, a, b, f, fm, df;
  long double pp = (max - min) / step;
  int k;
  xm = min, x = xm + (max - min) / 2., xp, a, b;

  for (k = 0; x >= min && x <= max && k < maximal_iterations; k++) {
    if (xm == x) {
      return x;
    }
    f = fn(x, arg)/fmax;
    fm = fn(xm, arg)/fmax;
    df = (f-fm);
    if (df == 0.0) { // worst case scenario....
      goto bisect;
    }

    //printf("->sz: {%.17g, %.17g, %.17g, %.17g}\n", f, fm, xm, x);
    xp = x - (x - xm)/df*f;
    //printf("sz: [%.17g;%.17g]:%.17g {%.17g, %.17g, %.17g}\n", min, max, step, xm, x, xp);
    //printf("sz: [%g,%g]:%g {%g, %g, %g}\n", min, max, step, xm, x, xp);
    if (xp > max || xp < min) {
      goto bisect;
    }
    if (x == xp) {
      return x;
    }
    xm = x;
    x = xp;
    continue;
bisect:
      // wenn wir hier landen, dann war x zu nah an einem Extrempunkt, versuche einen besseren
      // Startwert mit Bisektionsverfahren zu finden
      a = fn(min, arg) / fmax;
      b = fn(max, arg) / fmax;
      if (a * b < 0) {
        long double t = (min + max) / 2.;
        long double w = fn(t, arg) / fmax;
        if (w == 0) {
          return t;
        }
        if (a * w < 0) {
          max = t;
        } else {
          min = t;
        }
        xm = min;
        x = xm + (max - min) / pp;
        continue;
      }
      break;
  }
  printf("found no zero in: [%.17g:%.17g] : %.17g ~ %.17g\n", min, max, step, fmax);
  return 0.0/0.0; // nan - nichts gefunden
}

double getmaxabs(array_t * fn, int start, int end)
{
  double max = 0.0;
  double * f = fn->data;
  for (int k = start; k < end; k++) {
    double z = fabs(*(f++));
    if (z > max) {
      max = z;
    }
  }
  return max;
}

/*! \memberof array
 returns an array with r[i]=fn(z[i], arg) */
array_t * array_mapv(array_t * z, long double (*fn)(long double, void*), void * arg)
{
  int length = z->length;
  array_t * r = array_new(z->length);
  double * d = r->data, * s = z->data;
  for (int k = 0; k < length; k++) {
    *(d++) = fn(*(s++), arg);
  }
  return r;
}

array_t * numerov_energies(double dx, array_t * V, double min, double max)
{

  array_t * zp = array_new_sized(0, 100);

  int res = 4*4096;

  numerov_t num = {V, dx, 0, 1e-14};

  //double E = -1.44966298553678342e+02;
  //double E = -1.38753187223842588e+02;
  //double E = -2.01916470909939498e+01;
  /*double E = 102.5;//+10e-6;
  {
    int length = num.V->length;
    array_t * psia = array_new(length);
    double * psi = psia->data;
    double * V = num.V->data;
    psi[0] = num.psi1;
    psi[1] = num.psi2;
    double hh = num.h * num.h;
    for (int k = 2; k < length; k++) {
      double T1 = - hh * (E - V[k - 2]) / 12.;
      double T2 = - hh * (E - V[k - 1]) / 12.;
      double T3 = - hh * (E - V[k]) / 12.;
      psi[k] = numerov_step(T1, T2, T3, psi[k-2], psi[k-1]);
    }
    array_dump_to_file("scarep", " ", 1, psia);
    free(psia);
  }*/

  array_t * Epos = array_equipart(min, max, res);
  array_t * fn = array_mapv(Epos, numerov_integrate, &num);

  //array_dump_to_file("score", " ", 2, Epos, fn);

  double smin, smax, z = 0;
  array_t * sc = search_der_sign_change_3(Epos->data[1] - Epos->data[0], fn, 0, 0, 0);

  /*FILE * fp = fopen("sumo", "w");
  for (int k = 0; k < sc->length - 1; k++) {
    fprintf(fp, "%g\n", Epos->data[(int)sc->data[k]]);
  }
  fclose(fp);*/

  double eres = 0.00000001;

  if (sc->length > 0) {
    // we have changes of sign at sc[n], between each we must seek for zeros

    // search in range [min,s[0]]
    smin = min;
    smax = Epos->data[(int) sc->data[0]];
    z = search_zero(numerov_integrate, &num, smin, eres, smax, getmaxabs(fn, 0, (int) sc->data[0]));
    if (z == z) {
      zp = array_append(zp, z);
    }

    // suche im Bereich [sn[i],sn[i+1]]
    for (int k = 0; k < sc->length - 1; k++) {
      double smin = Epos->data[(int) sc->data[k]], smax = Epos->data[(int) sc->data[k+1]];
      //printf("[%g,%g]\n", smin, smax);
      z = search_zero(numerov_integrate, &num, smin, eres, smax,
                      getmaxabs(fn, (int) sc->data[k], (int) sc->data[k+1]));
      if (z == z) {
        zp = array_append(zp, z);
      }
    }

    // suche im Bereich [s[n],max]
    smin = Epos->data[(int) sc->data[sc->length-1]];
    smax = max;
    z = search_zero(numerov_integrate, &num, smin, eres, smax,
      getmaxabs(fn, (int) sc->data[sc->length-1], fn->length-1));
    if (z == z) {
      zp = array_append(zp, z);
    }
  } else {
    // wir haben keine Nullstellen: suche im Bereich [min,max]
    z = search_zero(numerov_integrate, &num, min, eres, max,
      getmaxabs(fn, 0, fn->length-1));
    if (z == z) {
      zp = array_append(zp, z);
    }
  }

  free(sc);

  free(Epos);
  free(fn);
  return zp;
}