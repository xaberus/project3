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

int maximal_iterations = 1000;

typedef struct {
  array_t * V;
  double h;
  double psi1;
  double psi2;
} numerov_t;

double numerov_integrate(double E, void * arg)
{
  numerov_t * num = arg;
  int length = num->V->length;
  double * V = num->V->data;
  double psi[3] = {num->psi1, num->psi2, 0};
  double hh = num->h * num->h;
  int k, p;
  for (k = 2, p = 0; k < length; k++, p++) {
    double f1 = psi[p % 3], f2 = psi[(p + 1) % 3];
    double s2 = 2 * (1 - 5 * hh / 12 * (E - V[k - 1])) * f2;
    double s1 = (1 + hh / 12 * (E - V[k - 2])) * f1;
    double d = (1 + hh / 12 * (E - V[k]));
    psi[(p + 2) % 3] = (s2 - s1) / d;
  }
  return psi[(p + 2) % 3];
}

/* newton: secant method
 * seek for zero in the intervall, with wt most one change of sign of first derrivative of fn */
double search_zero(double (*fn)(double, void*), void * arg, double min, double step, double max, double fmax) {
  double xm = min, x = xm + step, xp;
  int k;
  for (k = 0; x >= min && x <= max && k < maximal_iterations; k++) {
    double f = fn(x, arg)/fmax;
    double fm = fn(xm, arg)/fmax;
    double df = (f-fm);
    if (df == 0.0) { // worst case scenario....
      df = 10e-17;
    }
    //printf("->sz: {%.17g, %.17g, %.17g, %.17g}\n", f, fm, xm, x);
    xp = x - (x - xm)/df*f;
    //printf("sz: [%.17g;%.17g]:%.17g {%.17g, %.17g, %.17g}\n", min, max, step, xm, x, xp);
    //printf("sz: [%g,%g]:%g {%g, %g, %g}\n", min, max, step, xm, x, xp);
    if (xp > max || xp < min) {
      // wenn wir hier landen, dann war x zu nah an einem Extrempunkt, versuche einen besseren
      // Startwert mit Bisektionsverfahren zu finden
      double a = fn(min, arg)/fmax;
      double b = fn(max, arg)/fmax;
      if (a * b < 0) {
        double t = (min + max)/2;
        double w = fn(t, arg)/fmax;
        if (w == 0) {
          return w;
        }
        if (a * w < 0) {
          max = t;
        } else {
          min = t;
        }
        xm = min;
        x = xm + step;
        continue;
      }
      break;
    }
    if (x == xp) {
      return x;
    }
    xm = x;
    x = xp;
  }
  return 0.0/0.0; // nan - nichts gefunden
}

double getmaxabs(array_t * fn, int start, int end)
{
  double max = -1/0.0;
  double * f = fn->data;
  for (int k = start; k < end; k++) {
    double z = fabs(*(f++));
    if (z > max) {
      max = z;
    }
  }
  return max;
}

array_t * numerov_energies(preferences_t * prefs)
{

  array_t * zp = array_new_sized(0, 100);
  double smin, smax, z = 0;

  int res = 4096;

  numerov_t num = {prefs->potential, prefs->dx, 0, 10e-16};

  double min = prefs->enrgrange.min;
  double max = prefs->enrgrange.max;

  array_t * Epos = array_equipart(min, max, res);
  array_t * fn = array_mapv(Epos, numerov_integrate, &num);

  array_t * sc = search_der_sign_change_3(Epos->data[1] - Epos->data[0], fn, 0, 0, 0);

  array_dump_to_file("score", " ", 2, Epos, fn);
  FILE * fp = fopen("sumo", "w");
  for (int k = 0; k < sc->length - 1; k++) {
    fprintf(fp, "%g\n", Epos->data[(int)sc->data[k]]);
  }
  fclose(fp);

  if (sc->length > 0) {
    // we have changes of sign at sc[n], between each we must seek for zeros

    // search in range [min,s[0]]
    smin = min;
    smax = Epos->data[(int) sc->data[0]];
    z = search_zero(numerov_integrate, &num, smin, (smax - smin)/res, smax, getmaxabs(fn, 0, (int) sc->data[0]));
    if (z == z) {
      zp = array_append(zp, z);
    }

    // suche im Bereich [sn[i],sn[i+1]]
    for (int k = 0; k < sc->length - 1; k++) {
      double smin = Epos->data[(int) sc->data[k]], smax = Epos->data[(int) sc->data[k+1]];
      //printf("[%g,%g]\n", smin, smax);
      z = search_zero(numerov_integrate, &num, smin, (smax - smin)/res, smax,
                      getmaxabs(fn, (int) sc->data[k], (int) sc->data[k+1]));
      if (z == z) {
        zp = array_append(zp, z);
      }
    }

    // suche im Bereich [s[n],max]
    smin = Epos->data[(int) sc->data[sc->length-1]];
    smax = max;
    z = search_zero(numerov_integrate, &num, smin, (smax - smin)/res, smax,
      getmaxabs(fn, (int) sc->data[sc->length-1], fn->length-1));
    if (z == z) {
      zp = array_append(zp, z);
    }
  } else {
    // wir haben keine Nullstellen: suche im Bereich [min,max]
    z = search_zero(numerov_integrate, &num, min, (max - min)/res, max,
      getmaxabs(fn, 0, fn->length-1));
    if (z == z) {
      zp = array_append(zp, z);
    }
  }

  free(Epos);
  free(fn);
  free(sc);
  return zp;
}