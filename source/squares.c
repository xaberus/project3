#include <stdio.h>
#include <assert.h>

#include "squares.h"

double squares_sum_x(double x, sqparam_t * par)
{
  double sum = 0;
  double a0 = par->a0;
  int l = par->l;
  double * s = par->s;
  double * f = par->f;
  double (*fn)(double x, double a) = par->fn;
  for (int k = 0; k < l; k++) {
    double z = fn(s[k] - x, a0) - f[k];
    sum += z * z;
  }
  return sum;
}

double squares_sum_x_dx(double x, sqparam_t * par)
{
  double sum = 0;
  double a0 = par->a0;
  int l = par->l;
  double * s = par->s;
  double * f = par->f;
  double (*fn)(double x, double a) = par->fn;
  double (*dx)(double x, double a) = par->dx;
  for (int k = 0; k < l; k++) {
    double z = s[k] - x;
    sum += - 2 * (fn(z, a0) - f[k]) * dx(z, a0);
  }
  return sum;
}

double squares_sum_x_ddx(double x, sqparam_t * par)
{
  double sum = 0;
  double a0 = par->a0;
  int l = par->l;
  double * s = par->s;
  double * f = par->f;
  double (*fn)(double x, double a) = par->fn;
  double (*dx)(double x, double a) = par->dx;
  double (*ddx)(double x, double a) = par->ddx;
  for (int k = 0; k < l; k++) {
    double z = s[k] - x;
    sum += 2 * dx(z, a0) + 2 * (fn(z, a0) - f[k]) * ddx(z, a0);
  }
  return sum;
}

double squares_sum_a(double a, sqparam_t * par)
{
  double sum = 0;
  double x0 = par->x0;
  int l = par->l;
  double * s = par->s;
  double * f = par->f;
  double (*fn)(double x, double a) = par->fn;
  for (int k = 0; k < l; k++) {
    double z = fn(s[k] - x0, a) - f[k];
    sum += z * z;
  }
  return sum;
}

double squares_sum_a_da(double a, sqparam_t * par)
{
  double sum = 0;
  double x0 = par->x0;
  int l = par->l;
  double * s = par->s;
  double * f = par->f;
  double (*fn)(double x, double a) = par->fn;
  double (*da)(double x, double a) = par->da;
  for (int k = 0; k < l; k++) {
    double z = s[k] - x0;
    sum += 2 * (fn(z, a) - f[k]) * da(z, a);
  }
  return sum;
}

double squares_sum_a_dda(double a, sqparam_t * par)
{
  double sum = 0;
  double x0 = par->x0;
  int l = par->l;
  double * s = par->s;
  double * f = par->f;
  double (*fn)(double x, double a) = par->fn;
  double (*da)(double x, double a) = par->da;
  double (*dda)(double x, double a) = par->dda;
  for (int k = 0; k < l; k++) {
    double z = s[k] - x0;
    sum += 2 * da(z, a) + 2 * (fn(z, a) - f[k]) * dda(z, a);
  }
  return sum;
}

static const int maximal_iterations = 1000;

double squares_search_zero(
  sqparam_t * par,
  double fn(double a, sqparam_t * par),
  double dfn(double a, sqparam_t * par),
  double min, double max)
{
  double a, b;
  //printf("squares_search_zero([%g:%g])\n", min, max);
  /*if (fn == squares_sum_x_dx) {
    FILE * fp = fopen("suck", "w");
    for (double z = min; z < max; z += (max-min)/1000) {
      fprintf(fp, "%.17e %.17e %.17e %.17e\n", z,
        squares_sum_x(z, par),
        squares_sum_x_dx(z, par),
        squares_sum_x_ddx(z, par));
    }
    fclose(fp);
  } else {
    FILE * fp = fopen("sack", "w");
    for (double z = min; z < max; z += (max-min)/1000) {
      fprintf(fp, "%.17e %.17e %.17e %.17e\n", z,
        squares_sum_a(z, par),
        squares_sum_a_da(z, par),
        squares_sum_a_dda(z, par));
    }
    fclose(fp);
  }*/
  double x = min, xp = 0;
  int k;
  for (k = 0; x >= min && x <= max && k < maximal_iterations; k++) {
    xp = x - fn(x, par)/dfn(x, par);
    if (xp > max || xp < min) {
      //printf("bisect (range)\n");
      goto bisect;
    }
    if (x == xp) {
      return x;
    }
    x = xp;
    continue;
bisect:
      // wenn wir hier landen, dann war x zu nah an einem Extrempunkt, versuche einen besseren
      // Startwert mit Bisektionsverfahren zu finden
      a = fn(min, par);
      b = fn(max, par);
      if (a * b < 0) {
        double t = (min + max) / 2.;
        double w = fn(t, par);
        //printf("bis: [%.17e:%.17e:%.17e] ~ {%g, %g, %g}\n", min, t, max, a, w, b);
        if (w == 0) {
          return t;
        }
        if (a * w < 0) {
          max = t;
        } else {
          min = t;
        }
        x = min;
        continue;
      }
      break;
  }
  if ((x - xp) < 0.000000001) return x;
  return 0.0/0.0; // nan - nichts gefunden
}

double squares_search_min(
  sqparam_t * par,
  double fn(double z, sqparam_t * par),
  double dfn(double z, sqparam_t * par),
  double ddfn(double z, sqparam_t * par),
  double min, double max)
{
  (void) fn;
  double z;
  int next = 0;
  int iterations = 10;
  do {
    next = 0;
    z = squares_search_zero(par, dfn, ddfn, min, max);
    if (z == z) {
      if(ddfn(z, par) <= 0) {
        //printf("#### %.17e\n", z);
        min = z + (max - min) / 100000.;
        next = 1;
      }
    }
  } while (next && iterations-- && min < max);
  if (next)
    return 0.0/0.0;
  return z;
}

void straw(sqparam_t * par)
{
  FILE * fp = fopen("straw", "w"); assert(fp);
  double min = par->s[0];
  double max = par->s[par->l - 1];
  for (double z = min; z < max; z += (max-min)/500) {
    fprintf(fp, "%g %g\n", z, par->fn(z - par->x0, par->a0));
  }
  fclose(fp);
}

double squares_min_2d(sqparam_t * par)
{
  double xmin = par->xmin;
  double xmax = par->xmax;
  double amin = par->amin;
  double amax = par->amax;
  printf("squares_min_2d([%g:%g]) :: ", xmin, xmax);
  int iterations = 1;
  do {
    double a = squares_search_min(par,
      squares_sum_a,
      squares_sum_a_da,
      squares_sum_a_dda, amin, amax);
    if (a == a) {
      par->a0 = a;
      //printf("found a: %.17e->%.17e!\n", par->a0, a); straw(par);
      double x = squares_search_min(par,
        squares_sum_x,
        squares_sum_x_dx,
        squares_sum_x_ddx, xmin, xmax);
      if (x == x) {
        //printf("found x: %.17e->%.17e!\n", par->x0, x); straw(par);
        if (x == par->x0) {
          printf("converged after %d iterations -> %.17e\n", iterations, x);
          return x;
        }
        par->x0 = x;
        iterations++;
        continue;
      }
    }
    printf("failed after %d! ~ %.17e (probably not a peak)\n", iterations, par->x0);
    //fflush(stdout); assert(0);
    return 0.0/0.0;
  } while (iterations < 100);
  printf("(!) not converged after %d iterations -> %.17e\n", iterations, par->x0);
  return par->x0;
}