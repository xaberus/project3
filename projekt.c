#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include <assert.h>
#include <stdint.h>

#include <cairo.h>

/*************************************************************************************************/

typedef struct {
  int            bins;
  fftw_plan      fwd;
  fftw_plan      bwd;
  double       * V;
  fftw_complex * apsi;
  fftw_complex * eV;
  fftw_complex * ehV;
  fftw_complex * eT;
  fftw_complex * ehT;
  double         dt;
  double         dk;
  double         dx;
  double         L;
  double         min, max;
} splitop_t;

splitop_t * splitop_new(int bins, double dt, double min, double max, fftw_complex (*V)(double x)) {
  assert(bins > 1 && min < max);
  splitop_t * w = malloc(sizeof(splitop_t)); assert(w);
  w->bins = bins;
  w->min = min;
  w->max = max;
  w->L = max - min;
  w->dt = dt;
  w->dx = w->L / bins;
  w->dk = 2*M_PI/w->L;
  w->V = fftw_malloc(sizeof(fftw_complex) * bins); assert(w->V);
  for (int n = 0; n < bins; n++) {
    w->V[n] = V(min + w->dx*n);
  }
  w->eV = fftw_alloc_complex(bins); assert(w->eV);
  w->ehV = fftw_alloc_complex(bins); assert(w->ehV);
  for (int n = 0; n < bins; n++) {
    w->eV[n] = cexp(-I * w->V[n] * w->dt);
    w->ehV[n] = cexp(-I * w->V[n] * w->dt / 2);
  }
  w->ehT = fftw_alloc_complex(bins); assert(w->ehT);
  w->eT = fftw_alloc_complex(bins); assert(w->eT);
  for (int n = 0; n < bins; n++) {
    if (n < bins/2) {
      w->ehT[n] = cexp(-I * w->dk * w->dk * n * n * w->dt / 2);
      w->eT[n] = cexp(-I * w->dk * w->dk * n * n * w->dt);
    } else {
      int m = (n - bins);
      w->ehT[n] = cexp(-I * w->dk * w->dk * m * m * w->dt / 2);
      w->eT[n] = cexp(-I * w->dk * w->dk * m * m * w->dt);
    }
  }
  w->fwd = fftw_plan_dft_1d(bins, NULL, NULL, FFTW_FORWARD, FFTW_ESTIMATE);
  w->bwd = fftw_plan_dft_1d(bins, NULL, NULL, FFTW_BACKWARD, FFTW_ESTIMATE);
  return w;
}

void splitop_free(splitop_t * w) {
  fftw_destroy_plan(w->fwd);
  fftw_destroy_plan(w->bwd);
  fftw_free(w->V);
  fftw_free(w->eV);
  fftw_free(w->eT);
  fftw_free(w->ehT);
  fftw_free(w->apsi);
  free(w);
}

double cabssq(fftw_complex z) {
  double x = creal(z), y = cimag(z);
  return x * x + y * y;
}

double cfnormsq(int len, fftw_complex a[len]) {
  double A = 0;
  for (int k = 0; k < len; k++) A += cabssq(a[k]);
  return A;
}

#include <emmintrin.h>
#include <smmintrin.h>

/* a *= z */
inline static void cvect_mult_asign(int len, fftw_complex * a, fftw_complex * z) {
  for (int n = 0; n < len; n++, a++, z++) {
    double * va = (double *) a;
    double * vz = (double *) z;
    // (a0,a0)
    register __m128d tmp0 = _mm_load1_pd(va);
    // (b0,b1)
    register __m128d tmp1 = _mm_load_pd(vz);
    // (a1,a1)
    register __m128d tmp2 = _mm_load1_pd(va+1);
    // (b1,b0)
    register __m128d tmp3 = _mm_shuffle_pd(tmp1, tmp1, 1);
    // (a0 b0,a0 b1)
    register __m128d tmp4 = _mm_mul_pd(tmp0, tmp1);
    // (a1 b1,a1 b0)
    register __m128d tmp5 = _mm_mul_pd(tmp2, tmp3);
    register __m128d tmp6 = _mm_addsub_pd(tmp4, tmp5);
    _mm_store_pd(va, tmp6);
  }
}

void splitop_run(splitop_t * w, int times, fftw_complex * psi) {
  fftw_complex * psik = fftw_alloc_complex(w->bins); assert(psik);
  int bins = w->bins;
  fftw_complex * eV = w->eV;
  fftw_complex * ehV = w->eV;
  fftw_complex * ehT = w->ehT;
  fftw_complex * eT = w->eT;

  double c = sqrt(bins);

  /*for (int k = 0; k < times; k++) {
    //for (int n = 0; n < bins; n++) { psi[n] *= ehV[n]; }
    for (p = psi, as = ehV, ae = as + bins; as < ae; as++, p++) { (*p) *= *as; }
    fftw_execute_dft(w->fwd, psi, psik);
    //for (int n = 0; n < bins; n++) { psik[n] *= eT[n]; }
    for (p = psik, as = eT, ae = as + bins; as < ae; as++, p++) { (*p) *= *as; }
    fftw_execute_dft(w->bwd, psik, psi);
    //for (int n = 0; n < bins; n++) { psi[n] *= ehV[n] / bins; }
    for (p = psi, as = ehV, ae = as + bins; as < ae; as++, p++) { (*p) *= *as/bins; }
  }*/

  const double perm[2] = {-1,1};

  if (times > 0) {
    for (int n = 0; n < bins; n++) { psi[n] *= ehV[n]; }


    fftw_execute_dft(w->fwd, psi, psik);
    //for (int n = 0; n < bins; n++) { psik[n] *= eT[n]; }
    cvect_mult_asign(bins, psik, eT);

    for (int k = 0; k < times - 1; k++) {
      fftw_execute_dft(w->bwd, psik, psi);
      for (int n = 0; n < bins; n++) { psi[n] *= eV[n] / bins; }

      fftw_execute_dft(w->fwd, psi, psik);
      //for (int n = 0; n < bins; n++) { psik[n] *= eT[n]; }
      cvect_mult_asign(bins, psik, eT);
    }

    fftw_execute_dft(w->bwd, psik, psi);
    for (int n = 0; n < bins; n++) { psi[n] *= ehV[n] / bins; }
  }

  /*if (times > 0) {
    fftw_execute_dft(w->fwd, psi, psik);
    for (int n = 0; n < bins; n++) { psik[n] *= ehT[n] / c; }

    for (int k = 0; k < times - 1; k++) {
      fftw_execute_dft(w->bwd, psik, psi);
      for (int n = 0; n < bins; n++) { psi[n] *= eV[n] / c; }
      fftw_execute_dft(w->fwd, psi, psik);
      for (int n = 0; n < bins; n++) { psik[n] *= eT[n] / c; }
    }

    fftw_execute_dft(w->bwd, psik, psi);
    for (int n = 0; n < bins; n++) { psi[n] *= eV[n] / c; }
    fftw_execute_dft(w->fwd, psi, psik);
    for (int n = 0; n < bins; n++) { psik[n] *= ehT[n] / c; }

    fftw_execute_dft(w->bwd, psik, psi);
    for (int n = 0; n < bins; n++) { psi[n] /= c; }
  }*/
  fftw_free(psik);
}

fftw_complex * splitop_prepare(splitop_t * w, fftw_complex (*fn)(double x)) {
  fftw_complex * psi = fftw_alloc_complex(w->bins); assert(psi);
  for (int n = 0; n < w->bins; n++) {
    psi[n] = fn(w->min + w->dx * n);
  }
  double A = cfnormsq(w->bins, psi);
  for (int n = 0; n < w->bins; n++) {
    psi[n] *= 1/sqrt(A);
  }
  return psi;
}

void splitop_save(splitop_t * w, fftw_complex * spsi) {
  fftw_complex * psi = fftw_alloc_complex(w->bins); assert(psi);
  for (int n = 0; n < w->bins; n++) {
    psi[n] = spsi[n];
  }
  w->apsi = psi;
}

void splitop_draw(splitop_t * w, cairo_t * cr, cairo_rectangle_t rect, fftw_complex * psi) {
  cairo_save(cr);
  cairo_rectangle(cr, rect.x, rect.y, rect.width, rect.height);
  cairo_clip(cr);

  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_set_line_width(cr, 2);
  cairo_paint(cr);

  double y0 = rect.height/2;

  //cairo_set_line_width(cr, 10);
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_move_to(cr, rect.x, rect.y + y0);
  cairo_line_to(cr, rect.x + rect.width, rect.y + y0);
  cairo_stroke(cr);

  double xscale = rect.width/w->bins, yscale, max;

  max = 0;
  for (int n = 0; n < w->bins; n++) {
    double a = fabs(w->V[n]); if (max < a) { max = a; }
  }
  yscale = y0/max;

  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_move_to(cr, rect.x, y0);
  for (int n = 0; n < w->bins; n++) { cairo_line_to(cr, rect.x + n * xscale, rect.y + y0 - w->V[n] * yscale); }
  cairo_stroke(cr);

  fftw_complex * apsi = w->apsi;
  max = 0;
  for (int n = 0; n < w->bins; n++) {
    double a = cabs(apsi[n]); if (max < a) { max = a; }
  }
  yscale = y0/(5*max/4);

  cairo_set_line_width(cr, 1);
  cairo_set_source_rgb(cr, 0, 0, 1);
  cairo_move_to(cr, rect.x, y0);
  for (int n = 0; n < w->bins; n++) { cairo_line_to(cr, rect.x + n * xscale, rect.y + y0 - cabs(apsi[n]) * yscale); }
  cairo_stroke(cr);

  cairo_set_line_width(cr, 2);
  cairo_set_source_rgb(cr, 1, 0, 0);
  cairo_move_to(cr, rect.x, y0);
  for (int n = 0; n < w->bins; n++) { cairo_line_to(cr, rect.x + n * xscale, rect.y + y0 - cabs(psi[n]) * yscale); }
  cairo_stroke(cr);

  fftw_complex * psik = fftw_alloc_complex(w->bins); assert(psik);
  fftw_execute_dft(w->fwd, psi, psik);

  max = 0;
  for (int n = 0; n < w->bins; n++) {
    double a = cabs(psik[n]); if (max < a) { max = a; }
  }
  yscale = y0/(5*max/4);

  cairo_set_line_width(cr, 1);
  cairo_set_source_rgb(cr, 1, .5, 0);
  cairo_move_to(cr, rect.x, y0);
  for (int n = 0; n < w->bins; n++) {
    int l = (n+w->bins/2)%w->bins;
    cairo_line_to(cr, rect.x + (n-w->bins/3)*4 * xscale, rect.y + 2*y0 - cabs(psik[l]) * yscale);
  }
  cairo_stroke(cr);

  fftw_free(psik);

  cairo_restore(cr);
}

/*************************************************************************************************/

double pota = 0.9;
//double pota = 3;
double potb = 1;
double potc = 30;

/*double pota = 2;
double potb = 1;*/

fftw_complex zeropot(double x) {
  double b = potb * potb;
  double l = x + pota;
  double r = x - pota;
  return -potc*exp(-l*l/b)-potc*exp(-r*r/b);

  //return 40*(-1/sqrt(1+potb*pow(x-pota,2))-1/sqrt(1+potb*pow(x+pota,2)));

  /*double a = pota;
  double b = 2;
  double r = fabs(x);
  if (r > fabs(a/2-b/2) && r < a/2+b/2) return -100;
  if (r < fabs(a/2-b/2)) return -30;
  //if (r > 25) return 1000;
  return 0;*/
}

fftw_complex fn(double x) {
  double a = .5;
  double aa = a * a;

  double x0 = -pota;

  //double x0 = -1.941063416246160;

  double k0 = 0;
  double xx = (x-x0) * (x-x0);
  //return  exp(-xx/(aa))*cexp(I*k0*(x));
  return  exp(-xx/(aa));
}

int main(/*int argc, char *argv[]*/) {
  splitop_t * sop = splitop_new(4096, .0001, -25, 25, zeropot);

  fftw_complex * psi = splitop_prepare(sop, fn);

  splitop_save(sop, psi);

  //splitop_run(sop, 1, psi);

  int width = 1000, height = 400;
  char buf[256];

  cairo_surface_t * surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
  cairo_t *cr = cairo_create(surface);

  /*splitop_draw(sop, cr, (cairo_rectangle_t){0, 0, width, height}, psi);
  snprintf(buf, sizeof(buf), "image%05u.png", 0);
  cairo_surface_write_to_png(surface, buf);*/

  for (int k = 1; k < 220; k++) {
    //fprintf(stderr, "##### %g, %u\n", cfnormsq(sop->bins, psi), k);
    fprintf(stderr, "##### %u\r", k);
    splitop_run(sop, 150, psi);
    /*splitop_draw(sop, cr, (cairo_rectangle_t){0, 0, width, height}, psi);
    snprintf(buf, sizeof(buf), "image%05u.png", k);
    cairo_surface_write_to_png(surface, buf);*/
  }

  fprintf(stderr, "done---------\n");

  cairo_destroy(cr);
  cairo_surface_destroy(surface);

  fftw_free(psi);
  splitop_free(sop);
  fftw_cleanup();
  return 0;
}