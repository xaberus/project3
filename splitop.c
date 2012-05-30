#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "splitop.h"
#include "cvect.h"

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
  w->eVn = fftw_alloc_complex(bins); assert(w->eVn);
  w->ehV = fftw_alloc_complex(bins); assert(w->ehV);
  w->ehVn = fftw_alloc_complex(bins); assert(w->ehVn);
  for (int n = 0; n < bins; n++) {
    w->eV[n] = cexp(-I * w->V[n] * w->dt);
    w->eVn[n] = w->eV[n] / bins;
    w->ehV[n] = cexp(-I * w->V[n] * w->dt / 2);
    w->ehVn[n] = w->ehV[n] / bins;
  }
  w->eT = fftw_alloc_complex(bins); assert(w->eT);
  for (int n = 0; n < bins; n++) {
    if (n < bins/2) {
      w->eT[n] = cexp(-I * w->dk * w->dk * n * n * w->dt);
    } else {
      int m = (n - bins);
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
  fftw_free(w->eVn);
  fftw_free(w->ehV);
  fftw_free(w->ehVn);
  fftw_free(w->eT);
  fftw_free(w->apsi);
  free(w);
}

void splitop_run(splitop_t * w, int times, fftw_complex * psi) {
  fftw_complex * psik = fftw_alloc_complex(w->bins); assert(psik);
  int bins = w->bins;
  fftw_complex * ehV = w->ehV;
  fftw_complex * eVn = w->eVn;
  fftw_complex * ehVn = w->ehVn;
  fftw_complex * eT = w->eT;

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

  if (times > 0) {
    cvect_mult_asign(bins, psi, ehV);
    fftw_execute_dft(w->fwd, psi, psik);
    cvect_mult_asign(bins, psik, eT);
    for (int k = 0; k < times - 1; k++) {
      fftw_execute_dft(w->bwd, psik, psi);
      cvect_mult_asign(bins, psi, eVn);
      fftw_execute_dft(w->fwd, psi, psik);
      cvect_mult_asign(bins, psik, eT);
    }
    fftw_execute_dft(w->bwd, psik, psi);
    cvect_mult_asign(bins, psi, ehVn);
  }

  fftw_free(psik);
}

fftw_complex * splitop_prepare(splitop_t * w, fftw_complex (*fn)(double x)) {
  fftw_complex * psi = fftw_alloc_complex(w->bins); assert(psi);
  for (int n = 0; n < w->bins; n++) {
    psi[n] = fn(w->min + w->dx * n);
  }
  double A = cvect_normsq(w->bins, psi);
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