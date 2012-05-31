#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "splitop.h"
#include "cvect.h"

splitop_t * splitop_new(preferences_t * prefs) {
  splitop_t * w = malloc(sizeof(splitop_t)); assert(w);
  w->prefs = prefs;

  double * V = prefs->potential->data;
  int bins = prefs->bins;
  double dt = prefs->dt;
  double dk = prefs->dk;

  w->eV = fftw_alloc_complex(bins); assert(w->eV);
  w->eVn = fftw_alloc_complex(bins); assert(w->eVn);
  w->ehV = fftw_alloc_complex(bins); assert(w->ehV);
  w->ehVn = fftw_alloc_complex(bins); assert(w->ehVn);
  w->eT = fftw_alloc_complex(bins); assert(w->eT);
  w->apsi = fftw_alloc_complex(bins); assert(w->apsi);
  w->psi = fftw_alloc_complex(bins); assert(w->psi);
  w->psik = fftw_alloc_complex(bins); assert(w->psik);

  for (int n = 0; n < bins; n++) {
    w->eV[n] = cexp(-I * V[n] * dt);
    w->eVn[n] = w->eV[n] / bins;
    w->ehV[n] = cexp(-I * V[n] * dt / 2);
    w->ehVn[n] = w->ehV[n] / bins;
  }
  for (int n = 0; n < bins; n++) {
    if (n < bins/2) {
      w->eT[n] = cexp(-I * dk * dk * n * n * dt);
    } else {
      int m = (n - bins);
      w->eT[n] = cexp(-I * dk * dk * m * m * dt);
    }
  }


  /*w->fwd = fftw_plan_dft_1d(bins, NULL, NULL, FFTW_FORWARD, FFTW_ESTIMATE);
  w->bwd = fftw_plan_dft_1d(bins, NULL, NULL, FFTW_BACKWARD, FFTW_ESTIMATE);*/
  w->fwd = fftw_plan_dft_1d(bins, w->psi, w->psik, FFTW_FORWARD, FFTW_MEASURE);
  w->bwd = fftw_plan_dft_1d(bins, w->psik, w->psi, FFTW_BACKWARD, FFTW_MEASURE);

  for (int k = 0; k < bins; k++) {
    w->psi[k] = prefs->psi->data[k];
  }

  return w;
}

void splitop_free(splitop_t * w) {
  fftw_destroy_plan(w->fwd);
  fftw_destroy_plan(w->bwd);
  fftw_free(w->eV);
  fftw_free(w->eVn);
  fftw_free(w->ehV);
  fftw_free(w->ehVn);
  fftw_free(w->eT);
  fftw_free(w->apsi);
  fftw_free(w->psi);
  fftw_free(w->psik);
  free(w);
}

void splitop_run(splitop_t * w, int times) {
  int bins = w->prefs->bins;
  fftw_complex * ehV = w->ehV;
  fftw_complex * eVn = w->eVn;
  fftw_complex * ehVn = w->ehVn;
  fftw_complex * eT = w->eT;

  fftw_complex * psi = w->psi;
  fftw_complex * psik = w->psik;

  /*fftw_complex * p, * as, * ae;
  for (int k = 0; k < times; k++) {
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

}

void splitop_prepare(splitop_t * w) {
  int bins = w->prefs->bins;
  fftw_complex * psi = w->psi;
  double A = cvect_normsq(bins, psi);
  for (int n = 0; n < bins; n++) {
    psi[n] *= 1/sqrt(A);
  }
  return psi;
}

void splitop_save(splitop_t * w) {
  int bins = w->prefs->bins;
  fftw_complex * psi = w->psi;
  fftw_complex * apsi = w->apsi;
  for (int n = 0; n < bins; n++) {
    apsi[n] = psi[n];
  }
}

void splitop_restore(splitop_t * w) {
  int bins = w->prefs->bins;
  fftw_complex * psi = w->psi;
  fftw_complex * apsi = w->apsi;
  for (int n = 0; n < bins; n++) {
    psi[n] = apsi[n];
  }
}