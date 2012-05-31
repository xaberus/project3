#ifndef _SPLITOP_H_
#define _SPLITOP_H_

#include <complex.h>
#include <fftw3.h>

#include "simulation.h"

typedef struct {
  preferences_t * prefs;

  fftw_plan      fwd;
  fftw_plan      bwd;

  fftw_complex * apsi;

  fftw_complex * psi;
  fftw_complex * psik;

  fftw_complex * eV;
  fftw_complex * eVn;
  fftw_complex * ehV;
  fftw_complex * ehVn;

  fftw_complex * eT;
} splitop_t;

splitop_t * splitop_new(preferences_t * prefs);
void splitop_free(splitop_t * w);

void splitop_save(splitop_t * w);
void splitop_restore(splitop_t * w);
void splitop_prepare(splitop_t * w);
void splitop_run(splitop_t * w, int times);

#ifdef USE_CAIRO
#include <cairo.h>
void splitop_draw(splitop_t * w, cairo_t * cr, cairo_rectangle_t rect, fftw_complex * psi);
#endif

#endif /* _SPLITOP_H_ */