#ifndef _SPLITOP_H_
#define _SPLITOP_H_

#include <complex.h>
#include <fftw3.h>

#include "simulation.h"

/*! split step operator class */
typedef struct splitop {
  preferences_t * prefs; /**< preferences used to create this operator */

  fftw_plan      fwd;    /**< plan for forward fourier transform */
  fftw_plan      bwd;    /**< plan for bachward fourier transform */

  fftw_complex * apsi;   /**< reference wave function */

  fftw_complex * psi;    /**< current system state */
  fftw_complex * psik;   /**< temporary used for calculations in momentum space */

  fftw_complex * eV;    /**< U_V position space */
  fftw_complex * eVn;   /**< U_V/bins in position space */
  fftw_complex * ehV;   /**< U_{V/2} in position space */
  fftw_complex * ehVn;  /**< U_{V/2}/bins in position space */

  fftw_complex * eT;    /**< U_T in momentum space */
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