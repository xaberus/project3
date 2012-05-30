#ifndef _SPLITOP_H_
#define _SPLITOP_H_

#include <complex.h>
#include <fftw3.h>

typedef struct {
  int            bins;
  fftw_plan      fwd;
  fftw_plan      bwd;
  double       * V;
  fftw_complex * apsi;
  fftw_complex * eV;
  fftw_complex * eVn;
  fftw_complex * ehV;
  fftw_complex * ehVn;
  fftw_complex * eT;
  double         dt;
  double         dk;
  double         dx;
  double         L;
  double         min, max;
} splitop_t;

splitop_t * splitop_new(int bins, double dt, double min, double max, fftw_complex (*V)(double x));
void splitop_free(splitop_t * w);

void splitop_save(splitop_t * w, fftw_complex * spsi);
fftw_complex * splitop_prepare(splitop_t * w, fftw_complex (*fn)(double x));
void splitop_run(splitop_t * w, int times, fftw_complex * psi);

#endif /* _SPLITOP_H_ */