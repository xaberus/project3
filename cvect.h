#ifndef _CVECT_H_
#define _CVECT_H_

#include <complex.h>
#include <emmintrin.h>
#include <smmintrin.h>

inline static
double cvect_normsq(int len, complex double a[len]) {
  double A = 0;
  for (int k = 0; k < len; k++) {
    complex double z = a[k];
    double x = creal(z), y = cimag(z);
    A += x * x + y * y;
  }
  return A;
}

inline static
fftw_complex cvect_skp(int len, complex double * a, complex double * b) {
  fftw_complex c = 0;
  for (; len > 0 ; len--, a++, b++) {
    c += conj(*a)* *b;
  }
  return c;
}

/* a *= z */
inline static
void cvect_mult_asign(int len, complex double * a, complex double * z) {
  //for (int n = 0; n < len; n++) { a[n] *= z[n]; }
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

#endif /* _CVECT_H_ */