#ifndef _CVECT_H_
#define _CVECT_H_

#include <complex.h>
#include <emmintrin.h>
#include <smmintrin.h>

/*! \class cvect
 a pointer to a complex array with impicit length */

/*! \memberof cvect
 calculates the scalar product of \a a with itself */
inline static
double cvect_normsq(int len, complex double * a)
{
  double A = 0;
  for (int k = 0; k < len; k++) {
    complex double z = *(a++);
    double x = creal(z), y = cimag(z);
    A += x * x + y * y;
  }
  return A;
}

/*! \memberof cvect
 calculates the scalar product of \a a and \a b */
inline static
fftw_complex cvect_skp(int len, complex double * a, complex double * b)
{
  fftw_complex c = 0;
  for (; len > 0 ; len--, a++, b++) {
    c += conj(*a)* *b;
  }
  return c;
}

/*! \memberof cvect
 the operation of this function is to assign to each member of \a a
 the value of a[i]*z[i]

 as this function is one of the bottlenecks in the computation
 sse4.1 intrinsics are employes to help the compiler...
 */
inline static
void cvect_mult_asign(int len, complex double * a, complex double * z)
{
  /*! the plain code for this function would be
    \code
      for (int n = 0; n < len; n++) { a[n] *= z[n]; }
    \endcode
  */
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