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

#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <lua.h>

#include "array.h"
#include "carray.h"

/*! results struct */
typedef struct results {
  carray_t * co;    /**< unmodified correlation function */
  carray_t * c;     /**< correlation function with von-hann window applied */
  carray_t * ck;    /**< modulated spectrum (DFT of c) */
} results_t;

results_t * results_new();
void results_free(results_t * res);

typedef struct preferences {
  int        config;    /**< reference to config table in lua registry */
  int        bins;      /**< number of sampling points */
  double     dt;        /**< time delta */
  struct {
    double min;         /**< minimal x position */
    double max;         /**< maximal x position */
  } range;              /**< position space range */
  int        steps;     /**< steps to run between two evaluations of autocorrelation */
  int        runs;      /**< number of evaluations of autocorrelation */
  int        vstep;     /**< steps to run between two frames */
  int        vframes;   /**< number of rames to render */
  array_t  * xpos;      /**< position space sampling points */
  array_t  * potential; /**< pontential samples at position sampling points */
  carray_t * psi;       /**< system sate samples at position sampling points */

  struct {
    char * dir;         /**< directory to save outputs to */
    char * apsi;        /**< filename for initial system state sump */
    char * pot;         /**< filename for potential dump */
    char * corr;        /**< filename for autorcorrelation function dump */
    char * dftcorr;     /**< filename for dump DTF of corr */
    char * theoenrg;    /**< filename to dump theoretical energies to */
    char * spectrum;    /**< filename to dump spectrum to */
  } output;             /**< output configuration */

  int        tsteps;    /**< total number of iterations */
  double     dx;        /**< position delta */
  double     dk;        /**< momentum delta */
  double     dE;        /**< energy delty */

  struct {
    double min;         /**< minimal energy to consider */
    double max;         /**< maximal energy to consider */
    double win;         /**< maximal energy to consider */
  } enrgrange;          /**< hints for gnuplot */

  array_t  * theoenrg;  /**< theoretical energies from config */

  results_t * results;  /**< simulation results */
} preferences_t;

void preferences_free(preferences_t * prefs);
preferences_t * preferences_new();
int preferences_read(lua_State * L, preferences_t * prefs);
int start_simulation(preferences_t * prefs);
int dump_results(preferences_t * prefs);

#endif /* _SIMULATION_H_ */