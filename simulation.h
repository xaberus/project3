#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <lua.h>

#include "array.h"
#include "carray.h"

/*! results struct */
typedef struct retults {
  carray_t * co;
  carray_t * c;
  carray_t * ck;
} results_t;

results_t * results_new();
void results_free(results_t * res);

typedef struct preferences {
  int        config;
  int        bins;
  double     dt;
  struct {
    double min;
    double max;
  } range;
  int        steps;
  int        runs;
  int        vstep;
  int        vframes;
  array_t  * xpos;
  array_t  * potential;
  carray_t * psi;

  struct {
    char * dir;
    char * apsi;
    char * pot;
    char * corr;
    char * dftcorr;
    char * theoenrg;
  } output;

  int        tsteps;
  double     dx;
  double     dk;
  double     dE;

  struct {
    double min;
    double max;
  } enrgrange;

  array_t  * theoenrg;

  results_t * results;
} preferences_t;

void preferences_free(preferences_t * prefs);
preferences_t * preferences_new();
int preferences_read(lua_State * L, preferences_t * prefs);
int start_simulation(preferences_t * prefs);
int dump_results(preferences_t * prefs);

#endif /* _SIMULATION_H_ */