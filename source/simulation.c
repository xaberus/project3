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

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

#include <lauxlib.h>
#include <lualib.h>

#ifdef USE_CAIRO
# include <cairo.h>
#endif

#include "splitop.h"
#include "cvect.h"
#include "array.h"
#include "carray.h"
#include "peaks.h"

#include "numerov.h"

#include "simulation.h"

#include "akima.h"

/*! \memberof results
 creates an empty results struct */
results_t * results_new() {
  results_t * res = malloc(sizeof(results_t));
  memset(res, 0, sizeof(results_t ));
  return res;
}

/*! \memberof results
 releases memory used by results */
void results_free(results_t * res) {
  free(res->c);
  free(res->co);
  free(res->ck);
  free(res);
}

/*! \memberof preferences
 allocates a new preferences struct */
preferences_t * preferences_new() {
  preferences_t * prefs = malloc(sizeof(preferences_t));
  memset(prefs, 0, sizeof(preferences_t));
  return prefs;
}

/*! \memberof preferences
 releases memory used by preferences */
void preferences_free(preferences_t * prefs) {
  free(prefs->xpos);
  free(prefs->potential);
  free(prefs->psi);

  free(prefs->output.dir);
  free(prefs->output.apsi);
  free(prefs->output.pot);
  free(prefs->output.corr);
  free(prefs->output.dftcorr);
  free(prefs->output.spectrum);
  free(prefs->output.numen);
  free(prefs->output.splen);
  free(prefs->output.aken);

  if (prefs->results) {
    results_free(prefs->results);
  }

  free(prefs);
}

/*! \memberof preferences
 * this function reads in the whole configuration needed for the simulation
 */
int preferences_read(lua_State * L, preferences_t * prefs)
{
  // reset lua stack, get config table from registry and save its position on the stack
  lua_settop(L, 0);
  lua_rawgeti(L, LUA_REGISTRYINDEX, prefs->config);
  int config = lua_gettop(L);

  printf("configuration: \n");

  // get the number of bins for position space discretization */
  lua_getfield(L, config, "bins");
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, bins undefined\n"); return -1; }
  prefs->bins = lua_tointeger(L, -1); lua_pop(L, 1);
  printf("  bins:    %d\n", prefs->bins);

  // get time discretization delta
  lua_getfield(L, config, "dt");
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, dt undefined\n"); return -1; }
  prefs->dt = lua_tonumber(L, -1); lua_pop(L, 1);
  printf("  dt:      %g\n", prefs->dt);

  // get position space range [xmin:xmax]
  {
    lua_getfield(L, config, "range");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, range undefined\n"); return -1; }
    int range = lua_gettop(L);
    lua_rawgeti(L, range, 1);
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, range.min undefined\n"); return -1; }
    prefs->range.min = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_rawgeti(L, range, 2);
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, range.max undefined\n"); return -1; }
    prefs->range.max = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_pop(L, 1);
  }
  printf("  range:   [%g;%g]\n", prefs->range.min, prefs->range.max);

  // get number of steps to preform in each iteration
  lua_getfield(L, config, "steps");
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, steps undefined\n"); return -1; }
  prefs->steps = lua_tointeger(L, -1); lua_pop(L, 1);
  printf("  steps:   %d\n", prefs->steps);

  // get number of iterations to run
  lua_getfield(L, config, "runs");
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, runs undefined\n"); return -1; }
  prefs->runs = lua_tointeger(L, -1); lua_pop(L, 1);
  printf("  runs:    %d\n", prefs->runs);

#ifdef USE_CAIRO
  // get number of steps to preform for each frame
  lua_getfield(L, config, "vstep");
  if (!lua_isnil(L, -1)) {
    prefs->vstep = lua_tointeger(L, -1);
  } else {
    prefs->vstep = -1;
  }
  lua_pop(L, 1);
  // get number of frames to render
  lua_getfield(L, config, "vframes");
  if (!lua_isnil(L, -1)) {
    prefs->vframes = lua_tointeger(L, -1);
  } else {
    prefs->vframes = -1;
  }
  lua_pop(L, 1);
  if (prefs->vstep > 0 && prefs->vframes > 0) {
    printf("  will render video with:\n");
    printf("    vstep:       %d\n", prefs->vstep);
    printf("    vframes:     %d\n", prefs->vframes);
  }
#endif

  // equipatition of position space
  array_t * xpos = array_equipart(prefs->range.min, prefs->range.max, prefs->bins);
  prefs->xpos = xpos;

  // evaluate potential function from config for each position
  {
    lua_getfield(L, config, "potential");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, potential undefined\n"); return -1; }
    int potential = lua_gettop(L);

    array_t * V = array_new(xpos->length);

    for (int k = 0; k < xpos->length; k++) {
      lua_pushvalue(L, potential);
      lua_pushnumber(L, xpos->data[k]);
      lua_call(L, 1, 1);
      assert(lua_isnumber(L, -1));
      V->data[k] = lua_tonumber(L, -1);
      lua_pop(L, 1);
      //printf("potential(%g) = %g\n", xpos->data[k], V->data[k]);
    }
    lua_pop(L, 1);
    prefs->potential = V;
  }
  printf("  potential(x)\n");

  // evaluate wavefunction  from config for each position
  {
    lua_getfield(L, config, "psi");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, psi undefined\n"); return -1; }
    int psi = lua_gettop(L);

    carray_t * F = carray_new(xpos->length);

    for (int k = 0; k < xpos->length; k++) {
      lua_pushvalue(L, psi);
      lua_pushnumber(L, xpos->data[k]);
      lua_call(L, 1, 1);
      assert(lua_istable(L, -1));
      int tab = lua_gettop(L);
      lua_rawgeti(L, tab, 1);
      double x = lua_tonumber(L, -1);
      lua_pop(L, 1);
      lua_rawgeti(L, tab, 2);
      double y = lua_tonumber(L, -1);
      lua_pop(L, 2);
      F->data[k] = x + I * y;
      //printf("psi(%g) = %g + i %g\n", xpos->data[k], x, y);
    }
    lua_pop(L, 1);
    prefs->psi = F;
  }
  printf("  psi(x)\n");

  // get filenames for output
  {
    lua_getfield(L, config, "output");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, output undefined\n"); return -1; }
    int tab = lua_gettop(L);
    lua_getfield(L, tab, "dir");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, output.dir undefined\n"); return -1; }
    prefs->output.dir = strdup(lua_tostring(L, -1)); lua_pop(L, 1);
    lua_getfield(L, tab, "apsi");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, output.apsi undefined\n"); return -1; }
    prefs->output.apsi = strdup(lua_tostring(L, -1)); lua_pop(L, 1);
    lua_getfield(L, tab, "pot");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, output.pot undefined\n"); return -1; }
    prefs->output.pot = strdup(lua_tostring(L, -1)); lua_pop(L, 1);
    lua_getfield(L, tab, "corr");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, output.corr undefined\n"); return -1; }
    prefs->output.corr = strdup(lua_tostring(L, -1)); lua_pop(L, 1);
    lua_getfield(L, tab, "dftcorr");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, output.dftcorr undefined\n"); return -1; }
    prefs->output.dftcorr = strdup(lua_tostring(L, -1)); lua_pop(L, 1);
    lua_getfield(L, tab, "spectrum");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, output.spectrum undefined\n"); return -1; }
    prefs->output.spectrum = strdup(lua_tostring(L, -1)); lua_pop(L, 1);

    lua_getfield(L, tab, "numen");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, output.numen undefined\n"); return -1; }
    prefs->output.numen = strdup(lua_tostring(L, -1)); lua_pop(L, 1);
    lua_getfield(L, tab, "splen");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, output.splen undefined\n"); return -1; }
    prefs->output.splen = strdup(lua_tostring(L, -1)); lua_pop(L, 1);
    lua_getfield(L, tab, "aken");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, output.aken undefined\n"); return -1; }
    prefs->output.aken = strdup(lua_tostring(L, -1)); lua_pop(L, 1);


    lua_pop(L, 1);
  }
  printf("  output:\n");
  printf("    dir:         %s\n", prefs->output.dir);
  printf("    apsi:        %s/%s\n", prefs->output.dir, prefs->output.apsi);
  printf("    pot:         %s/%s\n", prefs->output.dir, prefs->output.pot);
  printf("    corr:        %s/%s\n", prefs->output.dir, prefs->output.corr);
  printf("    dftcorr:     %s/%s\n", prefs->output.dir, prefs->output.dftcorr);
  printf("    spectrum:    %s/%s\n", prefs->output.dir, prefs->output.spectrum);

  // number of steps
  printf("  derrived values are:\n");
  prefs->tsteps =  prefs->steps * prefs->runs;
  printf("    tsteps:      %d\n", prefs->tsteps);

  // position space delta
  prefs->dx = prefs->xpos->data[1] - prefs->xpos->data[0];
  printf("    dx:          %g\n", prefs->dx);

  // momentum space delta
  prefs->dk = 2 * M_PI / (prefs->range.max - prefs->range.min);
  printf("    dk:          %g\n", prefs->dk);

  // energy delty
  prefs->dE = 2 * M_PI / (prefs->dt * prefs->tsteps);
  printf("    dE:          %g\n", prefs->dE);

  double maxE = prefs->dE * prefs->runs / 2;
  printf("  maximal energy detectable is %g\n", maxE);

  // get the energy range for plot, this effectively just passed to gnuplot
  {
    lua_getfield(L, config, "enrgrange");
    if (!lua_isnil(L, -1)) {
      int range = lua_gettop(L);
      lua_rawgeti(L, range, 1);
      if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, enrgrange.min undefined\n"); return -1; }
      prefs->enrgrange.min = lua_tonumber(L, -1); lua_pop(L, 1);
      lua_rawgeti(L, range, 2);
      if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, enrgrange.max undefined\n"); return -1; }
      prefs->enrgrange.max = lua_tonumber(L, -1); lua_pop(L, 1);
      lua_rawgeti(L, range, 3);
      if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, enrgrange.win undefined\n"); return -1; }
      prefs->enrgrange.win = lua_tonumber(L, -1); lua_pop(L, 1);
      lua_rawgeti(L, range, 4);
      if (!lua_isnil(L, -1)) {
        prefs->enrgrange.sel = lua_tonumber(L, -1);
      } else {
        prefs->enrgrange.sel = 1.9;
      }
      lua_pop(L, 1);
    } else {
      prefs->enrgrange.min = 0;
      prefs->enrgrange.max = maxE;
    }
    lua_pop(L, 1);
  }
  printf("  enrgrange [%g;%g]\n", prefs->enrgrange.min, prefs->enrgrange.max);
  printf("  peak search window is %g\n", prefs->enrgrange.win);
  printf("  peak selector is %g sdev\n", prefs->enrgrange.sel);

  printf("\n");

  // we are done, release config table
  luaL_unref(L, LUA_REGISTRYINDEX, prefs->config);
  prefs->config = LUA_NOREF;

  return 0;
}

/*! \memberof preferences
 * this function runs the simulation according to preferences \a prefs
 */
int start_simulation(preferences_t * prefs)
{
  printf("starting simulation: \n");

  // create a new split operator usind values from prefs
  splitop_t * sop = splitop_new(prefs);

  // normalize the wavefunction
  splitop_prepare(sop);

  // save the normalized wavefunction for later plot
  fftw_complex * psi = sop->psi;
  carray_t * spsi = prefs->psi;
  for (int k = 0; k < spsi->length; k++) {
    spsi->data[k] = psi[k];
  }

  // save wavefunction so we have a reference
  splitop_save(sop);
  fftw_complex * apsi = sop->apsi;

  printf("  split operator created\n");

  // allocate new results struct
  results_t * res = results_new(); prefs->results = res;

  // cache some values (hopefully gcc will place them in registers)
  int bins = prefs->bins;
  int runs = prefs->runs;
  int steps= prefs->steps;
  int length = runs + 1;

  // allocate an array for the corraltion function and spectrum and the fftw plan
  carray_t * c = carray_new_sized(0, length);
  carray_t * co = carray_new_sized(0, length);
  carray_t * ck = carray_new_sized(0, length);
  fftw_plan p = fftw_plan_dft_1d(length, c->data, ck->data, FFTW_FORWARD, FFTW_ESTIMATE);

#ifdef USE_CAIRO
  // render video if enabled
  if (prefs->vstep > 0 && prefs->vframes > 0) {
    struct stat sb;
    int len = strlen(prefs->output.dir) + 100;
    int width = 1000, height = 400;
    char path[len];

    if (stat(prefs->output.dir, &sb)) {
      if (mkdir(prefs->output.dir, 0777)) {
        assert(0);
      }
    }

    snprintf(path, len, "%s/video", prefs->output.dir);
    if (stat(path, &sb)) {
      if (mkdir(path, 0777)) {
        assert(0);
      }
    }

    cairo_surface_t * surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
    cairo_t *cr = cairo_create(surface);
    for (int k = 0; k <= prefs->vframes ; k++) {
      fprintf(stderr, "  video (%g), frame %u\r", cvect_normsq(bins, psi), k);
      splitop_draw(sop, cr, (cairo_rectangle_t){0, 0, width, height}, psi);
      snprintf(path, len, "%s/video/image%05u.png", prefs->output.dir, k);
      cairo_surface_write_to_png(surface, path);
      splitop_run(sop, prefs->vstep);
    }
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    printf("  done...                 \n");
  }
  // restore initial wavefunction befora proceeding
  splitop_restore(sop);
#endif

  // this should always be 1, so just a safety check
  c = carray_append(c, cvect_skp(bins, psi, apsi));

  printf("  starting run\n");

  /******** start simulation **********/
  for (int k = 0; k < runs; k++) {
    splitop_run(sop, steps);
    if (k % 100 == 0) {
      printf("    %d/%d\r", k, runs);
      fflush(stdout);
    }
    c = carray_append(c, cvect_skp(bins, psi, apsi));
  }
  /********* simulation done **********/

  printf("  done...                 \n");

  printf("  hanning before dft\n");
  // von Hann Fenster:
  double hannfkt = 2 * M_PI/(length);
  for (int k = 0; k < length; k++) {
    complex double z = c->data[k];
    co->data[k] = z; // save original function
    c->data[k] = z * .5 * (1 - cos(hannfkt * k));
  }

  fftw_execute_dft(p, c->data, ck->data);
  printf("  calculated dft\n");

  // normalize spectrum, as fftw will not do this
  double nor = 1/sqrt(length);
  for (int k = 0; k < length; k++) { ck->data[k] *= nor; }

  res->c = c;
  res->co = co;
  res->ck = ck;
  ck->length = c->length;
  co->length = c->length;

  // dispose of now unneded ressources
  fftw_destroy_plan(p);
  splitop_free(sop);

  printf("  cleaned up\n");

  return 0;
}

void dump_results(preferences_t * prefs, FILE * fp)
{
  results_t * res = prefs->results;
#define WRITE_FIELD(_f) assert(fwrite(&res->_f, sizeof(res->_f), 1, fp) == 1)
  WRITE_FIELD(co->length);
  for (int k = 0; k < res->co->length; k++) {
    WRITE_FIELD(co->data[k]);
  }
  WRITE_FIELD(c->length);
  for (int k = 0; k < res->c->length; k++) {
    WRITE_FIELD(c->data[k]);
  }
  WRITE_FIELD(ck->length);
  for (int k = 0; k < res->ck->length; k++) {
    WRITE_FIELD(ck->data[k]);
  }
#undef WRITE_FIELD
}

void undump_results(preferences_t * prefs, FILE * fp)
{
#define READ_ARRAY(_f, _new) \
  do { \
    int len; \
    assert(fread(&len, sizeof(len), 1, fp) == 1); \
    res->_f = _new(len); assert(res->_f); \
    size_t sz = sizeof(res->_f->data[0]) * len; \
    assert(fread(res->_f->data, 1, sz, fp) == sz); \
  } while(0)
  results_t * res = malloc(sizeof(results_t)); assert(res);
  READ_ARRAY(co, carray_new);
  READ_ARRAY(c, carray_new);
  READ_ARRAY(ck, carray_new);
#undef READ_ARRAY
  prefs->results = res;
}

/*! \memberof preferences
 * this function just evaluates the results \a prefs
 */
int eval_results(preferences_t * prefs)
{
  int steps = prefs->steps;
  double dt = prefs->dt;
  double dE = prefs->dE;
  results_t * res = prefs->results;
  carray_t * c = res->co; // the original function without hanning
  carray_t * ck = res->ck;

  array_t * xpos = prefs->xpos;
  carray_t * apsi = prefs->psi;
  array_t * pot = prefs->potential;

  struct stat sb;

  if (stat(prefs->output.dir, &sb)) {
    if (mkdir(prefs->output.dir, 0777)) {
      assert(0);
    }
  }

  // dump initial wafevunction
  {
    int len = strlen(prefs->output.dir) + strlen(prefs->output.apsi) + 10;
    char path[len];
    snprintf(path, len, "%s/%s", prefs->output.dir, prefs->output.apsi);
    FILE * fp = fopen(path, "w"); assert(fp);
    for (int k = 0; k < apsi->length; k++) {
      complex double z = apsi->data[k];
      fprintf(fp, "%.17e %.17e %.17e %.17e\n", xpos->data[k], creal(z), cimag(z), cabs(z));
    }
    fclose(fp);
  }

  // dump potnetial used during simulation
  {
    int len = strlen(prefs->output.dir) + strlen(prefs->output.pot) + 10;
    char path[len];
    snprintf(path, len, "%s/%s", prefs->output.dir, prefs->output.pot);
    FILE * fp = fopen(path, "w"); assert(fp);
    for (int k = 0; k < apsi->length; k++) {
      fprintf(fp, "%.17e %.17e\n", xpos->data[k], pot->data[k]);
    }
    fclose(fp);
  }

  // dump the correlation function
  {
    int len = strlen(prefs->output.dir) + strlen(prefs->output.corr) + 10;
    char path[len];
    snprintf(path, len, "%s/%s", prefs->output.dir, prefs->output.corr);
    FILE * fp = fopen(path, "w"); assert(fp);
    for (int k = 0; k < c->length; k++) {
      complex double z = c->data[k];
      fprintf(fp, "%.17e %.17e %.17e %.17e\n", k * steps * dt, creal(z), cimag(z), cabs(z));
    }
    fclose(fp);
  }

  // dump the DTF of corr
  {
    int len = strlen(prefs->output.dir) + strlen(prefs->output.dftcorr) + 10;
    char path[len];
    snprintf(path, len, "%s/%s", prefs->output.dir, prefs->output.dftcorr);
    FILE * fp = fopen(path, "w"); assert(fp);
    /*for (int k = 0; k < ck->length; k++) {
      complex double z = ck->data[k];
      fprintf(fp, "%.17e %.17e %.17e %.17e\n", k * dE, creal(z), cimag(z), cabs(z));
    }*/
    for (int k = ck->length/2; k < ck->length; k++) {
      double x = (k - ck->length) * dE;
      complex double z = ck->data[k];
      fprintf(fp, "%.17e %.17e %.17e %.17e\n", x, creal(z), cimag(z), cabs(z));
    }
    for (int k = 0; k < ck->length/2; k++) {
      double x = k * dE;
      complex double z = ck->data[k];
      fprintf(fp, "%.17e %.17e %.17e %.17e\n", x, creal(z), cimag(z), cabs(z));
    }
    fclose(fp);
  }

  // dump spectra
  {
    int len = strlen(prefs->output.dir) + strlen(prefs->output.spectrum) + 10;
    char path[len];
    snprintf(path, len, "%s/%s", prefs->output.dir, prefs->output.spectrum);
    FILE * fp = fopen(path, "w"); assert(fp);

    /*array_t * data = array_new(ck->length);
    for (int k = 0; k < ck->length; k++) {
      data->data[k] = cabs(ck->data[k]);
    }*/

    int length = ck->length, o = 0;
    int * index = malloc(sizeof(int) * length); assert(index);

    /* create a map to place negative energies in the right place */
    for (int k = ck->length/2; k < ck->length; k++) { index[o++] = k; }
    for (int k = 0; k < ck->length/2; k++) { index[o++] = k; }

    array_t * data = carray_abs(ck, index);
    array_t * logdat = array_map(data, log);
    array_t * peaks = peaks_find(prefs->dE, logdat, prefs->enrgrange.win, prefs->enrgrange.sel);

    int odd = (length%2);

    {
      int m = prefs->enrgrange.min / dE + .5 + length/2.;
      int M = prefs->enrgrange.max / dE + .5 + length/2.;

      assert(M < length && m >= 0 && m < M);

      array_t * s = array_new(M - m + 1);
      array_t * f = array_new(s->length);
      for (int k = 0; k < s->length; k++) {
        s->data[k] = (m + k - ck->length/2 - odd) * dE;
        f->data[k] = logdat->data[m + k];
      }

      /* spline search */
      array_t * a = array_cspline_prepare(f, prefs->dE);
      snprintf(path, len, "%s/%s", prefs->output.dir, prefs->output.splen);
      array_t * d = array_cspline_zroots(s, f, a, prefs->dE, 1);
      array_dump_to_file(path, " ", 1, d);
      free(a);
      free(d);

      /* akima search */
      akima_t * ak = akima_new(s, f);
      array_t * ao = akima_zroots(ak, 1);
      snprintf(path, len, "%s/%s", prefs->output.dir, prefs->output.aken);
      array_dump_to_file(path, " ", 1, ao);
      akima_free(ak);
      free(ao);
    }

    printf("-------- %d\n", data->length);
    for (int k = 0; k < peaks->length; k++) {
      int i = peaks->data[k];
      if (i > 0 && i < data->length) {
        double z = (i - data->length/2 - odd) * dE;
        if (z >= prefs->enrgrange.min && z <= prefs->enrgrange.max) {
          fprintf(fp, "%.17e\n", z);
        }
      }
    }

    fclose(fp);
    free(peaks);
    free(data);
    free(logdat);
    free(index);
  }

  {
    array_t * numen = numerov_energies(prefs);
    int len = strlen(prefs->output.dir) + strlen(prefs->output.numen) + 10;
    char path[len];
    snprintf(path, len, "%s/%s", prefs->output.dir, prefs->output.numen);
    array_dump_to_file(path, " ", 1, numen);
    free(numen);
  }

  // dump variables for gnuplot
  {
    int len = strlen(prefs->output.dir) + strlen("stats") + 10;
    char path[len];
    snprintf(path, len, "%s/%s", prefs->output.dir, "stats");
    FILE * fp = fopen(path, "w"); assert(fp);

    fprintf(fp, "%d:", prefs->bins);
    fprintf(fp, "%.17e:", prefs->dt);
    fprintf(fp, "%.17e:", prefs->range.min);
    fprintf(fp, "%.17e:", prefs->range.max);
    fprintf(fp, "%d:", prefs->steps);
    fprintf(fp, "%d:", prefs->runs);

    {
      double min = 0, max=0;
      for (int k = 0; k < pot->length; k++) {
        double z = pot->data[k];
        if (z < min) {
          min = z;
        }
        if (z > max) {
          max = z;
        }
      }
      fprintf(fp, "%.17e:", min);
      fprintf(fp, "%.17e:", max);
    }

    {
      double min = 0, max=0;
      for (int k = 0; k < apsi->length; k++) {
        double z = cabs(apsi->data[k]);
        if (z < min) {
          min = z;
        }
        if (z > max) {
          max = z;
        }
      }
      fprintf(fp, "%.17e:", min);
      fprintf(fp, "%.17e:", max);
    }

    fprintf(fp, "%d:", prefs->tsteps);
    fprintf(fp, "%.17e:", prefs->dx);
    fprintf(fp, "%.17e:", prefs->dk);
    fprintf(fp, "%.17e:", prefs->dE);

    fprintf(fp, "%.17e:", prefs->enrgrange.min);
    fprintf(fp, "%.17e:", prefs->enrgrange.max);

    fclose(fp);
  }

  return 0;
}