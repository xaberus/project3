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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <complex.h>
#include <fftw3.h>
#include <pthread.h>

#include "simulation.h"

#include <assert.h>

int main(int argc, char * argv[argc])
{
  fftw_init_threads();
  fftw_plan_with_nthreads(4);

  if (argc != 3) {
    fprintf(stderr, "usage: %s config.lua result.dat\n", argv[0]);
    return -1;
  }

  lua_State * L = luaL_newstate();
  luaL_openlibs(L);

  preferences_t * prefs = preferences_new();

  if (luaL_dofile(L, argv[1])) {
    fprintf(stderr, "could not load '%s' : %s\n", argv[1], lua_tostring(L, 1));
  } else {
    // get config table
    lua_getfield(L, LUA_GLOBALSINDEX, "config");
    if (lua_isnil(L, -1)) {
      fprintf(stderr, "table config undefined\n");
    } else {
      // ref config table, so we can access it in preferences_read()
      prefs->config = luaL_ref(L, LUA_REGISTRYINDEX);
      if (!preferences_read(L, prefs)) {
        FILE * fp = fopen(argv[2], "rb"); assert(fp);
        undump_results(prefs, fp);
        fclose(fp);

        if (!prefs->cmp) {
          eval_results(prefs);
        } else {
          int len = strlen(prefs->output.dir) + 100;
          char path[len];

          carray_t * co = prefs->results->co;

          int num = prefs->runs / 10000;

          for (int length = 10000; length <= co->length; length += co->length / num) {
            carray_t * c = carray_pcopy(length, co->data);

            double hannfkt = 2 * M_PI/(length);
            for (int k = 0; k < length; k++) {
              complex double z = c->data[k];
              co->data[k] = z; // save original function
              c->data[k] = z * .5 * (1 - cos(hannfkt * k));
            }
            carray_t * ck = carray_new(c->length);
            fftw_plan p = fftw_plan_dft_1d(length, c->data, ck->data, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute_dft(p, c->data, ck->data);

            double nor = 1/sqrt(length);
            for (int k = 0; k < length; k++) { ck->data[k] *= nor; }

            double dE = 2 * M_PI / (prefs->dt * prefs->steps * length);

            spectra_t * spec = spectra_search(ck, dE,
              prefs->enrgrange.min, prefs->enrgrange.max,
              prefs->enrgrange.win, prefs->enrgrange.sel);

            snprintf(path, len, "%s/o-%d.dat", prefs->output.dir, length);
            array_dump_to_file(path, " ", 2, spec->s, spec->f);
            snprintf(path, len, "%s/d-%d.dat", prefs->output.dir, length);
            array_dump_to_file(path, " ", 1, spec->diren);
            snprintf(path, len, "%s/s-%d.dat", prefs->output.dir, length);
            array_dump_to_file(path, " ", 1, spec->splen);
            snprintf(path, len, "%s/a-%d.dat", prefs->output.dir, length);
            array_dump_to_file(path, " ", 1, spec->aken);

            spectra_free(spec);

            fftw_destroy_plan(p);
            free(c);
            free(ck);
          }
        }
      }
    }
  }

  preferences_free(prefs);

  lua_close(L);

  fftw_cleanup();
  fftw_cleanup_threads();
  pthread_exit(NULL);
}