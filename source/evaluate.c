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

          for (int length = 10000; length <= co->length; length += co->length / 60) {
            int odd = (length%2);
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

            int o = 0;
            int * index = malloc(sizeof(int) * length); assert(index);
            /* create a map to place negative energies in the right place */
            for (int k = ck->length/2; k < ck->length; k++) { index[o++] = k; }
            for (int k = 0; k < ck->length/2; k++) { index[o++] = k; }
            array_t * data = carray_abs(ck, index);

            double dE = 2 * M_PI / (prefs->dt * prefs->steps * length);

            int m = prefs->enrgrange.min / dE + .5 + length/2;
            int M = prefs->enrgrange.max / dE + .5 + length/2;

            assert(m >= 0 && m < M && M < co->length);

            array_t * s = array_new(M - m + 1);
            array_t * f = array_new(s->length);
            for (int k = 0; k < s->length; k++) {
              s->data[k] = (m + k - ck->length/2 - odd + 1) * dE;
              f->data[k] = log(data->data[m + k]);
            }

            snprintf(path, len, "%s/o-%d.dat", prefs->output.dir, length);
            array_dump_to_file(path, " ", 2, s, f);

            {
              snprintf(path, len, "%s/d-%d.dat", prefs->output.dir, length);
              array_t * peaks = direct_search(dE, s, f,
                prefs->enrgrange.win, prefs->enrgrange.sel);
              array_dump_to_file(path, " ", 1, peaks);
              free(peaks);
            }

            {
              snprintf(path, len, "%s/s-%d.dat", prefs->output.dir, length);
              array_t * peaks = spline_search(s, f);
              array_dump_to_file(path, " ", 1, peaks);
              free(peaks);
            }

            {
              snprintf(path, len, "%s/a-%d.dat", prefs->output.dir, length);
              array_t * peaks = akima_search(s, f);
              array_dump_to_file(path, " ", 1, peaks);
              free(peaks);
            }


            fftw_destroy_plan(p);
            free(c);
            free(index);
            free(data);
            free(s);
            free(f);
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