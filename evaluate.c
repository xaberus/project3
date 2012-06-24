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
    fprintf(stderr, "usage: %s config.ser results.dat\n", argv[0]);
    return -1;
  }

  preferences_t * prefs = preferences_new();

  FILE * fp = fopen(argv[1], "rb"); assert(fp);
  preferences_deserialize(prefs, fp);
  fclose(fp);

  fp = fopen(argv[2], "rb"); assert(fp);
  undump_results(prefs, fp);
  fclose(fp);

  eval_results(prefs);

  preferences_free(prefs);

  fftw_cleanup();
  fftw_cleanup_threads();
  pthread_exit(NULL);
}