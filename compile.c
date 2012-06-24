#include <stdio.h>
#include <stdlib.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <complex.h>
#include <fftw3.h>
#include <pthread.h>

#include <assert.h>

#include "simulation.h"

int main(int argc, char * argv[argc])
{
  fftw_init_threads();
  fftw_plan_with_nthreads(4);

  if (argc != 3) {
    fprintf(stderr, "usage: %s config.lua config.ser\n", argv[0]);
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
        FILE * fp = fopen(argv[2], "wb"); assert(fp);
        preferences_serialize(prefs, fp);
        fclose(fp);
      }
    }
  }

  preferences_free(prefs);

  lua_close(L);

  fftw_cleanup();
  fftw_cleanup_threads();
  pthread_exit(NULL);
}