#include <stdio.h>
#include <stdlib.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <complex.h>
#include <fftw3.h>
#include <pthread.h>

#include "simulation.h"

/*!
 this is a straight forward main function, that does the following
  - load lua config file
  - get configuration table from config
  - contruct a preferences class from them
  - start simulation according to the preferences
  - dump the result as precified in preferences
 */
int main(int argc, char * argv[argc])
{
  fftw_init_threads();
  fftw_plan_with_nthreads(4);

  if (argc < 2) {
    fprintf(stderr, "usage: %s config.lua\n", argv[0]);
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
        if (!start_simulation(prefs)) {
          dump_results(prefs);
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