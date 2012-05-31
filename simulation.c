#include <lauxlib.h>
#include <lualib.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "simulation.h"

results_t * results_new() {
  results_t * res = malloc(sizeof(results_t));
  memset(res, 0, sizeof(results_t ));
  return res;
}

void results_free(results_t * res) {
  free(res->c);
  free(res->ck);
  free(res);
}

preferences_t * preferences_new() {
  preferences_t * prefs = malloc(sizeof(preferences_t));
  memset(prefs, 0, sizeof(preferences_t));
  return prefs;
}

void preferences_free(preferences_t * prefs) {
  free(prefs->xpos);
  free(prefs->potential);
  free(prefs->psi);

  free(prefs->output.dir);
  free(prefs->output.apsi);
  free(prefs->output.pot);
  free(prefs->output.corr);
  free(prefs->output.dftcorr);

  if (prefs->results) {
    results_free(prefs->results);
  }

  free(prefs);
}

int preferences_read(lua_State * L, preferences_t * prefs) {
  lua_settop(L, 0);
  lua_rawgeti(L, LUA_REGISTRYINDEX, prefs->config);
  int config = lua_gettop(L);

  printf("configuration: \n");

  lua_getfield(L, config, "bins");
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, bins undefined\n"); return -1; }
  prefs->bins = lua_tointeger(L, -1); lua_pop(L, 1);
  printf("  bins:    %d\n", prefs->bins);

  lua_getfield(L, config, "dt");
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, dt undefined\n"); return -1; }
  prefs->dt = lua_tonumber(L, -1); lua_pop(L, 1);
  printf("  dt:      %g\n", prefs->dt);

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

  lua_getfield(L, config, "steps");
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, steps undefined\n"); return -1; }
  prefs->steps = lua_tointeger(L, -1); lua_pop(L, 1);
  printf("  steps:   %d\n", prefs->steps);

  lua_getfield(L, config, "vstep");
  if (!lua_isnil(L, -1)) {
    prefs->vstep = lua_tointeger(L, -1);
  } else {
    prefs->vstep = -1;
  }
  lua_pop(L, 1);
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

  lua_getfield(L, config, "runs");
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, runs undefined\n"); return -1; }
  prefs->runs = lua_tointeger(L, -1); lua_pop(L, 1);
  printf("  runs:    %d\n", prefs->runs);

  array_t * xpos = array_equipart(prefs->range.min, prefs->range.max, prefs->bins);
  prefs->xpos = xpos;

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
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, output.corr undefined\n"); return -1; }
    prefs->output.dftcorr = strdup(lua_tostring(L, -1)); lua_pop(L, 1);
    lua_pop(L, 1);
  }
  printf("  output:\n");
  printf("    dir:         %s\n", prefs->output.dir);
  printf("    apsi:        %s/%s\n", prefs->output.dir, prefs->output.apsi);
  printf("    pot:         %s/%s\n", prefs->output.dir, prefs->output.pot);
  printf("    corr:        %s/%s\n", prefs->output.dir, prefs->output.corr);
  printf("    dftcorr:     %s/%s\n", prefs->output.dir, prefs->output.dftcorr);

  printf("  derrived values are:\n");
  prefs->tsteps =  prefs->steps * prefs->runs;
  printf("    tsteps:      %d\n", prefs->tsteps);

  prefs->dx = prefs->xpos->data[1] - prefs->xpos->data[0];
  printf("    dx:          %g\n", prefs->dx);

  prefs->dk = 2 * M_PI / (prefs->range.max - prefs->range.min);
  printf("    dk:          %g\n", prefs->dk);

  prefs->dE = 2 * M_PI / (prefs->dt * prefs->tsteps);
  printf("    dE:          %g\n", prefs->dE);

  printf("\n");

  luaL_unref(L, LUA_REGISTRYINDEX, prefs->config);
  prefs->config = LUA_NOREF;

  return 0;
}