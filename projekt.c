#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include <assert.h>
#include <stdint.h>

//#include <cairo.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include "splitop.h"
#include "cvect.h"
#include "array.h"
#include "carray.h"

/*************************************************************************************************/

#if 0
void splitop_draw(splitop_t * w, cairo_t * cr, cairo_rectangle_t rect, fftw_complex * psi) {
  cairo_save(cr);
  cairo_rectangle(cr, rect.x, rect.y, rect.width, rect.height);
  cairo_clip(cr);

  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_set_line_width(cr, 2);
  cairo_paint(cr);

  double y0 = rect.height/2;

  //cairo_set_line_width(cr, 10);
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_move_to(cr, rect.x, rect.y + y0);
  cairo_line_to(cr, rect.x + rect.width, rect.y + y0);
  cairo_stroke(cr);

  int bins = w->bins;

  double xscale = rect.width/bins, yscale, max;

  max = 0;
  for (int n = 0; n < bins; n++) {
    double a = fabs(w->V[n]); if (max < a) { max = a; }
  }
  yscale = y0/max;

  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_move_to(cr, rect.x, y0);
  for (int n = 0; n < bins; n++) { cairo_line_to(cr, rect.x + n * xscale, rect.y + y0 - w->V[n] * yscale); }
  cairo_stroke(cr);

  fftw_complex * apsi = w->apsi;
  max = 0;
  for (int n = 0; n < bins; n++) {
    double a = cabs(apsi[n]); if (max < a) { max = a; }
  }
  yscale = y0/(5*max/4);

  cairo_set_line_width(cr, 1);
  cairo_set_source_rgb(cr, 0, 0, 1);
  cairo_move_to(cr, rect.x, y0);
  for (int n = 0; n < bins; n++) { cairo_line_to(cr, rect.x + n * xscale, rect.y + y0 - cabs(apsi[n]) * yscale); }
  cairo_stroke(cr);

  cairo_set_line_width(cr, 2);
  cairo_set_source_rgb(cr, 1, 0, 0);
  cairo_move_to(cr, rect.x, y0);
  for (int n = 0; n < bins; n++) { cairo_line_to(cr, rect.x + n * xscale, rect.y + y0 - cabs(psi[n]) * yscale); }
  cairo_stroke(cr);

  fftw_complex * psik = fftw_alloc_complex(bins); assert(psik);
  fftw_execute_dft(w->fwd, psi, psik);

  max = 0;
  for (int n = 0; n < bins; n++) {
    double a = cabs(psik[n]); if (max < a) { max = a; }
  }
  yscale = y0/(5*max/4);

  cairo_set_line_width(cr, 1);
  cairo_set_source_rgb(cr, 1, .5, 0);
  /*cairo_move_to(cr, rect.x, y0);
  for (int n = 0; n < bins; n++) {
    int l = (n+bins/2)%w->bins;
    cairo_line_to(cr, rect.x + (n-bins/3)*4 * xscale, rect.y + 2*y0 - cabs(psik[l]) * yscale);
  }*/
  cairo_move_to(cr, rect.x, y0);
  int dron = 1;
  double reg = 1.0 / bins;
  for (int k = bins/2; k < bins; k++) {
    double x = rect.x + rect.width/2 + (k-bins) * reg * rect.width;
    double y = rect.y + 2*y0 - cabs(psik[k]) * yscale;
    if (dron) { cairo_move_to(cr, x, y); dron = 0; } else { cairo_line_to(cr, x, y); }
  }
  for (int k = 0; k < bins/2; k++) {
    double x = rect.x + rect.width/2 + k * reg * rect.width;
    double y = rect.y + 2*y0 - cabs(psik[k]) * yscale;
    cairo_line_to(cr, x, y);
  }
  cairo_stroke(cr);

  fftw_free(psik);

  cairo_restore(cr);
}
#endif

/*************************************************************************************************/

#define HANN

//#define HARMOSZ
#define BOX

#ifdef HARMOSZ
double omega = 40;
fftw_complex zeropot(double x) {
  return omega*x*x;
}
fftw_complex fn(double x) {
  //double a = .2;
  double a = .5;
  double aa = a * a;
  double x0 = 0;
  double k0 = 10; // <-- good one
  double xx = (x-x0) * (x-x0);
  return  exp(-xx/(aa))*cexp(I*k0*(x));
}
#endif
#ifdef BOX
double pota = 5;
fftw_complex zeropot(double x) {
  double r = fabs(x);
  return r < pota ? 0 : 1000;
}
fftw_complex fn(double x) {
  return  exp(-x*x/(2)) * cexp(I * 4 * x);
}


/*double pota = 4;
fftw_complex zeropot(double x) {
  double r = fabs(x);
  return r < pota ? 0 : 10000;
}*/
/*fftw_complex fn(double x) {
  //double a = .2;
  double a = .1;
  double aa = a * a;
  double x0 = 0;
  //double k0 = 10;
  double k0 = 0;
  double xx = (x-x0) * (x-x0);
  return  exp(-xx/(aa))*cexp(I*k0*(x));
}*/
/*fftw_complex fn(double x) {
  return  //exp(-x*x/(.1))
    + I * exp(-x*x/(.2*.2)) * sin(x)
    ;
}*/
/*fftw_complex fn(double x) {
  //double a = .2;
  double a = .5;
  double aa = a * a;
  double x0 = 0;
  double k0 = 10; // <-- good one
  double xx = (x-x0) * (x-x0);
  return  exp(-xx/(aa))*cexp(I*k0*(x));
}*/
#endif

#if 0
double pota = 0.9;
//double pota = 3;
double potb = 1;
double potc = 30;

/*double pota = 2;
double potb = 1;*/
fftw_complex zeropot(double x) {
  //double r = fabs(x); return r < 6 ? 0 : 1000000000;
  //return 30*x*x;
  return 10*x*x; // <-- good one
  //return 10*x*x*x*x;
  //return 10*x*x;
  //return 40*cosh(x);

  /*double b = potb * potb;
  double l = x + pota;
  double r = x - pota;
  return -potc*exp(-l*l/b)-potc*exp(-r*r/b);*/

  //return 40*(-1/sqrt(1+potb*pow(x-pota,2))-1/sqrt(1+potb*pow(x+pota,2)));

  /*double a = pota;
  double b = 2;
  double r = fabs(x);
  if (r > fabs(a/2-b/2) && r < a/2+b/2) return -100;
  if (r < fabs(a/2-b/2)) return -30;
  //if (r > 25) return 1000;
  return 0;*/
}

fftw_complex fn(double x) {
  //double a = .2;
  double a = .5;
  double aa = a * a;

  double x0 = 0;

  //double x0 = -pota;

  //double x0 = -1.941063416246160;

  double k0 = 10; // <-- good one
  //double k0 = 1;
  double xx = (x-x0) * (x-x0);
  return  exp(-xx/(aa))*cexp(I*k0*(x));
}
#endif


#if 0
int nomain(/*int argc, char *argv[]*/) {
  fftw_init_threads();
  fftw_plan_with_nthreads(2);

  //splitop_t * sop = splitop_new(4096, .0001, -25, 25, zeropot);
#ifdef HARMOSZ
  splitop_t * sop = splitop_new(2*4096, .0001, -25, 25, zeropot); // <-- good one
#endif
#ifdef BOX
  splitop_t * sop = splitop_new(2*4096, .0001, -15, 15, zeropot); // <-- good one
  //splitop_t * sop = splitop_new(4096, .0001, -25, 25, zeropot); // <-- good one
#endif
  //splitop_t * sop = splitop_new(4096, .0001, -25, 25, zeropot);

  fftw_complex * psi = splitop_prepare(sop, fn);

  splitop_save(sop, psi);

  //splitop_run(sop, 1, psi);

#ifdef P4
  int width = 1000, height = 400;
  char buf[256];
  cairo_surface_t * surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
  cairo_t *cr = cairo_create(surface);

  splitop_draw(sop, cr, (cairo_rectangle_t){0, 0, width, height}, psi);
  snprintf(buf, sizeof(buf), "image%05u.png", 0);
  cairo_surface_write_to_png(surface, buf);

  for (int k = 1; k < 220; k++) {
    fprintf(stderr, "##### %g, %u\n", cvect_normsq(sop->bins, psi), k);
    //fprintf(stderr, "##### %u\r", k);
    splitop_run(sop, 150, psi);
    splitop_draw(sop, cr, (cairo_rectangle_t){0, 0, width, height}, psi);
    snprintf(buf, sizeof(buf), "image%05u.png", k);
    cairo_surface_write_to_png(surface, buf);
  }

  fprintf(stderr, "done---------\n");

  cairo_destroy(cr);
  cairo_surface_destroy(surface);
#endif

#ifdef P3
  //printf("### %g, %g\n", cvect_skp(sop->bins, psi, psi));

  int hli = 50000;
  int lli = 10;
  int mli = hli * lli;

  carray_t * c = carray_new_sized(0, hli);

  c = carray_append(c, cvect_skp(sop->bins, psi, sop->apsi));

  for (int k = 1; k < hli; k++) {
    splitop_run(sop, lli, psi);
    if (k % 100 == 0) fprintf(stderr, "%u/%u\r", k, hli);
    c = carray_append(c, cvect_skp(sop->bins, psi, sop->apsi));
  }
  fprintf(stderr, "done---------\n");

  {
    FILE * fp = fopen("plotc.txt", "w"); assert(fp);
    for (int k = 0; k < c->length; k++) {
      fprintf(fp, "%.17e %.17e %.17e %.17e\n",
        k * lli * sop->dt, creal(c->data[k]), cimag(c->data[k]), cabs(c->data[k]));
    }
    fclose(fp);
  }

  int length = c->length;

  {
    fftw_complex * ck = fftw_alloc_complex(length); assert(ck);

#ifdef HANN
    // Hann Fenster:
    double hannfkt = 2 * M_PI/(length);
    for (int k = 0; k < length; k++) { c->data[k] *= .5 * (1 - cos(hannfkt * k)); }
#endif

    fftw_plan p = fftw_plan_dft_1d(length, c->data, ck, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    {
      double nor = 1/sqrt(length);
      for (int k = 0; k < length; k++) { ck[k] *= nor; }
      FILE * fp = fopen("plotck.txt", "w"); assert(fp);
      double dom = 2*M_PI/(sop->dt * lli * hli);
      for (int k = 0; k < c->length; k++) {
        fprintf(fp, "%.17e %.17e %.17e %.17e\n",
          k * dom, creal(ck[k]), cimag(ck[k]), cabs(ck[k]));
      }
      fclose(fp);
    }
    {
      FILE * fp = fopen("plotE.txt", "w"); assert(fp);
#ifdef HARMOSZ
      double A = sqrt(4 * omega);
      for (int k = 0; k < c->length; k++) {
        fprintf(fp, "%.17e\n", A * (0.5 + k));
      }
#endif
#ifdef BOX
      double L = 2 * pota;
      for (int k = 0; k < c->length; k++) {
        int m = k * 10;
        fprintf(fp, "%.17e\n", M_PI * M_PI * m * m / (L * L));
      }
#endif
      fclose(fp);
    }

    fftw_destroy_plan(p);
    fftw_free(ck);
  }

#endif

  fftw_free(psi);
  splitop_free(sop);
  fftw_cleanup_threads();
  fftw_cleanup();
  return 0;
}
#endif

typedef struct {
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

  int        tsteps;
  double     dx;
  double     dk;
  double     dE;
} preferences_t;

void prepare_simulation(lua_State * L, preferences_t * prefs) {
  lua_settop(L, 0);
  lua_rawgeti(L, LUA_REGISTRYINDEX, prefs->config);
  int config = lua_gettop(L);

  printf("configuration: \n");

  lua_getfield(L, config, "bins");
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, bins undefined\n"); return; }
  prefs->bins = lua_tointeger(L, -1); lua_pop(L, 1);
  printf("  bins:    %d\n", prefs->bins);

  lua_getfield(L, config, "dt");
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, dt undefined\n"); return; }
  prefs->dt = lua_tonumber(L, -1); lua_pop(L, 1);
  printf("  dt:      %g\n", prefs->dt);

  {
    lua_getfield(L, config, "range");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, range undefined\n"); return; }
    int range = lua_gettop(L);
    lua_rawgeti(L, range, 1);
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, range.min undefined\n"); return; }
    prefs->range.min = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_rawgeti(L, range, 2);
    prefs->range.max = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  printf("  range:   [%g;%g]\n", prefs->range.min, prefs->range.max);

  lua_getfield(L, config, "steps");
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, steps undefined\n"); return; }
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
  if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, runs undefined\n"); return; }
  prefs->runs = lua_tointeger(L, -1); lua_pop(L, 1);
  printf("  runs:    %d\n", prefs->runs);

  array_t * xpos = array_equipart(prefs->range.min, prefs->range.max, prefs->bins);
  prefs->xpos = xpos;

  {
    lua_getfield(L, config, "potential");
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, potential undefined\n"); return; }
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
    if (lua_isnil(L, -1)) { fprintf(stderr, "aborting, psi undefined\n"); return; }
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

  printf("  derrived values are:\n");
  prefs->tsteps =  prefs->steps * prefs->runs;
  printf("    tsteps:      %d\n", prefs->tsteps);

  prefs->dx = prefs->xpos->data[1] - prefs->xpos->data[0];
  printf("    dx:          %g\n", prefs->dx);

  prefs->dE = 2 * M_PI / (prefs->dt * prefs->tsteps);
  printf("    dE:          %g\n", prefs->dE);

}

int main(int argc, char * argv[argc]) {

  if (argc < 2) {
    fprintf(stderr, "usage: %s config.lua\n", argv[0]);
    return -1;
  }

  preferences_t prefs = {0};

  lua_State * L = luaL_newstate();
  luaL_openlibs(L);

  if (luaL_dofile(L, argv[1])) {
    fprintf(stderr, "could not load '%s' : %s\n", argv[1], lua_tostring(L, 1));
  } else {
    lua_getfield(L, LUA_GLOBALSINDEX, "config");
    if (lua_isnil(L, -1)) {
      fprintf(stderr, "table config undefined\n");
    } else {
      prefs.config = luaL_ref(L, LUA_REGISTRYINDEX);
      prepare_simulation(L, &prefs);
    }
  }

  lua_close(L);

  return 0;
}