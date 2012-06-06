#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "peaks.h"

typedef struct {
  int init;
  int cols;
} state_t;

void parse_line(state_t * s, size_t len, char buf[len])
{
  char * end = buf + len;
  if (!s->init) {
    s->cols = 0;
    for (char * p = buf; p <end; p++) {
      if (*p == ';') {
        s->cols++;
      }
    }
    if (len && buf[len-1] != ';') {
      s->cols++;
    }
    fprintf(stdout, "\\begin{tabular}{|");
    for (int k = 0; k < s->cols; k++) {
      fprintf(stdout, "c|");
    }
    fprintf(stdout, "}\\hline\n");
    s->init = 1;
  }

  //fprintf(stdout, "%.*s\n", (int) len, buf);
  fprintf(stdout, "  ");
  int k = 0;
  for (char * tmp = buf, * tok, * sav; k < s->cols; tmp = NULL, k++) {
    tok = strtok_r(tmp, ";", &sav);
    if (k > 0) {
      fprintf(stdout, " & ");
    }
    double d = strtod(tok, NULL);
    fprintf(stdout, "%g", d);
  }
  fprintf(stdout, "\\\\ \\hline\n");
}

void parse_done(state_t * s)
{
  (void) s;
  fprintf(stdout, "\\end{tabular}");
}

int main(int argc, char * argv[])
{
  if (argc > 1) {
    state_t state = {.init = 0};
    FILE * fp = fopen(argv[1], "r"); assert(fp);
    size_t alen = 2, rlen = alen, rd;
    char * buf = malloc(alen); assert(buf);
    char * ap = buf;
    do {
      rd = fread(ap, 1, rlen, fp);
      for (char * p = ap, * end = ap + rd; p < end; p++) {
        if (*p == '\n') {
          size_t parsed = p - buf;
          buf[parsed] = 0;
          parse_line(&state, parsed, buf);

          size_t rparsed = p - ap + 1;
          size_t rest = rd  - rparsed;
          memmove(buf, ap + rparsed, rest);
          rlen = alen - rest;
          ap = buf;
          rd = rest;
        }
      }
      rlen -= rd;
      ap += rd;
      if (rlen == 0){
        size_t newa =  alen * 2;
        char * tmp = realloc(buf, newa); assert(tmp);
        ap = tmp + (ap - buf);
        buf = tmp;
        rlen = newa - alen;
        alen = newa;
      }
    } while(!feof(fp));
    parse_done(&state);
    fclose(fp);
    free(buf);
  }
}