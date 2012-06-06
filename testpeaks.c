#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "peaks.h"

int main(int argc, char * argv[])
{
  if (argc > 1) {
    FILE * fp = fopen(argv[1], "r"); assert(fp);
    array_t * data = array_new_sized(0, 100000);
    do {
      double a, b, c, d;
      fscanf(fp, "%le %le %le %le\n", &a, &b, &c, &d);
      data = array_append(data, d);
    } while(!feof(fp));
    fclose(fp);

    array_t * peaks = peaks_find(data, 6, 3);

    for (int k = 0; k < peaks->length; k++) {
      fprintf(stdout, "%g\n", 1.25663706143591736e-01 * (peaks->data[k] - data->length/2.0));
    }
    free(peaks);
    free(data);
  }
}