#include "array.h"

array_t * array_new(int length) {
  array_t * r = malloc(sizeof(array_t) + sizeof(double) * length);
  if (!r) { abort(); }
  r->length = length;
  r->alloc = length;
  return r;
}

// creates a array with an equipatition of the
// intervall [start;end] with length-1 subintevalls
array_t * array_equipart(double start, double end, int length) {
  int k;
  double range = end - start;
  array_t * r = array_new(length);

  for (k = 0; k < length; k++) {
    r->data[k] = start + range / (length - 1) * k;
  }

  return r;
}