#ifndef _ARRAY_H_
#define _ARRAY_H_

typedef struct array {
  int    alloc;
  int    length;
  double data[];
  // field elements are allocated after the array
} array_t;

array_t * array_new(int length);
array_t * array_new_sized(int length, int alloc);
array_t * array_equipart(double start, double end, int points);
array_t * array_append(array_t * z, double x);

#endif /* _CARRAY_H_ */