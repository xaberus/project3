#ifndef _ARRAY_H_
#define _ARRAY_H_

/*! A rezizeable vector of doubles */
typedef struct array {
  int    alloc;   /**< allocated length */
  int    length;  /**< used length */
  /*! c99 flexible array tail with \a alloc members
      (field elements are allocated after the header) */
  double data[];
} array_t;

array_t * array_new(int length);
array_t * array_new_sized(int length, int alloc);
array_t * array_equipart(double start, double end, int points);
array_t * array_append(array_t * z, double x);

#endif /* _CARRAY_H_ */
