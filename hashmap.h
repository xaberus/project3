#ifndef _HASHMAP_H
#define _HASHMAP_H

#include <stdlib.h>
#include <stdint.h>

typedef struct hashmap hashmap_t;

typedef void * (*hashmap_iterfn)(const char * key, void * value);

hashmap_t * hashmap_new(uint32_t bins);
void hashmap_insert(hashmap_t * h, const char * key, void * value);
void * hashmap_get(hashmap_t * h, char *key);
void * hashmap_delete(hashmap_t * h, char *key);
void hashmap_foreach(hashmap_t * h, hashmap_iterfn fn);
void hashmap_free(hashmap_t * h);


#endif /* _HASHMAP_H */