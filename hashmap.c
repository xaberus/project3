#include <string.h>
#include <assert.h>

#include "hashmap.h"
#include "hash.h"

typedef struct hashbin hashbin_t;

struct hashbin {
  uint32_t    hash;
  char      * key;
  void      * value;
  hashbin_t * next;
};

struct hashmap {
  uint32_t     size;
  hashbin_t ** bins;
};

hashmap_t * hashmap_new(uint32_t bins) {
  hashmap_t * h = malloc(sizeof(hashmap_t)); assert(h);
  h->size = bins;
  h->bins = malloc(sizeof(hashbin_t *) * bins); assert(h->bins);
  memset(h->bins, 0, sizeof(hashbin_t *) * bins);
  return h;
}

void hashmap_insert(hashmap_t * h, const char * key, void * value) {
  if (value == NULL) {
    return;
  }

  uint32_t o = hash(key, strlen(key), 777);
  uint32_t i = o % h->size;

  hashbin_t ** p = &h->bins[i], * b = *p;
  while (b) {
    if (o == b->hash && strcmp(key, b->key) == 0) {
      b->value = value;
      return;
    }
    p = &b->next;
    b = *p;
  }
  hashbin_t * n = malloc(sizeof(hashbin_t)); assert(n);
  n->hash = o;
  n->key = strdup(key); assert(n->key);
  n->value = value;
  n->next = NULL;
  *p = n;
}

void * hashmap_get(hashmap_t * h, char * key) {
  uint32_t o = hash(key, strlen(key), 777);
  uint32_t i = o % h->size;

  hashbin_t *b = h->bins[i];
  while (b) {
    if (o == b->hash && strcmp(key, b->key) == 0) {
      return b->value;
    }
    b = b->next;
  }
  return NULL;
}
void * hashmap_delete(hashmap_t * h, char * key) {
  uint32_t o = hash(key, strlen(key), 777);
  uint32_t i = o % h->size;

  hashbin_t ** p = &h->bins[i], * b = *p;
  while (b) {
    if (o == b->hash && strcmp(key, b->key) == 0) {
      void * r = b->value;
      *p = b->next;
      free(b->key);
      free(b);
      return r;
    }
    p = &b->next;
    b = *p;
  }
  return NULL;
}
void hashmap_free(hashmap_t * h) {
  for (uint32_t k = 0; k < h->size; k++) {
    hashbin_t * b = h->bins[k];
    while (b) {
      hashbin_t * n = b->next;
      free(b->key);
      free(b);
      b = n;
    }
  }
  free(h->bins);
  free(h);
}

void hashmap_foreach(hashmap_t * h, hashmap_iterfn fn) {
  for (uint32_t k = 0; k < h->size; k++) {
    hashbin_t ** p = &h->bins[k], * b = *p;
    while (b) {
      void * r = fn(b->key, b->value);
      if (!r) {
        *p = b->next;
        free(b->key);
        free(b);
        b = *p;
      } else {
        p = &b->next;
        b = *p;
      }
    }
  }
}

#if 0

void * free_value(const char * key, void * value) {
  free(value);
  return NULL;
}

int main() {
  char buf[128];
  hashmap_t * h = hashmap_new(256);
  for (uint32_t k = 0; k < 100; k++) {
    snprintf(buf, sizeof(buf), "%-15u", k);
    hashmap_insert(h, buf, (void *) malloc(123));
  }
  hashmap_foreach(h, free_value);
  hashmap_free(h);
}

#endif