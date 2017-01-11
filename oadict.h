#ifndef __OA_DICT_H__
#define __OA_DICT_H__

/* An open addressing dictionary implementation */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>

enum {DICT_OK, DICT_ERR, DICT_KEY_EXISTS, DICT_ADDED, DICT_REPLACED};

#define DICT_MINSIZE 8
#define PERTURB_SHIFT 5

#define DICT_SIZE(d) ((d)->sizemask + 1)
#define DICT_COMPARE_KEYS(d, key1, key2) \
	(((d)->type.keycomp) ? \
	 (d)->type.keycomp(key1, key2) : \
	 (key1) == (key2))

#define DICT_HASH_KEY(d, key) (d)->type.hashfunc(key)

typedef uint32_t dictsize_t;
typedef uint32_t hash_t;

typedef struct dicttype_ dicttype;
typedef struct dict_ dict;
typedef struct dictentry_ dictentry;

struct dictentry_ {
	void *key;
	void *value;
};

struct dicttype_ {
	hash_t (*hashfunc)(const void *key);
	int (*keycomp)(const void *key1, const void *key2);
};

struct dict_ {
	dictsize_t sizemask;
	dictsize_t used;
	dicttype type;

	dictentry *table;
};


/* ------ Functions ----- */
dict *dict_create(register int init_size, 
		hash_t (*hashfunc)(const void *key),
		int (*keycomp)(const void *key1, const void *key2));
void dict_destroy(dict *d);
void dict_reset(dict *d);
int dict_add(dict *d, void *key, void *value);
void dict_try_resize(dict *d);
dictentry *dict_add_find(dict *d, void *key, void *val);
dictentry *dict_find(dict *d, void *key);

hash_t dict_general_hash(const void *key, size_t len);

#endif
