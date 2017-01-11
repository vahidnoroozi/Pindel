#ifndef __PREMIER_HASHFUNC_H__
#define __PREMIER_HASHFUNC_H__

#include <stdint.h>

#include "premier.h"
#include "oadict.h"

hash_t dict_uint32_hash_raw(const void *_key);
hash_t dict_uint64_hash_raw(const void *_key);
hash_t dict_nbhd_key_hash(const void *_key);
int dict_compare_identity(const void *key1, const void *key2);
int dict_compare_nbhd_key(const void *key1, const void *key2);

#endif
