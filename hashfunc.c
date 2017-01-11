#include "hashfunc.h"

/* hash functions */
inline hash_t dict_nbhd_key_hash(const void *_key)
{
	nbhd_key_t *pkey = (nbhd_key_t *) _key;
	return (hash_t) dict_uint64_hash_raw((void *) pkey->pmask) ^ 
		dict_uint64_hash_raw((void *) pkey->masked_id);
}

inline hash_t dict_uint64_hash_raw(const void *_key)
{
	uint64_t key = (uint64_t) _key;
	key = (~key) + (key << 18); // key = (key << 18) - key - 1;
	key = key ^ (key >> 31);
	key = key * 21; // key = (key + (key << 2)) + (key << 4);
	key = key ^ (key >> 11);
	key = key + (key << 6);
	key = key ^ (key >> 22);
	return (hash_t) key;
}

inline hash_t dict_uint32_hash_raw(const void *_key)
{
	uint32_t key = ((uintptr_t) _key) & UINT32_MAX;

    key += ~(key << 15);
    key ^=  (key >> 10);
    key +=  (key << 3);
    key ^=  (key >> 6);
    key += ~(key << 11);
    key ^=  (key >> 16);
    return key;
}

/* comparison function */
inline int dict_compare_identity(const void *key1, const void *key2)
{
	return (key1 == key2);
}

inline int dict_compare_nbhd_key(const void *key1, const void *key2)
{
	register nbhd_key_t *pk1 = (nbhd_key_t *) key1;
	register nbhd_key_t *pk2 = (nbhd_key_t *) key2;

	return ((pk1->masked_id == pk2->masked_id) && (pk1->pmask == pk2->pmask));
}
