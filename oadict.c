#include "oadict.h"

static hash_t _dict_identity_hash_func(const void *key);
static int _dict_insert(dict *d, void *key, void *val);
static int _dict_resize(dict *d, dictsize_t minused);
static dictentry* _dict_lookup(dict *d, void *key, register hash_t hash);

static uint32_t dict_hash_function_seed = 5381;

dict *dict_create(register int init_size, 
		hash_t (*hashfunc)(const void *key),
		int (*keycomp)(const void *key1, const void *key2))
{
	register int _size;
	register dictentry *ptbl;

	if (init_size < DICT_MINSIZE) {
		init_size = DICT_MINSIZE;
	}

	for (_size = DICT_MINSIZE; _size < init_size && _size > 0; _size <<= 1);
	if (_size < 0){
		return NULL;
	}

	register dict *newd = (dict *) malloc(sizeof(*newd));
	if (newd == NULL){
		return NULL;
	}

	ptbl = (dictentry *) calloc(_size, sizeof(*ptbl));
	if (ptbl == NULL) {
		free(newd);
		return NULL;
	}

	newd->table = ptbl;
	newd->sizemask = _size - 1;
	newd->used = 0;
	/* if user doesn't specify a hash function, use the identity function,
	 * which use the address of the key pointer as the hash. */
	newd->type.hashfunc = hashfunc ? hashfunc : _dict_identity_hash_func;
	newd->type.keycomp = keycomp;

	return newd;
}

void dict_reset(dict *d)
{
	memset(d->table, 0, sizeof(*d->table) * (d->sizemask + 1));
	d->used = 0;
}

void dict_destroy(dict *d)
{
	free(d->table);
	free(d);
}

int dict_add(dict *d, void *key, void *val)
{
	register dictsize_t n_used = d->used;
	if (_dict_insert(d, key, val) != 0) {
		return -1;
	}

	if (!(d->used > n_used && d->used * 3 >= (d->sizemask+1) * 2)) {
		return 0;
	}
	return _dict_resize(d, d->used << (d->used < 1024 * 1024 ? 4 : 2));
}

dictentry *dict_add_find(dict *d, void *key, void *val)
{
	register hash_t hash = DICT_HASH_KEY(d, key);
	register dictentry *de;

	/* pre-check if the dictionary is on the verge of resizing,
	 * do it before inserting any entries to avoid data corruption. */
	if ((d->used + 1) * 3 >= (d->sizemask+1) * 2){
		_dict_resize(d, (d->used+1) << (d->used < 1024 * 1024 ? 2 : 1));
	}

	de = _dict_lookup(d, key, hash);
	if (de == NULL) {
		return NULL;
	}

	/* insert an new entry */
	if (de->key == NULL){
		de->key = key;
		de->value = val;

		d->used++;
		return de;
	}

	return de;
}

dictentry *dict_find(dict *d, void *key)
{
	register hash_t hash = DICT_HASH_KEY(d, key);
	register dictentry *de;

	de = _dict_lookup(d, key, hash);
	if (de == NULL || de->value == NULL) {
		return NULL;
	}

	return de;
}

void dict_try_resize(dict *d)
{
#pragma omp critical
{
	if ((d->used + 1) * 3 >= (d->sizemask+1) * 2) {
		_dict_resize(d, (d->used+1) << (d->used < 1024 * 1024 ? 2 : 1));
	}
}
}

/* general purpose hash function from Redis */
hash_t dict_general_hash(const void *key, size_t len)
{
    /* 'm' and 'r' are mixing constants generated offline.
     They're not really 'magic', they just happen to work well.  */
    uint32_t seed = dict_hash_function_seed;
    const uint32_t m = 0x5bd1e995;
    const int r = 24;

    /* Initialize the hash to a 'random' value */
    uint32_t h = seed ^ len;

    /* Mix 4 bytes at a time into the hash */
    const unsigned char *data = (const unsigned char *) key;

    while(len >= 4) {
        uint32_t k = *(uint32_t*) data;

        k *= m;
        k ^= k >> r;
        k *= m;

        h *= m;
        h ^= k;

        data += 4;
        len -= 4;
    }

    /* Handle the last few bytes of the input array  */
    switch(len) {
		case 3: h ^= data[2] << 16;
		case 2: h ^= data[1] << 8;
		case 1: h ^= data[0]; h *= m;
    };

    /* Do a few final mixes of the hash to ensure the last few
     * bytes are well-incorporated. */
    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;

    return (hash_t) h;
}

static int _dict_resize(dict *d, dictsize_t minused)
{
	register dictsize_t newsize;
	register dictsize_t i;
	register dictentry *newtable, *oldtable, *de;

	for (newsize = DICT_MINSIZE; newsize <= minused && newsize > 0; 
			newsize <<= 1);

	if (newsize < 0){
		/* dictionary size out of bound */
		return -1;
	}

	newtable = calloc(newsize, sizeof(*newtable));
	if (newtable == NULL) {
		return -1;
	}

	oldtable = d->table;
	i = d->used;
	/* replace hash table */
	d->table = newtable;
	d->used = 0;
	d->sizemask = newsize - 1;

	for (de = oldtable; i > 0; de++){
		if (de->value != NULL) {
			i--;
			_dict_insert(d, de->key, de->value);
		}
	}

	free(oldtable);
	return 0;
}

static int _dict_insert(dict *d, void *key, void *val)
{
	register hash_t hash = DICT_HASH_KEY(d, key);
	register dictentry *de;

	de = _dict_lookup(d, key, hash);
	if (de == NULL) {
		return -1;
	}

	if (de->value != NULL){
		/* should we free/release the memory of the old value ? */
		de->value = val;
	}
	else {
		if (de->key == NULL){
			d->used++;
			de->key = key;
			de->value = val;
		}
	}

	return 0;
}


static dictentry* _dict_lookup(dict *d, void *key, register hash_t hash)
{
	register dictsize_t i;
	register dictsize_t perturb;
	register dictsize_t mask = d->sizemask;
	//register dictentry *freeslot;
	register dictentry *de;

	dictentry *ptbl = d->table;

	i = (dictsize_t) hash & mask;
	de = &ptbl[i];
	
	if (de->key == NULL || de->key == key) {
		return de;
	}
	else {
		if (DICT_COMPARE_KEYS(d, de->key, key)) {
			return de;
		}
	}

	for (perturb = hash; ; perturb >>= PERTURB_SHIFT) {
		i = (i << 2) + i + perturb + 1;
		de = &ptbl[i & mask];

		if (de->key == NULL || de->key == key){
			return de; 
		}

		if (DICT_COMPARE_KEYS(d, de->key, key)) {
			return de;
		}
	}

	return NULL;
}

static hash_t _dict_identity_hash_func(const void *key)
{
	return (hash_t) ((intptr_t) key);
}
