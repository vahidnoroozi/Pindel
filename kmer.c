#include "kmer.h"
#include "data.h"
#include "mempool.h"
#include "oadict.h"

int kmer_incr_count_by_id(data *d, kbits_t kid, kbits_t next_base, 
		kbits_t next_2bases, kmer_t **kaddr, uint32_t label)
{
	dictentry *de;
	int count;
	kmer_t *k;

	de = dict_add_find(d->kmer_dict, kbits_cast_to_ptr(kid), NULL);
	if (de->value == NULL) {
		/* new kmer which can't be found in hash table */
		kmer_t *new_kmer = kmer_create(d, kid);
		*kaddr = new_kmer;	/* pass out the addr of newly created kmer */

		if (new_kmer == NULL) {
			return 0; 
		}

		if (next_base < 4) {
			new_kmer->trans_flag |= (1U << next_base);
		}

		if (next_2bases < 16) {
			new_kmer->trans_flag |= (1U << (4 + next_2bases));
		}

		new_kmer->label = label;
		new_kmer->dirty = 0;

		/* attach kmer to dictionary entry */
		de->value = new_kmer;

		count = 1;
	}
	else {
		*kaddr = de->value;

		k = (kmer_t *) de->value;
		count = ++(k->count);

		if (next_base < 4) {
			k->trans_flag |= (1U << next_base);
		}

		if (next_2bases < 16) {
			k->trans_flag |= (1U << (4 + next_2bases));
		}
	}

	return count; 
}

void kmer_count_trans_by_id(data *d, kbits_t kid, kbits_t next_base) 
{
	dictentry *de = dict_find(d->kmer_dict, kbits_cast_to_ptr(kid));
	if (de == NULL || de->value == NULL) return;

	kmer_t *k = de->value;
	/* convert base(s) indices into indices of the counter */
	int next_base_idx = _mm_popcnt_u64(((1UL << next_base) - 1) & k->trans_flag);

	++k->exp_count;
	++(d->total_kmer_exp_counts);
	++k->exp_trans_p[next_base_idx];
	//k->exp_trans_p[next_2bases_idx]++;
}
