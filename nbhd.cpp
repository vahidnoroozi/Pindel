#include <algorithm>
#include <list>
#include <iostream>

extern "C" {
#include "premier.h"
#include "kvec.h"
#include "kmer.h"
#include "numeric.h"
}

static int kmer_suffix_comp(const void *a, const void *b)
{
	register kmer_t *ka = *((kmer_t **) a);
	register kmer_t *kb = *((kmer_t **) b);
	register kbits_t sfx_a = ka->id, sfx_b = kb->id;

	return (sfx_a > sfx_b) ? 1 : -1;
}

void hmm_build_1st_nbhd(data *d, read_t *read, void *fdata) 
{
	int dmax = 2;

	int kmer_len = d->opt->kmer_length;
	double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);

	if (read->length <= kmer_len) return;

	kbits_t obs_kid = kmer_seq_to_id(read->sequence, kmer_len, 0);
	kmer_t *pk = (kmer_t *) dict_find(d->kmer_dict, 
				kbits_cast_to_ptr(obs_kid))->value;

	/* # of times a kmer and its reverse complement are observed */
	double k_cexp = EXP(pk->init_p) * d->total_kmer_exp_counts;

	kvec_t(kmer_t *) k_nbhd;
	kv_init(k_nbhd);

	/* observed kmer has zero Hamming distance */
	mut_kmer_t obs_mk = {obs_kid, 0};

	/* observed kmer is always in the first neighborhood */
	kv_push(kmer_t *, k_nbhd, pk);

	std::list<mut_kmer_t> k_mut_cands;

	/* expand the observed kmer if :
	 *   1) it has # occurence below threshold, or
	 *   2) it has no eligible transitions (removed in initialization) */
	if (k_cexp < thresh || pk->trans_flag == 0) {
		k_mut_cands.push_back(obs_mk);
	}

	int hamming_d = 0; /* hamming distance */
	int counter = 0;
	while (k_mut_cands.empty() == false) {
		mut_kmer_t &mk = k_mut_cands.front();

		int _d = _mm_popcnt_u64(mk.mut_flag);

		if (_d > hamming_d) {
			/* Hamming distance increased, check if we already have viable
			 * candidates in the neighborhood */
			if (kv_size(k_nbhd) > 1) break;
			hamming_d = _d;
		}

		for (int i = 0; i < kmer_len; i++) {
			int i2 = i << 1;

			/* to avoid duplicated mutations, we will only mutate position i,
			 * if no bits higher than i are set to 1. */
			if (mk.mut_flag & ((1UL << (i+1)) - 1UL)) continue;

			for (kbits_t base = 0; base < 4; ++base) {
				if (base == ((mk.id >> i2) & 3)) continue;
				if (base == ((obs_kid >> i2) & 3)) continue;

				kbits_t mut_kid = (mk.id & ~(3UL << i2)) | (base << i2);
				kbits_t mut_flag = mk.mut_flag | (1UL << i);

				dictentry *de_mut = dict_find(d->kmer_dict, 
						kbits_cast_to_ptr(mut_kid));

				if (de_mut == NULL && (_d+1) < dmax) {
					/* continue searching */
					mut_kmer_t new_mk = {mut_kid, mut_flag};
					k_mut_cands.push_back(new_mk);
				}
				else if (de_mut != NULL) {
					kmer_t *mut_pk = (kmer_t *) de_mut->value;

					double mut_k_cexp = EXP(mut_pk->init_p) * 
						d->total_kmer_exp_counts;

					if (mut_k_cexp > thresh && mut_pk->trans_flag != 0) {
						kv_push(kmer_t *, k_nbhd, mut_pk);
					}
					else if ((_d+1) < dmax) {
						mut_kmer_t new_mk = {mut_kid, mut_flag};
						k_mut_cands.push_back(new_mk);
					}
				}
			}
		}

		k_mut_cands.pop_front();
	}

	/* sort the nbhd by suffix */
	register int size_nbhd = kv_size(k_nbhd);
	kmer_t **new_nbhd = NULL;

	new_nbhd = (kmer_t **) mempool_alloc(d->mp, 
			size_nbhd * sizeof(kmer_t *));
	memcpy(new_nbhd, &(kv_A(k_nbhd, 0)), size_nbhd * sizeof(kmer_t *));

	qsort(new_nbhd, size_nbhd, sizeof(kmer_t *), kmer_suffix_comp);

	d->preconstructed_nbhds[read->id] = (kmer_nbhd_t) {size_nbhd, new_nbhd};

	kv_destroy(k_nbhd);
}
