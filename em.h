#ifndef __PREMIER_EM_H__
#define __PREMIER_EM_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>

#include "bitarray.h"
#include "premier.h"

/* FIXME: avoid branching */
/* FIXME: moving to numeric.h */
#define update_max_quantity(val, idx, max_val, amax_val)	\
	do {													\
		if (isnan(max_val) || (max_val < val)) {			\
			max_val = val;									\
			amax_val = idx;									\
		}													\
	} while (0);

/* FIXME: consider using lookup table instead */
/*
#define nbhd_alloc_size(n)	\
		(((n * ( + sizeof(double) + sizeof(kbits_t))) << 1) +  \
		 (n * sizeof(intptr_t) * 3) +                          \
		 (((BITS_TO_U64(n << 1)) + 2 * (BITS_TO_U64(n << 2)) + \
		  (BITS_TO_U64(n << 3))) << 3) + (n * sizeof(int)))
*/
static inline size_t nbhd_alloc_size(int n)
{
	return ((((BITS_TO_U64(n << 5)) + (BITS_TO_U64(n << 6))) << 3) +
			(sizeof(double) * 4 + sizeof(kbits_t) * 2 +
			 sizeof(int) + sizeof(uint64_t) + sizeof(intptr_t) * 2) * n);
}

void EM(data *d, struct timeval *ts);
void hmm_e_step_wrapper(data *d, read_t *read, void *fdata);
void hmm_e_step(data *d, read_t *read, void *fdata);
char *convert_quality_scores(const char *qscore, int len, int offset,
		mempool_t *mp);

void hmm_viterbi_wrapper(data *d, read_t *read, void *fdata);

#if 0
#define hmm_iterate_forward(n_sfx, par, algo_name)						\
	do {																\
		state_nbhd_t *curr_nbhd = par.curr_nbhd;						\
		state_nbhd_t *next_nbhd = par.next_nbhd;						\
		int kmer_idx = 0, idx_pfx = 0;									\
		kbits_t *t_states_sorted_sfx = curr_nbhd->states_sorted_sfx;	\
		double **tnext_kmer_trans_p = next_nbhd->kmer_trans_p;			\
																		\
		for (int i = 0; i < n_sfx; i++) {								\
			int n_common_sfx_kmers = bitarray_get_pwr2(					\
					curr_nbhd->ba_distinct_sfx, i, 1) + 1;				\
			kbits_t mut_flag = ((t_states_sorted_sfx[kmer_idx] ^		\
						par.obs_kmer) | par.obs_n_flag) & 3UL;			\
			kbits_t sfx_hd = bitarray_get_pwr2(							\
					curr_nbhd->ba_hamming_dist, i, 2);					\
			kbits_t common_sfx = t_states_sorted_sfx[kmer_idx] >> 2;	\
																		\
			kbits_t _tf = bitarray_get_pwr2(curr_nbhd->ba_pfx_sfx_flags,\
					i, 3);												\
			kbits_t ups_trans_packed = trans_flag2packed(_tf & 15);		\
			kbits_t dns_trans_packed = trans_flag2packed(_tf >> 4);		\
																		\
			int j, k;													\
			int index_pfx = 0;											\
																		\
			register kbits_t mismatch_flag = ~(1UL <<					\
					par.next_obs_base) | ~(par.next_obs_n_flag - 1UL);	\
																		\
			for (j = dns_trans_packed & 7, dns_trans_packed >>= 3;		\
					j > 0; --j, dns_trans_packed >>= 2) {				\
				kbits_t dns_base = dns_trans_packed & 3;				\
				kbits_t tnext_kid = common_sfx |						\
					(dns_base << par.shift);							\
				/* compute the position sorted by suffix */				\
				int index_sfx = cum_actg_trans_count[dns_base] -		\
					actg_trans_count[dns_base];							\
				actg_trans_count[dns_base]--;							\
																		\
				double alpha = hmm_compute_##algo_name##(				\
						ups_trans_packed,								\
						&t_state_alphas[kmer_idx],						\
						&t_kmer_trans_p[kmer_idx],						\
						dns_base, par.emit_dens);						\
																		\
				update_max_quantity(alpha, index_sfx,					\
						tnext_max_alpha, tnext_argmax_alpha);			\
																		\
				kmer_t *tnext_kmer = dict_find(kdict,					\
						CAST_TO_PTR(tnext_kid))->value;					\
																		\
				tnext_kmer_ptrs[index_sfx] = tnext_kmer;				\
				tnext_kmer_trans_p[index_sfx] = tnext_kmer->trans_p;	\
				/* set up transition flag */							\
				bitarray_set_pwr2(tnext_kmer_trans_flag, index_sfx,		\
						tnext_kmer->trans_flag, 2);						\
				/* set up Hamming distance */							\
				bitarray_set_pwr2(tnext_hamming_dist, index_sfx,		\
						sfx_hd + ((mismatch_flag >> j) & 1), 2);		\
																		\
				tnext_sfx_order[index_pfx] = index_sfx;					\
				tnext_states_sorted_pfx[index_pfx++] = tnext_kid;		\
				tnext_states_sorted_sfx[index_sfx] = tnext_kid;			\
				tnext_alphas[index_sfx] = alpha;						\
			}															\
																		\
			kmer_idx += n_common_sfx_kmers;								\
		}																\
	} while (0);
#endif

typedef struct candidate_s {
	int pfx_order;
	int is_insertion;
	kmer_t *kmer_ptr;
	kbits_t id;
	uint64_t edist;
	double emit_dens;
} candidate_t;

#endif
