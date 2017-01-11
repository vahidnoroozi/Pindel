#ifndef __EM_NATURAL_SCALE_H__
#define __EM_NATURAL_SCALE_H__

#include "em.h"

/* Create a pair for a given kmer k, and its following downstream base,
 * dns_base. Then insert this pair into the dictionary.
 */
static void kmer_add_to_pair(kmer_t *k, kbits_t kid, 
		kbits_t *hamming_dist_v, DOUBLE *norm_trans_p, 
		DOUBLE inv_sum_tp, DOUBLE alph, 
		kmer_pair_t *link_kpair, dict *kp_dict, 
		kbits_t dns_path_flag, mempool_t *mp) 
{
	INIT_TRANS_LOOKUP_TABLE(_kpair)

	dictentry *de;
	kmer_pair_t *kpair;
	kbits_t k_sfx = KBITS_SUFFIX(kid);
	register kbits_t ups_base = kid & 3;

	de = dict_add_find(kp_dict, CAST_TO_PTR(k_sfx), NULL);
	if (de->value == NULL){
		/* --- create a kmer pair --- */
		/* kmer_pair_t must be aligned by 16 bytes! */
		kpair = mempool_nalloc(mp, sizeof(*kpair), 16);

		memset(kpair, 0, sizeof(*kpair));
		memset(kpair->dns_hamming_dist, 0xFF, sizeof(kbits_t) << 2); 
		kpair->common_substr = k_sfx;

		de->value = kpair;
	}
	else {
		kpair = (kmer_pair_t *) de->value;
	}

	kpair->dns_path_flag |= dns_path_flag;
	/* record upstream references/values */
	kpair->ups_alpha[ups_base] = alph;
	kpair->ups_ref_kmers[ups_base] = k;
	kpair->links[ups_base] = link_kpair;
	kpair->ups_path_flag |= (1 << ups_base);

	/* set up downstream values */
	DOUBLE *_dns_trans_p = kpair->dns_st_trans_p;

	register int i;
	register kbits_t nonzero_trans = _kpair_trans_lookup_table[dns_path_flag];

	for (i = nonzero_trans & 7, nonzero_trans >>= 3; i > 0;
			i--, nonzero_trans >>= 2) {
		int dns_base = nonzero_trans & 3;
		kpair->dns_hamming_dist[dns_base] = hamming_dist_v[dns_base];

		_dns_trans_p[(dns_base << 2) + ups_base] = norm_trans_p[dns_base];
	}
}


static inline DOUBLE compute_loglik(kmer_pair_t *kpair, int shift, 
		kbits_t dmax, dict *kdict, DOUBLE *ptr_emit_p)
{
	register kbits_t k_sfx = kpair->common_substr;

	DOUBLE alpha __attribute__ ((aligned (16)));
	/* log-likelihood of current HMM */

	register DOUBLE sum_ll = 0.0;

	DOUBLE *ptr_ups_alpha = kpair->ups_alpha;
	DOUBLE *ptr_dns_st_trans_p = kpair->dns_st_trans_p;
	DOUBLE *ptr_dns_gamma = kpair->dns_gamma;
	DOUBLE *ptr_dns_beta = kpair->dns_beta;
	kmer_t **ptr_dns_ref_kmers = (kmer_t **) kpair->dns_ref_kmers;
	kbits_t *ptr_dns_hd = kpair->dns_hamming_dist;

	/* downstream states */
	for (int i = 0; i < 4; i++) {
		if (ptr_dns_hd[i] > dmax) continue;

		DOUBLE *ptr_trans_p = &ptr_dns_st_trans_p[i << 2];

		/* compute forward recursion relation */
		compute_alpha(ptr_ups_alpha, ptr_trans_p, ptr_emit_p[i], &alpha);
		sum_ll += alpha; 

		/* set up values/references for backward iteration */
		kbits_t nb_kbits = k_sfx | ((kbits_t) i << shift);
		kmer_t *nb_kmer = (kmer_t *) (dict_find(kdict, 
					CAST_TO_PTR(nb_kbits))->value);

		ptr_dns_ref_kmers[i] = nb_kmer;
		ptr_dns_beta[i] = 1.0;
		ptr_dns_gamma[i] = alpha;
	}

	return sum_ll;
}

static inline DOUBLE backward_iteration(kmer_pair_t *kpair, 
		DOUBLE *ptr_emit_p, int km1_shift, int dmax)
{
	register kbits_t k_sfx = kpair->common_substr;

	DOUBLE *ptr_dns_beta = kpair->dns_beta;
	DOUBLE *ptr_ups_alpha = kpair->ups_alpha;
	kbits_t *ptr_dns_hd = kpair->dns_hamming_dist;
	kmer_pair_t **ptr_up_links = kpair->links;
	kmer_t **ptr_ups_kmers = (kmer_t **) kpair->ups_ref_kmers;

	/* normalized transition probability */
	static DOUBLE norm_trans_p[4] __attribute__ ((aligned (16)));

	/*
	 *  A -\   #/-  A
	 *  C --TTGC--  C
	 *  G -/    \-  G
	 *  T -|    |-  T
	 *
	 *  the C base with # on top of it is the downstream nucleotide at previous
	 *  step, which uses [ACGT]TTG as the (k-1)-substring.
	 */
	register kbits_t sfx_kth_base = (k_sfx >> km1_shift);
	register kbits_t sfx_kth_base_sll2 = sfx_kth_base << 2;

	register DOUBLE sum_exp_trans = 0.0;
	DOUBLE _beta;

	/* upstream states */
	kbits_t nonzero_trans = _bwd_trans_lookup_table[kpair->ups_path_flag];
	int n_nz_trans = nonzero_trans & 7;
	nonzero_trans >> 3;

	for (int i = 0; i < 4; i++){
		kmer_pair_t *_up_link = ptr_up_links[i];
		if (_up_link == NULL) continue;

		DOUBLE *ptr_norm_tp;

		DOUBLE *ptr_trans_p = ptr_ups_kmers[i]->transition_p;
		kbits_t utf = ptr_ups_kmers[i]->trans_flag;

		/* if no transition is crossed out, skip normalization */
		if (utf == kpair->dns_path_flag) {
			ptr_norm_tp = ptr_trans_p;
		}
		else {
			normalize_trans_prob_hd(ptr_trans_p, ptr_dns_hd, dmax,
					(DOUBLE *) norm_trans_p);
			ptr_norm_tp = norm_trans_p;
		}

		DOUBLE *exp_trans = &ptr_up_links[i]->exp_trans_p[
			sfx_kth_base_sll2];

		register DOUBLE _expcnt = compute_beta(ptr_dns_beta, ptr_norm_tp, 
				ptr_emit_p, ptr_ups_alpha[i], 
				&_beta, exp_trans);	

		sum_exp_trans += _expcnt;

		_up_link->dns_beta[sfx_kth_base] = _beta; 
		_up_link->dns_gamma[sfx_kth_base] = _expcnt; 
		/* compute beta_{t} recursively */
	}

	return sum_exp_trans;
}

/* compute the last beta iteration (at the first kmer) */
static inline DOUBLE final_backward_iteration(kmer_pair_t *kpair, 
		DOUBLE *ptr_emit_p, int dmax)
{
	DOUBLE *ptr_dns_beta = kpair->dns_beta;
	kbits_t *ptr_dns_hd = kpair->dns_hamming_dist;
	DOUBLE *ptr_ups_alpha = kpair->ups_alpha;
	DOUBLE *ptr_exp_trans_p = kpair->exp_trans_p;
	kmer_t **ptr_ups_kmers = (kmer_t **) kpair->ups_ref_kmers;

	/* normalized transition probability */
	static DOUBLE norm_trans_p[4] __attribute__ ((aligned (16)));

	/*
	 *  A -\   #/-  A
	 *  C --TTGC--  C
	 *  G -/    \-  G
	 *  T -|    |-  T
	 *
	 *  the C base with # on top of it is the downstream nucleotide at previous
	 *  step, which uses [ACGT]TTG as the (k-1)-substring.
	 */

	/* set all expected transition prob to zero to avoid contamination from
	 * previous iteration */
	memset(ptr_exp_trans_p, 0, sizeof(*ptr_exp_trans_p) << 4);

	DOUBLE sum_exp_trans = 0.0;
	DOUBLE _beta;
	/* upstream states */
	for (int i = 0; i < 4; i++){
		/* hidden state not allowed */
		if (ptr_ups_kmers[i] == NULL) continue;

		//kbits_t nb_kid = k_sfx_sll2 | i;

		DOUBLE *ptr_trans_p = ptr_ups_kmers[i]->transition_p;
		DOUBLE *ptr_norm_tp;
		/* normalize probablity using betas as mask */
		if (ptr_ups_kmers[i]->unique_trans) {
			ptr_norm_tp = ptr_trans_p;
		}
		else {
			normalize_trans_prob_hd(ptr_trans_p, ptr_dns_hd, dmax,
					(DOUBLE *) norm_trans_p);
			ptr_norm_tp = norm_trans_p;
		}

		DOUBLE *exp_trans = &ptr_exp_trans_p[i << 2];

		/* FIXME:
		 * the usage of array norm_trans_p is totally redundant and
		 * inefficient, because one has to store the normalized probabilities
		 * to some memory block and later load them again in compute_beta. 
		 *
		 * consider normalizing the probabilities inside compute_beta by
		 * keeping all variables in the XMM registers necessarily.
		 **/
		DOUBLE _expcnt = compute_beta(ptr_dns_beta, ptr_norm_tp, 
				ptr_emit_p, ptr_ups_alpha[i], 
				&_beta, exp_trans);	

		ptr_ups_alpha[i] = _expcnt;
		sum_exp_trans += _expcnt;

		/* compute beta_{t} recursively */
	}
	/* compute gamma, the expected occurrence of each state/kmer */

	return sum_exp_trans;
}

static void normalize_trans_prob_hd(DOUBLE *raw_tp, kbits_t *hd, 
		int dmax, DOUBLE *norm_tp)
{
#  ifdef __SSE3__
	register __m128d x_dmax = _mm_set1_pd((DOUBLE) dmax);
	register __m128d x_mask = _mm_cmple_pd(
			_mm_set_pd((DOUBLE) hd[1], (DOUBLE) hd[0]), x_dmax);
	register __m128d x_lo = _mm_and_pd(x_mask, _mm_load_pd(raw_tp));

	x_mask = _mm_cmple_pd(
			_mm_set_pd((DOUBLE) hd[3], (DOUBLE) hd[2]), x_dmax);
	register __m128d x_hi = _mm_and_pd(x_mask, _mm_load_pd(raw_tp + 2));
	register __m128d x_sum = _mm_add_pd(x_lo, x_hi);
	x_sum = _mm_hadd_pd(x_sum, x_sum);

	register __m128d x_inv_sum = _mm_div_pd(_mm_set1_pd(1.0), x_sum);

	/* in case inv_sum is infinity */
	register __m128d x_inf = _mm_set1_pd(INFINITY);
	register __m128d x_inf_mask = _mm_cmpneq_pd(x_inf, x_inv_sum);
	x_inv_sum = _mm_and_pd(x_inf_mask, x_inv_sum);

	x_lo = _mm_mul_pd(x_lo, x_inv_sum);
	x_hi = _mm_mul_pd(x_hi, x_inv_sum);
	_mm_store_pd(norm_tp, x_lo);
	_mm_store_pd(norm_tp + 2, x_hi);
#  endif
}

static inline void compute_alpha(DOUBLE *ptr_ups_alpha, DOUBLE *norm_trans_p,
		DOUBLE emit_dens, DOUBLE *alpha)
{
#ifdef __SSE3__
		/* transition prob */
		register __m128d x_tp = _mm_load_pd(norm_trans_p); 
		/* upstream alpha's */
		register __m128d x_pa = _mm_load_pd(ptr_ups_alpha); 
		/* compute 1st part of recursion */
		register __m128d x_tmp = _mm_mul_pd(x_tp, x_pa);
		register __m128d x_emit = _mm_set1_pd(emit_dens);

		x_tp = _mm_load_pd(norm_trans_p + 2);
		x_pa = _mm_load_pd(ptr_ups_alpha + 2);
		/* compute 2nd part of recursion */
		x_tp = _mm_mul_pd(x_tp, x_pa);
		/* horizontal add */
		x_pa = _mm_add_pd(x_tmp, x_tp);
		/* now xmm_pa hold two partial sum, sum them together to get the
		 * sum over all bases */
		x_pa = _mm_hadd_pd(x_pa, x_pa);
		/* multiply by the emission probability */
		x_pa = _mm_mul_pd(x_pa, x_emit);
		_mm_storel_pd(alpha, x_pa);
#else
		/* legacy code without speedup */
#endif
}

inline DOUBLE compute_beta(DOUBLE *ptr_dns_beta, DOUBLE *ptr_trans_p, 
		DOUBLE *ptr_emit_p, DOUBLE alpha, DOUBLE *ptr_beta,
		DOUBLE *ptr_exp_trans)
{
#ifdef __SSE3__
	__m128d x_alpha = _mm_set1_pd(alpha);
	/* A/C first */
	__m128d x_trans_p = _mm_load_pd(ptr_trans_p);
	__m128d x_emit_p = _mm_load_pd(ptr_emit_p);
	__m128d x_dns_beta = _mm_load_pd(ptr_dns_beta);
	/* upstream alpha's */
	__m128d x_tail = _mm_mul_pd(x_trans_p, _mm_mul_pd(x_emit_p,
				x_dns_beta));
	/* expected transition */
	__m128d x_exp_trans = _mm_mul_pd(x_alpha, x_tail);
	_mm_store_pd(ptr_exp_trans, x_exp_trans);

	/* T/G */
	x_dns_beta = _mm_load_pd(ptr_dns_beta + 2);
	x_emit_p = _mm_load_pd(ptr_emit_p + 2);
	x_trans_p = _mm_load_pd(ptr_trans_p + 2);

	__m128d x_tail_tg = _mm_mul_pd(x_trans_p, _mm_mul_pd(x_emit_p,
						x_dns_beta));
	__m128d x_exp_trans_tg = _mm_mul_pd(x_alpha, x_tail_tg);
	_mm_store_pd(ptr_exp_trans + 2, x_exp_trans_tg);

	x_tail = _mm_add_pd(x_tail, x_tail_tg);
	x_tail = _mm_hadd_pd(x_tail, x_tail);
	_mm_storel_pd(ptr_beta, x_tail);
	
	return (*ptr_beta) * alpha;
#else
	/* TODO: non-SSE-optimized version */
#endif /* __SSE3__ */
}

/* NOTE: Although this function claims to return a pointer, it is NOT a REAL
 * pointer. Instead, the returned value is a compound value consists of the 
 * address of the best kmer pair (pointer), plus an integer (offset) between 
 * 0 and 3 (inclusive) to indicate which state is the optimal within the kmer pair.
 *
 * Because kmer_pair_t is dynamically allocated on memory with an address at
 * least 16-bytes aligned, the pointer and the offset can be easily decoupled
 * by,
 *		offset = p & 15
 *		pointer = p & ~((uinptr_t) 15)
 * if p is the value returned by this function.
 */
static inline kmer_pair_t *compute_delta(DOUBLE *ptr_ups_alpha, 
		DOUBLE *norm_trans_p, DOUBLE emit_dens, DOUBLE *delta, 
		kmer_pair_t *kpair)
{
#ifdef __SSE2__
	DOUBLE idx_argmax;

	/* transition prob */
	register __m128d x_tp = _mm_load_pd(norm_trans_p); 
	/* upstream alpha's */
	register __m128d x_pa = _mm_load_pd(ptr_ups_alpha); 
	/* compute 1st part of recursion */
	register __m128d x_delta = _mm_mul_pd(x_tp, x_pa);
	register __m128d x_emit = _mm_set1_pd(emit_dens);

	x_tp = _mm_load_pd(norm_trans_p + 2);
	x_pa = _mm_load_pd(ptr_ups_alpha + 2);
	/* compute 2nd part of recursion */
	__m128d x_delta_tg = _mm_mul_pd(x_tp, x_pa);
	/* find the two largest delta's */
	register __m128d x_ord_mask = _mm_cmplt_pd(x_delta, x_delta_tg);
	register __m128d x_max = _mm_or_pd(
			_mm_andnot_pd(x_ord_mask, x_delta), 
			_mm_and_pd(x_ord_mask, x_delta_tg));
	register __m128d x_argmax = _mm_or_pd(
			_mm_andnot_pd(x_ord_mask, _mm_set_pd(1.0, 0.0)), /* A/C */
			_mm_and_pd(x_ord_mask, _mm_set_pd(3.0, 2.0)));
	register __m128d x_max_sr = (__m128d) _mm_srli_si128((__m128i) x_max, 8);
	register __m128d x_argmax_sr = (__m128d) _mm_srli_si128(
		(__m128i) x_argmax, 8);

	x_ord_mask = _mm_cmplt_pd(
		(__m128d) _mm_srli_si128((__m128i) x_max, 8), x_max);
	x_max = _mm_or_pd(_mm_andnot_pd(x_ord_mask, x_max_sr),
		_mm_and_pd(x_ord_mask, x_max));
	/* finally, multiply by emission probability */
	x_max = _mm_mul_pd(x_max, x_emit);
	x_argmax = _mm_or_pd(_mm_andnot_pd(x_ord_mask, x_argmax_sr),
		_mm_and_pd(x_ord_mask, x_argmax));

	_mm_storel_pd(delta, x_max);
	_mm_storel_pd(&idx_argmax, x_argmax);

	return (kmer_pair_t *) ((uintptr_t) kpair + 
			(uintptr_t) idx_argmax);
#else
		/* legacy code without speedup */
#endif /* SSE2 */
}


static inline void normalize_gamma(DOUBLE *ptr_gamma, DOUBLE inv_exp_state)
{
#ifdef __SSE3__
	register __m128d x_numert = _mm_load_pd(ptr_gamma);
	register __m128d x_denomt = _mm_set1_pd(inv_exp_state);

	_mm_store_pd(ptr_gamma, _mm_mul_pd(x_numert, x_denomt));
	x_numert = _mm_load_pd(ptr_gamma + 2);
	_mm_store_pd(ptr_gamma + 2, _mm_mul_pd(x_numert, x_denomt));
#endif
}

static inline void update_gamma(kmer_t **ptr_ref_kmers, DOUBLE *ptr_gamma,
		kbits_t obs_base, DOUBLE *ptr_qual_exp_c, 
		DOUBLE *ptr_qual_exp_e, DOUBLE *ptr_emit_exp)
{
	kmer_t *k;
	for (int i = 0; i < 4; i++){
		if ((k = (ptr_ref_kmers[i])) == NULL) continue;

		DOUBLE _gm = ptr_gamma[i];

		*ptr_qual_exp_c += _gm;
		ptr_emit_exp[i] += _gm;
		if (obs_base != i){
			*ptr_qual_exp_e += _gm;
		}
	}
}

static inline void update_gamma_and_xi(kmer_t **ptr_ref_kmers, DOUBLE *ptr_gamma,
		DOUBLE *ptr_xi, DOUBLE inv_exp_trans, kbits_t obs_base, 
		DOUBLE *ptr_qual_exp_c, DOUBLE *ptr_qual_exp_e, 
		DOUBLE *ptr_emit_exp)
{
	kmer_t *k;
	register __m128d x_denomt = _mm_set1_pd(inv_exp_trans);

	register DOUBLE _gm_sum = 0.0;
	register DOUBLE _gm_sum_e = 0.0;

	for (int i = 0; i < 4; i++){
		if ((k = (ptr_ref_kmers[i])) == NULL) continue;
		register int _isll2 = i << 2;

#ifdef __SSE3__
		DOUBLE *ptr_kmer_xi = k->xi;
		register dbl_t _gm = {.dbl = ptr_gamma[i] * inv_exp_trans};

		register __m128d x_op_b = _mm_mul_pd(
				_mm_load_pd(&ptr_xi[_isll2]), x_denomt);
		register __m128d x_op_a = _mm_load_pd(ptr_kmer_xi);

		x_op_a = _mm_add_pd(x_op_a, x_op_b);
		_mm_store_pd(ptr_kmer_xi, x_op_a);

		x_op_a = _mm_load_pd(ptr_kmer_xi + 2);
		x_op_b = _mm_mul_pd(_mm_load_pd(&ptr_xi[_isll2 + 2]), x_denomt);
		x_op_a = _mm_add_pd(x_op_a, x_op_b);
		_mm_store_pd(ptr_kmer_xi + 2, x_op_a);

		k->gamma += _gm.dbl; 
		_gm_sum += _gm.dbl;
		ptr_emit_exp[i] += _gm.dbl;

		_gm.u64 &= ~((uint64_t) (obs_base != i) - 1UL);
		_gm_sum_e += _gm.dbl;
#endif	/* SSE3 */
	}

	*ptr_qual_exp_c += _gm_sum;
	*ptr_qual_exp_e += _gm_sum_e;
}

static inline void update_gamma_and_xi_1st(kmer_t **ptr_ref_kmers, 
		DOUBLE *ptr_gamma, DOUBLE *ptr_xi, DOUBLE inv_exp_trans, 
		kbits_t obs_kid, kbits_t obs_N_flag, DOUBLE *ptr_qual_exp_c, 
		DOUBLE *ptr_qual_exp_e, DOUBLE *ptr_emit_exp, 
		const char *conv_q, int kmer_len)
{
	kmer_t *k;

	register __m128d x_denomt = _mm_set1_pd(inv_exp_trans);
	register __m128d x_inf = _mm_set1_pd(INFINITY);
	x_denomt = _mm_and_pd(_mm_cmpneq_pd(x_inf, x_denomt), x_denomt);

	for (int i = 0; i < 4; i++){
		register int _isll2 = i << 2;

		if ((k = (ptr_ref_kmers[i])) == NULL) continue;
#ifdef __SSE3__
		DOUBLE *ptr_kmer_xi = k->xi;
		kbits_t kid = k->id;
		DOUBLE _gm = ptr_gamma[i] * inv_exp_trans;

		register __m128d x_op_a = _mm_load_pd(ptr_kmer_xi);
		register __m128d x_op_b = _mm_mul_pd(
				_mm_load_pd(&ptr_xi[_isll2]), x_denomt);

		x_op_a = _mm_add_pd(x_op_a, x_op_b);
		_mm_store_pd(ptr_kmer_xi, x_op_a);

		x_op_a = _mm_load_pd(ptr_kmer_xi + 2);
		x_op_b = _mm_mul_pd(_mm_load_pd(&ptr_xi[_isll2 + 2]), x_denomt);
		x_op_a = _mm_add_pd(x_op_a, x_op_b);
		_mm_store_pd(ptr_kmer_xi + 2, x_op_a);

		k->gamma += _gm; 
		for (int j = 0; j < kmer_len; j++){
			register int _jsll1 = j << 1;
			kbits_t bj = (kid >> _jsll1) & 3;

			kbits_t _nnmask = ((obs_N_flag >> _jsll1) & 1) - 1;
			kbits_t obs_base = (((obs_kid >> _jsll1) & 3) & _nnmask) |
				(~_nnmask & 4);

			int qual = (int) conv_q[j];
			ptr_qual_exp_c[qual] += _gm;
			ptr_emit_exp[(j << 2) + bj] += _gm;

			if (obs_base != bj) {
				ptr_qual_exp_e[qual] += _gm;
			}
		}
#endif /* SSE3 */
	}
}

#endif
