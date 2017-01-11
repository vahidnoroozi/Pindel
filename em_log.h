#ifndef __EM_LOG_SCALE_H__
#define __EM_LOG_SCALE_H__

#include "em.h"

static DOUBLE compute_alpha(DOUBLE *ptr_ups_alpha, DOUBLE *norm_trans_p,
		kbits_t nonzero_ups_path, DOUBLE emit_dens);
static DOUBLE compute_beta(DOUBLE *ptr_dns_beta, DOUBLE *ptr_trans_p, 
		kbits_t nonzero_paths, DOUBLE alpha, 
		DOUBLE *ptr_emit_p, DOUBLE *ptr_exp_trans); 
static void normalize_emit_prob_mask(DOUBLE *raw_ep,
		kbits_t nonzero_emits, DOUBLE *norm_ep);

/* the final stage of forward algorithm, as well as the first stage of the 
 * backward algorithm (combined to reuse the same loop), computed on log
 * scale.
 */
int_pair_t compute_loglik(int tmax, kbits_t obs_base,
		DOUBLE *ptr_emit_p_list, dict *kpair_dict, 
		int shift, kbits_t dmax, dict *kdict, kmer_t **list_kmers, 
		DOUBLE *nx_e_logcounts, argmax_t *nx_amax_e_logcnt,
		DOUBLE *e_logcounts, argmax_t *amax_e_logcnt,
		DOUBLE *base_tp)
{
	int i;
	dictentry *de;
	int km1_shift = shift - 2;
	int fwd_kidx = 0, bwd_kidx = 0;

	DOUBLE _nx_max_exp_logcnt = NAN;
	int _nx_argmax_id = 0;
	DOUBLE _max_exp_logcnt = NAN;
	int _argmax_id = 0;

	/* normalized transition probability */
	//static DOUBLE norm_trans_p[4] __attribute__ ((aligned (16)));

	for (i = kpair_dict->used, de = kpair_dict->table; i > 0; de++) {
		if (de->value == NULL) continue;

		i--;

		kmer_pair_t *kpair = (kmer_pair_t *) de->value;

		kbits_t *ptr_dns_hd = kpair->dns_hamming_dist;
		DOUBLE *ptr_ups_alpha = kpair->ups_alpha;
		DOUBLE *ptr_dns_st_trans_p = kpair->dns_st_trans_p;
		DOUBLE *ptr_dns_gamma = kpair->dns_gamma;
		DOUBLE *ptr_dns_beta = kpair->dns_beta;
		kmer_t **ptr_dns_ref_kmers = (kmer_t **) kpair->dns_ref_kmers;
		kmer_pair_t **ptr_up_links = kpair->links;
		kmer_t **ptr_ups_kmers = (kmer_t **) kpair->ups_ref_kmers;
		kbits_t dns_pflag = kpair->dns_path_flag;

		kbits_t emit_flag = dns_pflag | ((1 << obs_base) &
				(((obs_base >> 2) & 1UL) - 1UL));
		DOUBLE *ptr_emit_p = &ptr_emit_p_list[(emit_flag - 1) << 2];

		kbits_t k_sfx = kpair->common_substr;
		kbits_t sfx_kth_base = (kpair->common_substr >> km1_shift);
		kbits_t sfx_kth_base_sll2 = sfx_kth_base << 2;

		kbits_t _upspf = trans_flag2packed(kpair->ups_path_flag);
		for (int j = 0; j < 4; j++) {
			if (ptr_dns_hd[j] > dmax) continue;
			/* ----- forward iteration ---- */
			DOUBLE *ptr_trans_p = &ptr_dns_st_trans_p[j << 2];

			DOUBLE _alpha = compute_alpha(ptr_ups_alpha, ptr_trans_p, 
					_upspf, ptr_emit_p[j]);

			/* record the maximum expected count within the last 
			 * neighborhood. */
			if (isnan(_nx_max_exp_logcnt) || (!isnan(_alpha) &&
						_nx_max_exp_logcnt < _alpha)) {
				_nx_max_exp_logcnt = _alpha; 
				_nx_argmax_id = fwd_kidx;
			}
			nx_e_logcounts[fwd_kidx] = _alpha;

			kbits_t nb_kbits = k_sfx | ((kbits_t) j << shift);
			kmer_t *nb_kmer = (kmer_t *) (dict_find(kdict, 
						CAST_TO_PTR(nb_kbits))->value);
			
			ptr_dns_ref_kmers[j] = nb_kmer;
			ptr_dns_beta[j] = 0.0;
			ptr_dns_gamma[j] = _alpha;

			list_kmers[fwd_kidx] = nb_kmer;

			fwd_kidx++;
		}

		for (int j = 0; j < 4; j++) {
			/* ----- backward iteration ----- */
			kmer_pair_t *_up_link = ptr_up_links[j];
			if (_up_link == NULL) continue;

			DOUBLE *ptr_trans_p = ptr_ups_kmers[j]->transition_p;

			kbits_t _tran_nzp = trans_flag2packed(
				ptr_ups_kmers[j]->trans_flag);

			DOUBLE *ptr_exp_trans = &_up_link->exp_trans_p[sfx_kth_base_sll2];
			DOUBLE _beta = compute_beta(ptr_dns_beta, ptr_trans_p, 
					_tran_nzp, ptr_ups_alpha[j], ptr_emit_p, ptr_exp_trans);

			/* expected count of kmer on log scale */
			DOUBLE _logcnt = PRODUCT(ptr_ups_alpha[j], _beta);

			/* set up variables needed for next iteration and expectation 
			 * normalization */
			_up_link->dns_beta[sfx_kth_base] = _beta;
			_up_link->dns_gamma[sfx_kth_base] = _logcnt;

			e_logcounts[bwd_kidx] = _logcnt;

			/* find largest expected log count for follow-up normalization */
			if (isnan(_max_exp_logcnt) || (!isnan(_logcnt) && 
						_max_exp_logcnt < _logcnt)) {
				_max_exp_logcnt = _logcnt;
				_argmax_id = bwd_kidx;
			}

			bwd_kidx++;
		}
	}

	nx_amax_e_logcnt->maxval = _nx_max_exp_logcnt;
	nx_amax_e_logcnt->argmax = _nx_argmax_id; 
	amax_e_logcnt->maxval = _max_exp_logcnt;
	amax_e_logcnt->argmax = _argmax_id; 

	return (int_pair_t) {fwd_kidx, bwd_kidx};
}

int final_backward_iteration(dict *kpair_dict, DOUBLE *ptr_emit_p_list,
		DOUBLE inv_sum_e_logcnt, kbits_t obs_base, 
		DOUBLE *ptr_qual_exp_c, DOUBLE *ptr_qual_exp_e, 
		DOUBLE *ptr_emit_exp, kmer_t **list_kmers, DOUBLE *k_exp_logcnts, 
		DOUBLE *k_exp_logtrans, argmax_t *amax_e_logcnt, kbits_t dmax,
		DOUBLE *t_marg_prob, int update_exp, DOUBLE *base_tp)
{
	register int i, j;
	int kidx = 0;
	dictentry *de;

	register dbl_t _max_exp_logcnt = {.dbl = NAN};
	register int _argmax_id = 0;

	register DOUBLE sum_qec = 0.0, sum_qee = 0.0;

	DOUBLE t_err_rate = 0.0;
	for (i = kpair_dict->used, de = kpair_dict->table; i > 0; de++) {
		if (de->value == NULL) continue;

		i--;

		kmer_pair_t *kpair = (kmer_pair_t *) de->value;

		DOUBLE *ptr_dns_gamma = kpair->dns_gamma;

		DOUBLE *ptr_dns_beta = kpair->dns_beta;
		DOUBLE *ptr_ups_alpha = kpair->ups_alpha;
		kmer_t **ptr_ups_kmers = (kmer_t **) kpair->ups_ref_kmers;
		kmer_t **ptr_dns_kmers = (kmer_t **) kpair->dns_ref_kmers;

		kbits_t dns_pflag = kpair->dns_path_flag;

		kbits_t emit_flag = dns_pflag | ((1 << obs_base) &
				(((obs_base >> 2) & 1UL) - 1UL));
		DOUBLE *ptr_emit_p = &ptr_emit_p_list[(emit_flag - 1) << 2];

		/* (1) Normalize the expected kmer counts/transitions
		 * and update the summary expected values needed by maximization. */
		//kbits_t _dnspf = _bwd_trans_lookup_table[dns_pflag]; 
		kbits_t _dnspf = trans_flag2packed(dns_pflag);
		kbits_t nonzero_paths = _dnspf;
		for (j = nonzero_paths & 7, nonzero_paths >>= 3; j > 0; 
				j--, nonzero_paths >>= 2) {
			int dns_base = nonzero_paths & 3;

			kmer_t *ref_k = ptr_dns_kmers[dns_base];
			DOUBLE *ptr_xi = ref_k->xi;

			register DOUBLE _nlc = EXP(PRODUCT(ptr_dns_gamma[dns_base], 
							inv_sum_e_logcnt));
			register dbl_t _norm_logcnt = {.dbl = _nlc};

			/* only update this when dns_base equals obs_base */
			uint64_t _equal_mask = ((uint64_t) (dns_base != obs_base) - 1UL);
			t_marg_prob[dns_base] += _nlc;

			if (update_exp) {
				ref_k->gamma += _norm_logcnt.dbl;
				ptr_emit_exp[dns_base] += _norm_logcnt.dbl;

				_norm_logcnt.u64 &= _equal_mask; 
				sum_qec += _norm_logcnt.dbl;

				/* only update this when dns_base differs from obs_base */
				_norm_logcnt.dbl = _nlc;
				_norm_logcnt.u64 &= ~_equal_mask;
				sum_qee += _norm_logcnt.dbl;

				DOUBLE *ptr_exp_trans = &kpair->exp_trans_p[dns_base << 2];
				kbits_t nonzero_xi = trans_flag2packed(
						ref_k->trans_flag);
				/*
				kbits_t nonzero_xi = _bwd_trans_lookup_table[
					ref_k->trans_flag];
				*/
				register int k;
				for (k = nonzero_xi & 7, nonzero_xi >>= 3; k > 0; 
						k--, nonzero_xi >>= 2) {
					int xidx = nonzero_xi & 3;
					ptr_xi[xidx] += EXP(PRODUCT(ptr_exp_trans[xidx], 
								inv_sum_e_logcnt));
				}
			}
		}

		/* (2) Compute one iteration of backward algorithm, and prepare 
		 * values for normalization in next iteration. */
		kbits_t _upspf = trans_flag2packed(kpair->ups_path_flag);
		//kbits_t _upspf = _bwd_trans_lookup_table[kpair->ups_path_flag];
		nonzero_paths = _upspf;
		for (j = nonzero_paths & 7, nonzero_paths >>= 3; j > 0; 
				j--, nonzero_paths >>= 2) {
			int ups_base = nonzero_paths & 3;
			kmer_t *_ups_kmer = ptr_ups_kmers[ups_base];

			/* dns_pflag contains all possible states can be reached
			 * from upstream states, CONSIDERING THE DMAX CONSTRAINT;
			 * where as _ups_kmer->trans_flag indicates which downstream
			 * nucleotides _ups_kmer can transition into, WITHOUT DMAX
			 * CONSTRAINT!.
			 *
			 * Hence, as long as _ups_kme->trans_flag is a subset of 
			 * dns_pflag, we will not normalize the transition probability.
			 */
			kbits_t _tran_nzp = trans_flag2packed(_ups_kmer->trans_flag);

			DOUBLE *ptr_exp_trans = &k_exp_logtrans[kidx << 2]; 
			DOUBLE _beta = compute_beta(ptr_dns_beta, _ups_kmer->transition_p,
					_tran_nzp, ptr_ups_alpha[ups_base], 
					ptr_emit_p, ptr_exp_trans);
			dbl_t _logcnt = {.dbl = PRODUCT(ptr_ups_alpha[ups_base], _beta)};

			k_exp_logcnts[kidx] = _logcnt.dbl;
			list_kmers[kidx] = _ups_kmer;

			/* find largest expected log count for follow-up normalization */
			uint64_t ge_mask = (uint64_t) (isnan(_max_exp_logcnt.dbl) ||  
						(_max_exp_logcnt.dbl < _logcnt.dbl)) - 1UL;
			uint64_t inv_ge_mask = ~ge_mask;

			_max_exp_logcnt.u64 = (ge_mask & _max_exp_logcnt.u64) | 
				(inv_ge_mask & _logcnt.u64);
			_argmax_id = (ge_mask & _argmax_id) | (inv_ge_mask & kidx);

			kidx++;
		}
	}

	*ptr_qual_exp_c += sum_qec;
	*ptr_qual_exp_e += sum_qee;
	
	amax_e_logcnt->maxval = _max_exp_logcnt.dbl;
	amax_e_logcnt->argmax = _argmax_id;

	return kidx;
}

/* one iteration of backward algorithm computed on log scale */
int backward_iteration(dict *kpair_dict, DOUBLE *ptr_emit_p_list,
		DOUBLE inv_sum_e_logcnt, kbits_t obs_base, 
		DOUBLE *ptr_qual_exp_c, DOUBLE *ptr_qual_exp_e,
		DOUBLE *ptr_emit_exp, DOUBLE *k_exp_logcnts, 
		argmax_t *amax_e_logcnt, int km1_shift, kbits_t dmax, 
		DOUBLE *t_marg_prob, int update_exp, DOUBLE *base_tp)
{
	register int i, j;
	int kidx = 0;
	int avail_kpairs = kpair_dict->used;
	dictentry *de;

	register dbl_t _max_exp_logcnt = {.dbl = NAN};
	register int _argmax_id = 0;

	register DOUBLE sum_qec = 0.0, sum_qee = 0.0;

	DOUBLE t_err_rate = 0.0;
	for (i = avail_kpairs, de = kpair_dict->table; i > 0; de++) {
		if (de->value == NULL) continue;

		i--;

		kmer_pair_t *kpair = (kmer_pair_t *) de->value;

		DOUBLE *ptr_dns_gamma = kpair->dns_gamma;

		DOUBLE *ptr_dns_beta = kpair->dns_beta;
		DOUBLE *ptr_ups_alpha = kpair->ups_alpha;
		kmer_pair_t **ptr_up_links = kpair->links;
		kmer_t **ptr_ups_kmers = (kmer_t **) kpair->ups_ref_kmers;
		kmer_t **ptr_dns_kmers = (kmer_t **) kpair->dns_ref_kmers;

		kbits_t dns_pflag = kpair->dns_path_flag;
		kbits_t emit_flag = dns_pflag | ((1 << obs_base) &
				(((obs_base >> 2) & 1UL) - 1UL));
		DOUBLE *ptr_emit_p = &ptr_emit_p_list[(emit_flag - 1) << 2];

		kbits_t sfx_kth_base = (kpair->common_substr >> km1_shift);
		kbits_t sfx_kth_base_sll2 = sfx_kth_base << 2;

		/* (1) Normalize the expected kmer counts/transitions
		 * and update the summary expected values needed by maximization. */
		kbits_t _dnspf = trans_flag2packed(dns_pflag);
		kbits_t nonzero_paths = _dnspf;

		for (j = nonzero_paths & 7, nonzero_paths >>= 3; j > 0; 
				j--, nonzero_paths >>= 2) {
			int dns_base = nonzero_paths & 3;

			kmer_t *ref_k = ptr_dns_kmers[dns_base];
			DOUBLE *ptr_xi = ref_k->xi;

			register DOUBLE _nlc = EXP(PRODUCT(ptr_dns_gamma[dns_base], 
							inv_sum_e_logcnt));
			register dbl_t _norm_logcnt = {.dbl = _nlc};

			uint64_t _equal_mask = ((uint64_t) (dns_base != obs_base) - 1UL);

			t_marg_prob[dns_base] += _nlc;

			if (update_exp) {
				ref_k->gamma += _norm_logcnt.dbl;
				ptr_emit_exp[dns_base] += _norm_logcnt.dbl;

				/* only update this when dns_base equals obs_base */
				_norm_logcnt.u64 &= _equal_mask; 
				sum_qec += _norm_logcnt.dbl;

				/* only update this when dns_base differs from obs_base */
				_norm_logcnt.dbl = _nlc;
				_norm_logcnt.u64 &= ~_equal_mask; 
				sum_qee += _norm_logcnt.dbl;

				DOUBLE *ptr_exp_trans = &kpair->exp_trans_p[dns_base << 2];

				kbits_t nonzero_xi = trans_flag2packed(
					ref_k->trans_flag);

				register int k;
				for (k = nonzero_xi & 7, nonzero_xi >>= 3; k > 0; 
						k--, nonzero_xi >>= 2) {
					int xidx = nonzero_xi & 3;
					ptr_xi[xidx] += EXP(PRODUCT(ptr_exp_trans[xidx], 
								inv_sum_e_logcnt));
				}
			}
		}


		/* (2) Compute one iteration of backward algorithm, and prepare 
		 * values for normalization in next iteration. */
		kbits_t _upspf = trans_flag2packed(kpair->ups_path_flag);
		nonzero_paths = _upspf;
		for (j = nonzero_paths & 7, nonzero_paths >>= 3; j > 0; 
				j--, nonzero_paths >>= 2) {
			int ups_base = nonzero_paths & 3;
			kmer_t *_ups_kmer = ptr_ups_kmers[ups_base];

			kbits_t _tran_nzp = trans_flag2packed(_ups_kmer->trans_flag);

			kmer_pair_t *_up_link = ptr_up_links[ups_base];
			DOUBLE *ptr_exp_trans = &_up_link->exp_trans_p[sfx_kth_base_sll2];

			DOUBLE _beta = compute_beta(ptr_dns_beta, _ups_kmer->transition_p, 
					_tran_nzp, ptr_ups_alpha[ups_base], 
					ptr_emit_p, ptr_exp_trans);
			dbl_t _logcnt = {.dbl = PRODUCT(ptr_ups_alpha[ups_base], _beta)};

			/* set up variables needed for next iteration and expectation 
			 * normalization */
			_up_link->dns_beta[sfx_kth_base] = _beta;
			_up_link->dns_gamma[sfx_kth_base] = _logcnt.dbl;

			k_exp_logcnts[kidx] = _logcnt.dbl;

			/* find largest expected log count for follow-up normalization */
			uint64_t ge_mask = (uint64_t) (isnan(_max_exp_logcnt.dbl) ||  
						(_max_exp_logcnt.dbl < _logcnt.dbl)) - 1UL;
			uint64_t inv_ge_mask = ~ge_mask;

			_max_exp_logcnt.u64 = (ge_mask & _max_exp_logcnt.u64) | 
				(inv_ge_mask & _logcnt.u64);
			_argmax_id = (ge_mask & _argmax_id) | (inv_ge_mask & kidx);

			kidx++;
		}
	}

	*ptr_qual_exp_c += sum_qec;
	*ptr_qual_exp_e += sum_qee;
	
	amax_e_logcnt->maxval = _max_exp_logcnt.dbl;
	amax_e_logcnt->argmax = _argmax_id;

	return kidx;
}

/* This function is fully optimized to remove all branching conditions.
 * Notice that this function does NOT normalize transition probability
 * ON SITE when on log scale, because eventually we need an extra loop
 * for either computation or setting up transition probabilities 
 * downstream, we can do the the actual normalization there */
static inline void normalize_emit_prob_mask(DOUBLE *raw_ep,
		kbits_t nonzero_emits, DOUBLE *norm_ep)
{
	int i;

	memset(norm_ep, 0, sizeof(*norm_ep) * 5);

	register kbits_t nzef = nonzero_emits;
	DOUBLE sum_emit_dens = raw_ep[4];
	for (i = nzef & 7, nzef >>= 3; i > 0; i--, nzef >>= 2) {
		kbits_t obs_base = nzef & 3;
		norm_ep[obs_base] = raw_ep[obs_base];

		sum_emit_dens += norm_ep[obs_base];
	}

	/* copy N emission probability */
	norm_ep[4] = raw_ep[4] / sum_emit_dens;

	nzef = nonzero_emits;
	for (i = nzef & 7, nzef >>= 3; i > 0; i--, nzef >>= 2) {
		kbits_t obs_base = nzef & 3;
		norm_ep[obs_base] /= sum_emit_dens; 
	}
}

void kmer_add_to_pair(kmer_t *k, kbits_t kid, 
		kbits_t *hamming_dist_v, DOUBLE *norm_trans_p, 
		DOUBLE alph, kmer_pair_t *link_kpair, dict *kp_dict, 
		kbits_t dns_path_flag, kbits_t nonzero_trans, mempool_t *mp) 
{
	static kmer_pair_t _kpair_tpl = {{NAN, NAN, NAN, NAN}, /* ups_alpha*/
		{NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
			NAN, NAN, NAN, NAN}, /* dns_st_trans_p */
		{NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN,
			NAN, NAN, NAN, NAN}, /* exp_trans_p */
		{NULL, NULL, NULL, NULL}, /* dns_ref_kmers */
		{KBITS_MASK_MAX, KBITS_MASK_MAX, KBITS_MASK_MAX, KBITS_MASK_MAX},
		{NAN, NAN, NAN, NAN}, /* dns_beta */
		{NAN, NAN, NAN, NAN}, /* dns_gamma */
		{NULL, NULL, NULL, NULL}, /* ups_ref_kmers */
		{NULL, NULL, NULL, NULL}, /* links */
		0, 0, 0};

	dictentry *de;
	kmer_pair_t *kpair;
	kbits_t k_sfx = KBITS_SUFFIX(kid);

	register kbits_t ups_base = kid & 3;

	de = dict_add_find(kp_dict, CAST_TO_PTR(k_sfx), NULL);
	if (de->value == NULL){
		/* --- create a kmer pair --- */
		/* kmer_pair_t must be aligned by 16 bytes! */
		kpair = mempool_nalloc(mp, sizeof(*kpair), 16);
		memcpy(kpair, &_kpair_tpl, sizeof(*kpair));
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

	for (i = nonzero_trans & 7, nonzero_trans >>= 3; i > 0;
			i--, nonzero_trans >>= 2) {
		int dns_base = nonzero_trans & 3;
		kpair->dns_hamming_dist[dns_base] = hamming_dist_v[dns_base];

		_dns_trans_p[(dns_base << 2) + ups_base] = norm_trans_p[dns_base];
	}
}

static inline DOUBLE compute_alpha(DOUBLE *ptr_ups_alpha, DOUBLE *norm_trans_p,
		kbits_t nonzero_ups_path, DOUBLE emit_dens)
{
	static DOUBLE summands[4];

	register kbits_t nzpf = nonzero_ups_path;
	int n_nz_ups_path = nzpf & 7;

	/* if current kmer only has one non-zero transition/pathway upstream,
	 * we can easily compute alpha, granted that transition probability is
	 * 1.0 (or 0.0 on log-scale) */
	if (n_nz_ups_path == 1) {
		int ups_base = nonzero_ups_path >> 3;
		return PRODUCT(PRODUCT(ptr_ups_alpha[ups_base],
					norm_trans_p[ups_base]), emit_dens);
	}
	else {
		dbl_t max_summand = {.dbl = NAN};
		int argmax_idx = 0, i;
		for (i = n_nz_ups_path, nzpf >>= 3; i > 0; i--, nzpf >>= 2) {
			int ups_base = nzpf & 3;

			register dbl_t _summand = {.dbl = PRODUCT(ptr_ups_alpha[ups_base], 
					norm_trans_p[ups_base])};

			/* (greater or equal) ge_mask is set to 0, 
			 * if max_summand < _summand */
			uint64_t ge_mask = (uint64_t) (isnan(max_summand.dbl) || 
					(max_summand.dbl < _summand.dbl)) - 1UL;
			uint64_t inv_ge_mask = ~ge_mask;

			max_summand.u64 = (ge_mask & max_summand.u64) | 
				(inv_ge_mask & _summand.u64);
			argmax_idx = (ge_mask & argmax_idx) | (inv_ge_mask & ups_base);

			summands[ups_base] = _summand.dbl;
		}

		/* compute the sum of all summands gives us alpha */
		DOUBLE _expsum = 0.0;
		for (i = n_nz_ups_path, nzpf = nonzero_ups_path >> 3; i > 0; 
				i--, nzpf >>= 2) {
			int ups_base = nzpf & 3;

			register uint64_t argmax_mask = ((ups_base == argmax_idx) - 1UL);
			register dbl_t _summand = {
				.dbl = EXP(summands[ups_base] - max_summand.dbl)};
			_summand.u64 &= argmax_mask;

			_expsum += _summand.dbl; 
		}

		return PRODUCT(emit_dens,
				PRODUCT(max_summand.dbl, log1p(_expsum)));
	}
}

static inline DOUBLE compute_beta(DOUBLE *ptr_dns_beta, DOUBLE *ptr_trans_p, 
		kbits_t nonzero_paths, DOUBLE alpha,
		DOUBLE *ptr_emit_p, DOUBLE *ptr_exp_trans)
{
	static dbl_t tails[4];
	dbl_t max_tail = {.dbl = NAN};
	int argmax_tail = 0;

	register int i;
	register kbits_t nzpf = nonzero_paths;

	for (i = nzpf & 7, nzpf >>= 3; i > 0; i--, nzpf >>= 2) {
		int dns_base = nzpf & 3;

		tails[dns_base].dbl = PRODUCT(
				ptr_trans_p[dns_base],
				PRODUCT(ptr_emit_p[dns_base], ptr_dns_beta[dns_base]));

		/* if max_tail < tails[dns_base], ge_mask is set to 0, 
		 * otherwise UINT64_MAX */
		uint64_t ge_mask = (uint64_t) (isnan(max_tail.dbl) || 
				(max_tail.dbl < tails[dns_base].dbl)) - 1UL;
		uint64_t inv_ge_mask = ~ge_mask;

		max_tail.u64 = (ge_mask & max_tail.u64) | (inv_ge_mask & tails[dns_base].u64);
		argmax_tail = (argmax_tail & ge_mask) | (dns_base & inv_ge_mask);

		ptr_exp_trans[dns_base] = PRODUCT(alpha, tails[dns_base].dbl);
	}

	DOUBLE _expsum = 0.0; 
	for (i = nonzero_paths & 7, nzpf = nonzero_paths >> 3; i > 0; 
			i--, nzpf >>= 2) {
		int dns_base = nzpf & 3;

		register uint64_t argmax_mask = ((dns_base == argmax_tail) - 1UL);
		register dbl_t summand = {.dbl = EXP(tails[dns_base].dbl - max_tail.dbl)};
		summand.u64 &= argmax_mask;

		_expsum += summand.dbl; 
	}

	return PRODUCT(max_tail.dbl, log1p(_expsum));
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
		DOUBLE *norm_trans_p, kbits_t nonzero_ups_path, 
		DOUBLE emit_dens, DOUBLE *delta, kmer_pair_t *kpair)
{
	dbl_t max_delta = {.dbl = NAN};
	int argmax_delta = 0, i;
	
	int n_nz_ups_path = nonzero_ups_path & 7;
	if (n_nz_ups_path == 1) {
		argmax_delta = nonzero_ups_path >> 3;
		max_delta.dbl = PRODUCT(ptr_ups_alpha[argmax_delta],
				norm_trans_p[argmax_delta]);
	}
	else {
		for (i = nonzero_ups_path & 7, nonzero_ups_path >>= 3; i > 0;
				i--, nonzero_ups_path >>= 2) {
			int ups_base = nonzero_ups_path & 3;

			dbl_t _delta = {.dbl = PRODUCT(ptr_ups_alpha[ups_base], 
					norm_trans_p[ups_base])};

			/* is max_delta > _delta ? */
			uint64_t ge_mask = (uint64_t) (isnan(max_delta.dbl) ||
					(max_delta.dbl < _delta.dbl)) - 1UL;
			uint64_t inv_ge_mask = ~ge_mask;

			argmax_delta = (argmax_delta & ge_mask) | (ups_base & inv_ge_mask);
			max_delta.u64 = (max_delta.u64 & ge_mask) | (_delta.u64 & inv_ge_mask);
		}
	}

	*delta = PRODUCT(max_delta.dbl, emit_dens);

	return (kmer_pair_t *) ((uintptr_t) kpair + 
			(uintptr_t) argmax_delta);
}

static DOUBLE update_expected_values_t1(kmer_t **list_kmers, DOUBLE *e_logcnts,
		DOUBLE *e_logtrans, argmax_t *amax_e_logcnt, int n_kmers, 
		kbits_t obs_kid, kbits_t obs_N_flag, DOUBLE *ptr_qual_exp_c, 
		DOUBLE *ptr_qual_exp_e, DOUBLE *ptr_emit_exp, 
		const char *conv_q, int kmer_len)
{
	DOUBLE _expsum = 0.0;
	for (int i = 0; i < n_kmers; i++) {
		if (i == amax_e_logcnt->argmax) continue;

		_expsum += EXP(e_logcnts[i] - amax_e_logcnt->maxval);
	}
	DOUBLE inv_sum_e_logcnt = -PRODUCT(amax_e_logcnt->maxval,
			log1p(_expsum));

	/* update */
	for (int i = 0; i < n_kmers; i++) {
		register int _isll2 = i << 2;

		/* convert back to original scale before updating the log counts */
		DOUBLE _gm = EXP(PRODUCT(e_logcnts[i], inv_sum_e_logcnt));
		kmer_t *refk = list_kmers[i];
		kbits_t kid = refk->id;

		/* expected transitions */
		//kbits_t nonzero_xi = _ev_trans_lookup_table[refk->trans_flag];
		
		kbits_t nonzero_xi = trans_flag2packed(refk->trans_flag);
		register int k;
		for (k = nonzero_xi & 7, nonzero_xi >>= 3; k > 0; 
				k--, nonzero_xi >>= 2) {
			int xidx = nonzero_xi & 3;

			DOUBLE _elt = EXP(PRODUCT(e_logtrans[_isll2 + xidx],
					inv_sum_e_logcnt));
			refk->xi[xidx] += _elt;
		}

		refk->gamma += _gm;
		for (int j = 0; j < kmer_len; j++){
			register int _jsll1 = j << 1;
			register int _jsll2 = j << 2;

			kbits_t bj = (kid >> _jsll1) & 3;
			kbits_t _nnmask = ((obs_N_flag >> _jsll1) & 1) - 1;
			kbits_t obs_base = (_nnmask & ((obs_kid >> _jsll1) & 3)) |
				(~_nnmask & 4);

			int qual = (int) conv_q[j];
			ptr_emit_exp[_jsll2 + bj] += _gm;

			if (obs_base != bj) {
				ptr_qual_exp_e[qual] += _gm;
			}
			else {
				ptr_qual_exp_c[qual] += _gm;
			}
		}
	}

	return inv_sum_e_logcnt;
}

static void adjust_qual_score(char *qs, int pos, DOUBLE err_rate)
{
	int qual;
	if (err_rate < 1e-4) {
		qual = 73;
	}
	else {
		qual = (int) (-10 * LOG(err_rate) / LOG(10)) + 33;
	}

	qs[pos] = qual;
}

#endif
