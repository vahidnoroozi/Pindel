#include "data.h"
#include "trans.h"

static const char* file_name = "data.c";

static int cmp_kmer_ptr_sfx(const void *ka, const void *kb)
{
	if ((*(kmer_t **) ka)->id < (*(kmer_t **) kb)->id) {
		return -1;
	}
	else {
		return 1;
	}
}


static int masked_kmer_id_comp(const void *a, const void *b, const void *pmask)
{
	register kbits_t mask = ((mask_t *) pmask)->pmask;
	int shift = ((mask_t *) pmask)->mask_shift;
	register kbits_t id_a = (*(kbits_t *) a);
	register kbits_t id_b = (*(kbits_t *) b);
	register kbits_t mid_a = id_a & mask, mid_b = id_b & mask;

	kbits_t id_diff = ((~mask & id_a) >> shift) - ((~mask & id_b) >> shift);

	/* set to all 1 if id_1 != id_b, all 0 otherwise */
	register kbits_t non_equal_mask = (kbits_t) (mid_a == mid_b) - 1UL;
	register kbits_t _masked_sign = (int) (((mid_a - mid_b) >> 63) << 31) + 1;
	return (non_equal_mask & _masked_sign) | (~non_equal_mask & id_diff);
}

int compute_kmer_hamming_nbhds(data *d)
{
	int dmax = d->opt->preconstructed_nbhd_hd;
	int num_perms = d->n_choose_k[d->opt->kmer_length][dmax];
	int n_kmers = d->kmer_dict->used;

	kbits_t *kmer_ids = malloc(n_kmers * sizeof(kbits_t));

	dictentry *de;
	register int i, n = 0;
	for (i = n_kmers, de = d->kmer_dict->table; i > 0; de++) {
		if (de->value == NULL) continue;
		i--;

		kmer_t *pk = (kmer_t *) de->value;
		pk->n_neighbors = 1;
		pk->alloc_nbhd_size = 2;
		pk->data = malloc(sizeof(kmer_t *) << 2);
		/* a kmer is always a neighbor of itself */
		*(kmer_t **) pk->data = pk;

		kmer_ids[n++] = pk->id;
	}

	mask_t pos_mask = {0, 0};
	for (i = 0; i < num_perms; i++) {
		INFO("Neighborhood construction iteration %2d/%2d...\n", i+1, num_perms);

		int shift = (i << 1);
		kbits_t pmask = d->position_masks[i];
		pos_mask.mask_shift = shift;
		pos_mask.pmask = pmask;

		/* sort masked kmer ids */
		qsort_r(kmer_ids, n_kmers, sizeof(kbits_t), 
				masked_kmer_id_comp, (void *) &pos_mask);

		kbits_t *base_loc = kmer_ids;

		/* identify neighborhood for all first kmers */
		for (int j = 0; j < n_kmers; j++) {
			kmer_t *pk = dict_find(d->kmer_dict, kbits_cast_to_ptr(kmer_ids[j]))->value;

			kbits_t kid = pk->id;
			register kbits_t masked_kid = pk->id & pmask;

			/* find the location of the kmer id in sorted list */
			kbits_t *id_loc = base_loc;
			for (; *id_loc != kid; id_loc++);
			base_loc = id_loc;
		
			/* identify the neighborhood by exploring kmers within locality */
			kbits_t *ptr_lb, *ptr_ub;
			kbits_t masked_base = (~pmask & kid) >> shift;
			for (ptr_lb = id_loc - masked_base; ptr_lb <= id_loc + (3UL - masked_base) && 
					masked_kid != (pmask & *(ptr_lb)); ptr_lb++);
			for (ptr_ub = ptr_lb + 1; masked_kid == (pmask & *(ptr_ub)); 
					ptr_ub++);

			int n_masked_nbhd_size = ptr_ub - ptr_lb - 1;
			if (n_masked_nbhd_size > 0) {
				int pnn = pk->n_neighbors;
				pk->n_neighbors += n_masked_nbhd_size;

				int size_nbhd = pk->n_neighbors;
				int realloc_size = 1 << pk->alloc_nbhd_size;
				if (realloc_size < size_nbhd) {
					for (; realloc_size < size_nbhd; realloc_size <<= 1, 
							pk->alloc_nbhd_size++);
					pk->data = realloc(pk->data, 
							sizeof(kmer_t *) << pk->alloc_nbhd_size);
				}

				kmer_t **nb_ptrs = (kmer_t **) pk->data;
				for (kbits_t *pid = ptr_lb; pid < ptr_ub; pid++) {
					register kbits_t _id = *pid;
					if (_id != kid) {
						nb_ptrs[pnn++] = dict_find(d->kmer_dict, 
								kbits_cast_to_ptr(_id))->value;
					}
				}
			}
		}
	}

	free(kmer_ids);
}

void kmer_id2seq(kbits_t kid, int size)
{
	int i;
	for (i = 0; i < size; i++){
		switch((kid >> (i * 2)) & 3){
			case 0:
				putchar('A');
				break;
			case 1:
				putchar('C');
				break;
			case 2:
				putchar('T');
				break;
			case 3:
				putchar('G');
				break;
		}
	}
}

void adjust_kmer_labeling(data *d)
{
	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) << 1;
	double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);

	dictentry *de = d->kmer_dict->table;
	int i = d->kmer_dict->used;
	for (; i > 0; de++) {
		if (de->value == NULL) continue;

		i--;

		kmer_t *pk = de->value;
		kbits_t kid = pk->id;

		// skip (k+1)-mers, or kmers that have been labeled already
		if ((kid >> 62) || pk->dirty) continue;

		uint32_t graph_nodes = 0;
		uint32_t uv_trans_flags = 0;
		int uv_kcov[8] = {0};
		int uv_nonuniq[8] = {0};
		int nU = 0, nV = 0; // the two disjoint sets of kmers in the bipartite graph
		kmer_t *lbl_cands[8] = {NULL};

		/* set U: */
		for (kbits_t base = 0; base < 4; ++base) {
			kbits_t _ku = (kid & ~3UL) | base;
			dictentry *_deu = dict_find(d->kmer_dict,
					kbits_cast_to_ptr(_ku));

			if (_deu != NULL) {
				++nU;
				kmer_t* _pku = _deu->value;
				lbl_cands[base] = _pku;

				graph_nodes |= (1 << base);
				uv_trans_flags |= (_pku->trans_flag & 15) << (base * 4);
				int n_nonzero = 0;
				for (int b = 0; b < _pku->n_ktrans; ++b) {
					n_nonzero += (_pku->exp_trans_p[b] > 0);
				}
				uv_nonuniq[base] = n_nonzero > 1;
				uv_kcov[base] = _pku->exp_count;

				/*
				uv_trans_flags |= (kmer_counter2flag((uint64_t) _deu->value, 
						&uv_kcov[base], &uv_nonuniq[base]) << (base * 4));
				*/
			}
		}

		/* set V: */
		for (kbits_t base = 0; base < 4; ++base) {
			kbits_t _kv = kmer_reverse_complement((kid >> 2) | (base << shift), 
					kmer_len);

			kbits_t cbase = complementary_base(base);

			dictentry *_dev = dict_find(d->kmer_dict,
					kbits_cast_to_ptr(_kv));

			if (_dev != NULL) {
				++nV;
				kmer_t* _pkv = _dev->value;
				lbl_cands[4 + cbase] = _pkv;

				graph_nodes |= (1 << (4 + cbase));
				uv_trans_flags |= (_pkv->trans_flag & 15) << ((4 + cbase) * 4);
				int n_nonzero = 0;
				for (int b = 0; b < _pkv->n_ktrans; ++b) {
					n_nonzero += (_pkv->exp_trans_p[b] > 0);
				}
				uv_nonuniq[4 + cbase] = n_nonzero > 1;
				uv_kcov[4 + cbase] = _pkv->exp_count;

				/*
				uv_trans_flags |= (kmer_counter2flag((uint64_t) _dev->value,
						&uv_kcov[4 + cbase], &uv_nonuniq[4 + cbase]) << ((4 + cbase) * 4));
				*/
			}
		}

		/* --- graph partitioning --- */
		while (graph_nodes) {
			// we arbitrarily call the first four bases/kmers as "forward"
			int n_fwd_nonuniq = 0, n_rev_nonuniq = 0;
			int fwd_cov = 0, rev_cov = 0;

			uint32_t visited = 0, pending = 0;

			int seed = 0;
			for (seed = 0; seed < 8; ++seed) if (graph_nodes & (1 << seed)) break;

			if (seed < 4) {
				fwd_cov += uv_kcov[seed];
				n_fwd_nonuniq += uv_nonuniq[seed];
			}
			else {
				rev_cov += uv_kcov[seed];
				n_rev_nonuniq += uv_nonuniq[seed];
			}

			visited |= (1 << seed);
			pending |= (complementary_trans_flag(((uv_trans_flags >> (4 * seed)) & 15)) 
				<< ((seed < 4) * 4)) & graph_nodes;

			do {
				for (seed = 0; seed < 8; ++seed) if (pending & (1 << seed)) break;
				pending &= ~(1UL << seed);
				visited |= (1 << seed);

				if (seed < 4) {
					fwd_cov += uv_kcov[seed];
					n_fwd_nonuniq += uv_nonuniq[seed];
				}
				else {
					rev_cov += uv_kcov[seed];
					n_rev_nonuniq += uv_nonuniq[seed];
				}

				pending |= (complementary_trans_flag(((uv_trans_flags >> (4 * seed)) & 15)) 
						<< ((seed < 4) * 4)) & graph_nodes & (~visited);
			} while(pending);

			/* visited now contains a connected component in the bipartite graph */
			int label = (n_fwd_nonuniq < n_rev_nonuniq);

			if (n_fwd_nonuniq == n_rev_nonuniq) {
				label = nU > nV;
			}

			if ((n_fwd_nonuniq == n_rev_nonuniq) && (nU == nV)) {
				label = fwd_cov > rev_cov;
			}

			/*
			if (n_fwd_nonuniq == n_rev_nonuniq) {
				label = fwd_cov > rev_cov;
			}
			*/

			for (int base = 0; base < 8; ++base) {
				if ((visited & (1 << base)) == 0) continue;

				int set_label = (base < 4)? label : 1-label;

				kmer_t *pkb = lbl_cands[base];
				pkb->label = set_label;
				pkb->dirty = 1;
			}

			graph_nodes &= ~visited;
		}
	}
}

void adjust_kmer_strandedness(data *d)
{
	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) << 1;
	kmer_t *adj_cands[8];

	dictentry *de = d->kmer_dict->table;
	int i = d->kmer_dict->used;
	for (; i > 0; de++) {
		if (de->value == NULL) continue;

		int n_adj = 0;

		i--;
		kmer_t *pk = de->value;
		kbits_t kid = pk->id;

		if ((kid >> 62) || pk->dirty) continue;

		adj_cands[n_adj++] = pk;

		/* number of forward/reverse non-unique transitions.
		 * note that forward/reverse is deemed relative to the kmer identified
		 * by `kid`. */
		int n_fwd_nonuniq = 0; 
		int n_rev_nonuniq = 0;

		/* downstream kmer transition flag */
		uint64_t dns_tf = pk->trans_flag & 15;
		/* coverage on forward/reverse strand */
		uint32_t fwd_cov = pk->exp_count, rev_cov = 0;
		/* kmers sharing the same k-1 suffix with kid */
		uint64_t ups_tf = 0;

		/* keep the downstream transition flag for the "seed" kmer, i.e.
		 * the first kmer seen in the pair */
		//kbits_t seed_dns_tf = dns_tf | alt_trans_flag;

		/* create this detailed balance using transition flag obtained from the
		 * other strand as well. but do NOT count supplemented transitions
		 * toward total transition counts (so that we can tell on which strand
		 * the errors orginate) */
		kbits_t tpak_3p = trans_flag2packed(dns_tf);
		int n_trans = tpak_3p & 7, j = 0;
		int n_nonzero= 0;
		for (int base = 0; base < pk->n_ktrans; ++base) {
			n_nonzero += (pk->exp_trans_p[base] > 0);
		}
		n_fwd_nonuniq += (n_nonzero > 1);

		//if (n_trans > 1) n_fwd_lt_thresh += _nlt;
		for (j = 0, tpak_3p >>= 3; j < n_trans; ++j, tpak_3p >>= 2) {
			kbits_t b3p = tpak_3p & 3;
			kbits_t dns_rc_kid = kmer_reverse_complement(
					kbits_suffix(kid) | (b3p << shift), kmer_len);

			kmer_t *k_rc_dns = dict_find(d->kmer_dict,
					kbits_cast_to_ptr(dns_rc_kid))->value;

			// complementary transition flag
			uint32_t comp_tflag = k_rc_dns->trans_flag & 15;
			comp_tflag = ((comp_tflag & 3) << 2) | (comp_tflag >> 2);

			ups_tf |= comp_tflag; 
		}

		kbits_t tf_jnt_dns = dns_tf;
		if (ups_tf > (1 << (kid & 3))) {
			kbits_t tpak_j5p = trans_flag2packed(ups_tf & ~(1U << (kid & 3)));
			n_trans = tpak_j5p & 7;
			for (j = 0, tpak_j5p >>= 3; j < n_trans; ++j, tpak_j5p >>= 2) {
				kbits_t b5p = tpak_j5p & 3;
				kbits_t ups_kid = (kid & ~3UL) | b5p;
				kmer_t *k_ups = dict_find(d->kmer_dict,
						kbits_cast_to_ptr(ups_kid))->value;
				fwd_cov += k_ups->exp_count;

				tf_jnt_dns |= k_ups->trans_flag & 15; 
				n_nonzero= 0;
				for (int base = 0; base < k_ups->n_ktrans; ++base) {
					n_nonzero += (k_ups->exp_trans_p[base] > 0);
				}
				n_fwd_nonuniq += (n_nonzero > 1);
				/*
				if ((trans_flag2packed(_tf) & 7) > 1) {
					n_fwd_lt_thresh += _nlt;
				}
				*/

				adj_cands[n_adj++] = k_ups;
			}
		}

		int n_ups = n_adj;

		kbits_t tpak_j3p = trans_flag2packed(tf_jnt_dns);
		n_trans = tpak_j3p & 7;
		for (j = 0, tpak_j3p >>= 3; j < n_trans; ++j, tpak_j3p >>= 2) {
			kbits_t b3p = tpak_j3p & 3;
			kbits_t dns_rc_kid = kmer_reverse_complement(
					kbits_suffix(kid) | (b3p << shift), kmer_len);

			kmer_t *k_dns_rc = dict_find(d->kmer_dict,
					kbits_cast_to_ptr(dns_rc_kid))->value;

			n_nonzero = 0;
			for (int base = 0; base < k_dns_rc->n_ktrans; ++base) {
				n_nonzero += (k_dns_rc->exp_trans_p[base] > 0);
			}
			rev_cov += k_dns_rc->exp_count;
			n_rev_nonuniq += (n_nonzero > 1); 

			adj_cands[n_adj++] = k_dns_rc; 
		}

		/* if there are more forward transitions than reverse transitions,
		 * it is likely that there is an error in the forward transition. 
		 * hence, put the upstream kmer in set 0 (treated as forward) */
		uint64_t label = (n_fwd_nonuniq < n_rev_nonuniq);
		if (n_fwd_nonuniq == n_rev_nonuniq) {
			label = n_ups > (n_adj - n_ups);
		}

		if (n_ups == (n_adj - n_ups)) {
			// check coverage on both strands and favors the strand with higher
			// coverage.
			label = fwd_cov > rev_cov;
		}
		//uint64_t label = (n_fwd_lt_thresh < n_rev_lt_thresh);

		for (j = 0; j < n_adj; ++j) {
			kmer_t* kj = adj_cands[j];

			int set_label = (j < n_ups)? label : 1-label;
			kj->label = set_label;
			kj->dirty = 1;

			//d->total_kmer_exp_counts += kj->exp_count;
		}
	}
}

/* --- precomputed tables --- */

/* Table that converts an (offset) Phred quality score to the corresponding 
 * error probability of base calling. 
 * For display, the quality score is usually added to an offset (33 by
 * default).
 */
int compute_phred_to_probability(data *d)
{
	int i; 
	d->phred_to_prob = malloc(sizeof(*(d->phred_to_prob)) * 
			(d->opt->max_qual_score + 1));

	if (d->phred_to_prob == NULL)
		return message(stderr, file_name, __LINE__, ERROR_MSG,
				MEMORY_ALLOCATION, NULL);

	for (i = 0; i <= d->opt->max_qual_score; i++)
			d->phred_to_prob[i] = POW(10., -((double)i/10));

	return NO_ERROR;
}	/* compute_phred_to_probablity */

/** 
 * Compute binomial coefficients, i.e. N choose K using Pascal triangle.
 * The largest N to compute is kmer_length + 1.
 */
int compute_binomial_coefficients(data *d)
{
	int n = d->opt->kmer_length + 1;
	int i, j;

	d->n_choose_k = malloc(n * sizeof *d->n_choose_k);
	if (d->n_choose_k == NULL)
		return message(stderr, file_name, __LINE__,
				ERROR_MSG, MEMORY_ALLOCATION, NULL);

	for (i = 0; i < n; i++){
		d->n_choose_k[i] = calloc(i+1, sizeof **d->n_choose_k);
		if (d->n_choose_k[i] == NULL)
			goto ERROR_RETURN;

		/* Pascal's triangle */
		for (j = 0; j <= i; j++){
			if (j == 0 || j == i)
				d->n_choose_k[i][j] = 1;
			else
				d->n_choose_k[i][j] = d->n_choose_k[i-1][j] +
					d->n_choose_k[i-1][j-1];
		}
	}

	return NO_ERROR;

ERROR_RETURN:
	for (n = 1; n < i; n++)
		free(d->n_choose_k[n]);
	free(d->n_choose_k);
	return message(stderr, file_name, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, NULL);
}	/* compute_binomial_coefficients */

/* Permute all (k choose d) combinations for d erroneous positions in any kmer,
 * d = 1, 2, ..., d_max.
 * 
 * All computed results are stored in two-dimenional array: d->position_perms.
 * The 1st dimension of this array corresponds to d = 1, 2, ..., d_max.
 * The 2nd dimension for given d, is of length [(k choose d) * d]. 
 * Every d consecutive positions store a combination of d positions.
 */
int compute_error_positions(data *d)
{
	int ne, i, j;	/* number of errors */
	int num_perms;
	char *ptr_slots;	/* pointer to current d slots */
	int kmer_len = d->opt->kmer_length;
	int dmax = d->opt->preconstructed_nbhd_hd;
	kbits_t mask;

	d->position_perms = malloc(dmax * sizeof *d->position_perms);
	d->position_masks = malloc(d->n_choose_k[kmer_len][dmax] * 
			sizeof *d->position_masks);

	if (d->position_perms == NULL || d->position_masks == NULL)
		return message(stderr, file_name, __LINE__, ERROR_MSG,
				MEMORY_ALLOCATION, NULL);

	for (ne = 1; ne <= dmax; ne++){
		num_perms = d->n_choose_k[kmer_len][ne];
		d->position_perms[ne-1] = malloc(num_perms * ne * 
				sizeof **d->position_perms);

		if (d->position_perms[ne-1] == NULL)
			goto ERROR_RETURN;

		ptr_slots = d->position_perms[ne-1];

		/* initialization */
		mask = KBITS_MASK_MAX;
		for (i = 0; i < ne; i++){
			ptr_slots[i] = i + 1;
			mask &= ~(3UL << (i << 1));
		}

		if (ne == dmax)
			d->position_masks[0] = mask;

		memcpy(ptr_slots + ne, ptr_slots, ne * sizeof(*ptr_slots));
		ptr_slots += ne;	/* move pointer to next permutation */

		/* In total we have (n choose k) permutations, among which the initial
		 * permutation (1, 2, ..., d) counts 1. 
		 * Following code permute the remaining (n choose k) - 1 combinations.
		 */
		for (i = 1; i < num_perms; i++){
			ptr_slots[ne-1]++;	/* increase the rightmost position by 1 */
			for (j = ne - 1; j >= 0; j--){
				if (ptr_slots[j] > kmer_len - (ne - 1 - j)){
					/* carray on */
					ptr_slots[j-1]++;
					ptr_slots[j]  = ptr_slots[j-1] + 1;
					/* check if we can reset the next radix */
					if (j < ne - 1 && ptr_slots[j+1] > (kmer_len-(ne-2-j))){
						ptr_slots[j+1] = ptr_slots[j] + 1;
					}
				}
				else
					break;
			}

			if (ne == dmax){
				mask = KBITS_MASK_MAX; 
				for (j = 0; j < ne; j++){
					mask &= ~(3UL << ((ptr_slots[j] - 1) << 1));
				}
				d->position_masks[i] = mask;
			}

			if (ptr_slots[0] > (kmer_len - ne + 1))
				break;

			if (i < num_perms - 1){
				memcpy(ptr_slots + ne, ptr_slots, ne * sizeof(*ptr_slots));
				ptr_slots += ne;
			}
		}
	}

	return NO_ERROR;

ERROR_RETURN:
	for (i = 0; i < ne - 1; i++)
		free(d->position_perms[i]);
	free(d->position_perms);

	return message(stderr, file_name, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, NULL);
}	/* compute_error_positions */

/* compute all (k choose d) x 4^d error patterns, 
 * for d = 1, 2, ..., dmax.
 *
 */
int compute_all_errors(data *d)
{
	int i, j, k;
	int offset, num_perms, num_err_combs, num_states;
	kbits_t error_mask, dmer_mask, pos_shift;
	int kmer_len = d->opt->kmer_length;
	int dmax = d->opt->preconstructed_nbhd_hd;

	num_perms = d->n_choose_k[kmer_len][dmax];
	num_err_combs = pow(4, dmax);
	num_states = num_err_combs * num_perms;

	d->error_perms = malloc(sizeof(*d->error_perms) * num_states);
	if (d->error_perms == NULL)
		return MEMORY_ALLOCATION;

	for (i = 0; i < num_perms; i++){
		offset = i * dmax;
		for (j = 0; j < num_err_combs; j++){
			dmer_mask = j;
			error_mask = 0;

			for (k = 0; k < dmax; k++){
				pos_shift = (d->position_perms[dmax-1][offset+k] - 1) << 1;
				error_mask |= ((dmer_mask & 3) << pos_shift);

				dmer_mask >>= 2;
			}

			d->error_perms[i*num_err_combs+j] = error_mask;
		}
	}

	return NO_ERROR;
}
