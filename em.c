#include "em.h"
#include "read.h"
#include "mstep.h"
#include "kmer.h"
#include "iodata.h"
#include "trans.h"
#include "numeric.h"
#include "atomic.h"

#include <assert.h>
#include <stdbool.h>

//static int max_nbhd_size = 0;

/*
 * compute_per_read_nbhd_size()
 *
 * Computes kmer neighborhood size N_t, as well as heuristic neighborhood
 * Hamming distance constraint dmax_t,  for every position t, of any given 
 * read. 
 * */

static double triplet_del_exp[8];

static void hmm_init_model_params(data *d);
static void hmm_init_e_step(data *d);
static void hmm_fwd_trans_prob(double *tp, kbits_t ups_trans_flag, int ups_kmer_idx, 
		kmer_t **t_kmers, uint64_t *edist, double *ins_prop, kbits_t dns_bases, 
		kbits_t tnext_kid, double c_pi, double c_pd, double c_pn, int shift, int dmax, 
		dict *kdict, int unlimit_ins);
static void hmm_bwd_trans_prob(double *tp, kmer_t *pk, 
		double c_pi, double c_pi_self,  double c_pd, double c_pn, int shift, dict *kdict);
/*
static double hmm_compute_kmer_emission(data *d, kbits_t obs_kid, 
		kbits_t obs_n_flag, kmer_t *pk, int kmer_len, const char *conv_q);
*/
static double hmm_compute_alpha(int n_ups_trans,
	kbits_t ups_trans_flag, int ups_kmer_idx, double *t_state_alphas, 
	double *summands, double emit_dens);
static double hmm_compute_beta(kbits_t dns_trans_packed,
		int dns_kmer_idx, int pfx_idx, double *tnext_state_betas, double *kmer_trans_p,
		double *kmer_exp_trans, kbits_t k_unc_tf, kbits_t k_dmax_tf,
		int *tnext_sfx_order, double *norm_emit_dens, double alpha,
		kbits_t kid, kbits_t *tnext_states);
inline static int hmm_compute_delta(int n_ups_trans,
		kbits_t ups_trans_flag, int ups_kmer_idx, double *t_state_deltas,
		double *tp, double emit_dens, double *delta);
static int hmm_load_preconstructed_nbhd(data *d, int rid,
		kbits_t obs_1st_kid, kbits_t obs_n_flag, state_nbhd_t *first_nbhd, 
		int kmer_len, int compute_emit_dens, int observed_only,
		dict *kdict, mempool_t *mp);
static void hmm_setup_nbhd_ptrs(state_nbhd_t *nbhd, int n, 
		char *pmem);
static int hmm_count_distinct_suffixes(state_nbhd_t *nbhd, 
		int *tnext_nbhd_size, kbits_t t_dmax, 
		kbits_t obs_n_flag, kbits_t obs_next_base, double *unique_trans_p,
		int kmer_len, int unlimit_ins);
static double hmm_compute_emit_dens(
		data *d, kbits_t st, kbits_t xt, int yt, int yctx);
static void hmm_update_expected_vals(data *d, state_nbhd_t *curr_nbhd, 
		int t, kbits_t xt, int yt, double *t_exp_counts, double *t_exp_trans, 
		double inv_sum_exp_count, int kmer_len);
static int *compute_qual_contexts(const char *conv_q, int rlen, 
		int klen, int ctxlen, int qcnt, mempool_t *mp);
static void hmm_m_step(data *d);

void hmm_e_step(data *d, read_t *read, void *fdata);
void hmm_viterbi(data *d, read_t *read, void *fdata);
static void hmm_correct_errors(data *d, read_t *read, int tmax,
		kbits_t st, kbits_t xt, kbits_t xt_n_flag, int kmer_len);
static void hmm_correct_one_error(data *d, read_t *read, int t,
		kbits_t st, kbits_t xt);

static int kmer_cmp_suffix(const void *a, const void *b);

void EM(data *d, struct timeval *ts)
{
	double pll = 0.0, delta = 0.0;
	int iter = 0;

	/* initialize HMM parameters */
	hmm_init_model_params(d);

	if (d->opt->skip_em == 0)
		INFO("Running iterative Baum-Welch algorithms...\n");

	do {
		if (d->opt->skip_em) break;

		d->loglik = 0.0;
		d->n_reads_zero_likelihood = 0;

		memset(triplet_del_exp, 0, 8 * sizeof(*triplet_del_exp));

		hmm_init_e_step(d);
		INFO("E-step / iteration: %2d\n", iter);
		io_iterate_reads_fastq(d, d->opt->num_threads, true, hmm_e_step_wrapper, NULL);

		INFO("Reads that have zero likelihood: %d\n",
				d->n_reads_zero_likelihood);

		INFO("M-step / iteration: %2d\n", iter);
		hmm_m_step(d);

		delta = d->loglik - pll;

		INFO("[EM/ITER %2d] LL: " DOUBLE_E_FMT " PLL: " DOUBLE_E_FMT "("
				DOUBLE_E_FMT ")\n", 
				iter++, d->loglik, pll, delta);

		pll = d->loglik;
	} while (FABS(delta / d->loglik) > d->opt->tol);

	INFO("Viterbi decoding.\n");

	d->n_reads_cannot_decode = 0;
	io_iterate_reads_fastq(d, d->opt->num_threads, true, hmm_viterbi_wrapper, NULL);

	INFO("Reads that cannot be decoded: %d\n.",
			d->n_reads_cannot_decode);
}

/***
 * Maximization step in EM algorithm to update model parameters.
 */
void hmm_m_step(data *d)
{
	d->ins_error_rate = LOG(d->ins_sum_exp_trans / d->all_sum_exp_trans);
	d->del_error_rate = LOG(d->del_sum_exp_trans / d->all_sum_exp_trans);
	//d->ins_error_rate_self = LOG(d->ins_sum_exp_trans_self / d->all_sum_exp_trans);
	d->error_free_p = LOG((d->all_sum_exp_trans - d->ins_sum_exp_trans - 
				d->del_sum_exp_trans) / d->all_sum_exp_trans);

	for (int si = 0; si < 4; ++si) {
		double _sum = 0.0;
		for (int xi = 0; xi < 5; ++xi) {
			_sum += d->base_emit_exp[si * 5 + xi];
		}

		for (int xi = 0; xi < 5; ++xi) {
			d->base_emission_p[si * 5 + xi] = 
				LOG(d->base_emit_exp[si * 5 + xi] / _sum);
		}
	}

	int qmax = d->opt->max_qual_score + 1;
	for (int qctx = 0; qctx < 4; ++qctx) {
		double _sum = 0.0;
		for (int q = 0; q < qmax; ++q) {
			_sum += d->qual_emit_exp[qctx * qmax + q];
		}
		for (int q = 0; q < qmax; ++q) {
			d->qual_emission_p[qctx * qmax + q] = 
				LOG(d->qual_emit_exp[qctx * qmax + q] / _sum);
		}
	}

	INFO("Estimated errors rates: pi=%lg / pd=%lg\n",
			d->ins_sum_exp_trans / d->all_sum_exp_trans, 
			d->del_sum_exp_trans / d->all_sum_exp_trans);

	mstep_transition(d);
	mstep_initial(d);
	//mstep_emission(d);
}

void hmm_e_step_wrapper(data *d, read_t *read, void *fdata)
{
	int rerun;
	do {
		rerun = false;
		// rerun E-step on the same read if the rerun flag is set.
		hmm_e_step(d, read, (void *) &rerun);
	} while(rerun);
}

/**:
 * Some general algorithm in forward-backward algorithm.
 *
 */
void hmm_e_step(data *d, read_t *read, void *fdata)
{
	static char *BASES[] = {"A", "C", "T", "G",
		"AA", "CA", "TA", "GA",
		"AC", "CC", "TC", "GC",
		"AT", "CT", "TT", "GT",
		"AG", "CG", "TG", "GG",
		"self"};


	/*
	int actg_trans_count[4];
	int cum_actg_trans_count[4];
	*/
	kbits_t tnext_kid_list[4];
	kmer_t* tnext_kmer_list[4];

	const int kmer_len = d->opt->kmer_length;
	const int num_kmers = read->length - kmer_len + 1;
	const kbits_t klen_flag = (1UL << (kmer_len << 1)) - 1;
	const int shift = (kmer_len - 1) << 1;
	const int qmax = d->opt->max_qual_score + 1;
	const int qoff = d->opt->qual_score_offset;
	const double c_pi = d->ins_error_rate;
	const double c_pi_self = d->ins_error_rate_self;
	const double c_pd = d->del_error_rate;
	/* 1 - c_pi - c_pd */
	const double c_pn = d->error_free_p; 

	int *rerun = fdata;

	if (read->length <= kmer_len) return;

	mempool_t *mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

	dict *kdict = d->kmer_dict;
	const kbits_t dmax = d->opt->max_edit_dist;

	/* sequence */
	char *rseq = read->sequence;

	state_nbhd_t *T_nbhds = mempool_nalloc(mp,
			sizeof(*T_nbhds) * num_kmers, 16);

	/* flag for N (not determined) base, first kmer only */
	kbits_t obs_n_flag = kmer_n_base_flag(rseq, kmer_len);
	kbits_t next_obs_n_flag = kmer_n_base_flag(rseq + 1, kmer_len);

	double tnext_max_alpha = NAN;
	int tnext_argmax_alpha = 0;

	/* convert string literal to numeric representation */
	/*
	kbits_t observed_kmers[BITS_TO_U64(read->length)];
	read_seq_to_numeric(rseq, observed_kmers, read->length, kmer_len);
	*/

	/* set up first kmer neighborhood */
	
	int tmin = d->opt->adjust_read_start_pos & 0x1 ? d->read_start_pos[read->id] : 0;
	T_nbhds[tmin].size = hmm_load_preconstructed_nbhd(
			d, read->id, kmer_seq_to_id(rseq + tmin, kmer_len, 0),
			obs_n_flag, &T_nbhds[tmin], kmer_len, true, 
			tmin > 0,	// if the starting position is shifted toward the 3' end,
						// use the observed kmer singleton as the neighborhood.
			kdict, mp);
	
	/*
	kbits_t obs_kid = observed_kmers[0] & ((1UL << (kmer_len << 1)) - 1);
	*/
	kbits_t obs_base = kmer_effective_base(rseq[tmin+kmer_len-1]);
	kbits_t obs_next_base = kmer_effective_base(rseq[tmin+kmer_len]);

	/* initial state distribution */
#ifdef PMR_RABINER_SCALING
	double rabiner_loglik = 0.0;
#endif

	int tmax = read->length - kmer_len;
	for (int t = tmin; t < tmax; t++) {
		int tnext_kp = t + kmer_len;
		int tnext_nbhd_size;

		obs_next_base = kmer_effective_base(rseq[t + kmer_len]);
		/*
		obs_kid = bitarray_get(observed_kmers, t << 1, kmer_len << 1);
		register kbits_t _new_nf = ((obs_next_base >> 2) << 1) | 
			(obs_next_base >> 2);
		next_obs_n_flag = (next_obs_n_flag >> 2) | (_new_nf << shift);
		*/

		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];


		int t_nbhd_size = curr_nbhd->size;

		kbits_t *t_states_sorted_sfx = T_nbhds[t].states_sorted_sfx;
		double *t_states_alphas = T_nbhds[t].alphas;
		double **t_kmer_trans_p = T_nbhds[t].kmer_trans_prob;

		/* ----- II. compute one iteration of forward algorithm ----- */

		/* --- II a. set up kmer neighborhood at next position --- */
		/* one linear scan to determine # of distinct suffix, 
		 * compute downstream actg nucleotide counts etc. */
		int n_distinct_suffixes = hmm_count_distinct_suffixes(curr_nbhd, 
				&tnext_nbhd_size, dmax, obs_n_flag,
				obs_next_base, d->transition_p, kmer_len,
				d->opt->unlimited_insertions & INS_MODE_EM);

		curr_nbhd->n_uniq_sfx = n_distinct_suffixes;

		if (tnext_nbhd_size == 0) {
			if ((d->opt->adjust_read_start_pos & 0x1) && tmin < tmax - 1) {
				// try to shift the starting position to the right...
				double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);
				int valid_kmer_located = false;
				int s = tmin + 1;
				for (; s < tmax; ++s) {
					kbits_t _kid = kmer_seq_to_id(rseq + s, kmer_len, 0);
					kmer_t *_pks = dict_find(kdict, kbits_cast_to_ptr(_kid))->value;
					if (_pks->count > thresh) {
						valid_kmer_located = true;
						break;
					}
				}

				if (valid_kmer_located) {
					*rerun = true;
					d->read_start_pos[read->id] = s;
				}
			}
		
			if (!*rerun) 
				atomic_add_u32(&d->n_reads_zero_likelihood, 1);
			// ++d->n_reads_zero_likelihood;
			goto imm_return;
		}

		/* allocate memory given tnext_nbhd_size */
		char *pmem = mempool_nalloc(mp,
				nbhd_alloc_size(tnext_nbhd_size), 16);
		candidate_t *tnext_candidates = mempool_alloc(mp,
			(tnext_nbhd_size + 1) * sizeof(*tnext_candidates));
		double *tnext_alphas = mempool_alloc(mp,
			tnext_nbhd_size * sizeof(double));

		next_nbhd->size = tnext_nbhd_size;
		hmm_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);

		//printf("[t=%2d / size:%d]\n", t+2, tnext_nbhd_size);

		/* FIXME: is there any means to wrap following code block into a 
		 * macro template? */

		int kmer_idx = 0, index_pfx = 0;

		tnext_max_alpha = NAN;
		tnext_argmax_alpha = 0;

		register kbits_t mismatch_flag = ~(1UL << obs_next_base) |
			~((obs_next_base >> 2) - 1UL);

		for (int i = 0; i < n_distinct_suffixes; i++) {
			/*
			int n_common_sfx_kmers = bitarray_get_pwr2(
					curr_nbhd->ba_distinct_sfx, i, 1) + 1;
			*/
			int n_common_sfx_kmers = 0;

			kbits_t _repr_kid = t_states_sorted_sfx[kmer_idx];
			//kbits_t mut_flag = ((_repr_kid ^ obs_kid) | obs_n_flag) & 3UL;
			/*
			kbits_t sfx_hd = bitarray_get_pwr2(curr_nbhd->ba_hamming_dist, 
					i, 2);
			*/
			kbits_t common_sfx = kbits_trans_suffix(_repr_kid);

			/* v 4 bits                   v 4 bits (higher)
			 *
			 * A  -                       - A
			 * C  |  -- k-1 substring --  | C
			 * T  |                       | T
			 * G  -                       - G
			 */
			kbits_t _tf = curr_nbhd->ba_pfx_sfx_flags[i];
			kbits_t ups_trans_flag = _tf & UINT32_MAX;
			/* which kmers/nucleotides current suffix transition into */
			kbits_t dns_trans_flag = _tf >> 32;

			int n_dns_trans = _mm_popcnt_u64(dns_trans_flag);
			int n_ups_trans = _mm_popcnt_u64(ups_trans_flag);

			/*
			printf("[T] ");
			kmer_id_to_seq(_repr_kid, kmer_len);
			printf(" --> ");
			*/

			uint64_t k_sfx_edist = UINT64_MAX;
			int k_sfx_n_err = kmer_len + 1;
			for (int j = 0; j < n_ups_trans; j++) {
				kbits_t _id = t_states_sorted_sfx[kmer_idx + j];
				uint64_t _edist = kbits_suffix_edist(_id,
						curr_nbhd->edit_distances[kmer_idx + j]);
				int _ne = _mm_popcnt_u64(_edist);
				if (_ne < k_sfx_n_err) {
					k_sfx_n_err = _ne;
					k_sfx_edist = _edist;
				}
			}

			for (int j = 0; j < 21; j++) {
				if (((dns_trans_flag >> j) & 1) == 0) continue;

				//printf(" %s ", BASES[j]);
				
				kbits_t next_bases = (j < 4) ? j : j - 4;
				int len = 0;

				if (j < 20) {
					tnext_kid_list[0] = (common_sfx & KBITS_KLEN_REMOVAL_FLAG) | 
						(next_bases << shift) | ((uint64_t) (j > 3) << 62);
					tnext_kmer_list[0] = dict_find(kdict, 
						kbits_cast_to_ptr(tnext_kid_list[0]))->value;
					
					len = 1;
				}
				else {
					len = 0;
					kbits_t _ktp = trans_flag2packed(ups_trans_flag & 15);
					int _nk = _ktp & 7;
					//_ktp >>= 3;
					for (int k = 0; k < _nk; k++, _ktp >>= 2) {
						//kbits_t ups_base = _ktp & 3;

						int _idx = kmer_idx + k;
						int _edist = _mm_popcnt_u64(
								curr_nbhd->edit_distances[_idx]);
						// allow the insertions to extend, regardless of the
						// current edit distance, under the "unlimited
						// insertions" mode.
						if (_edist < dmax || (
									(d->opt->unlimited_insertions & INS_MODE_EM) &&
									curr_nbhd->ins_prop[_idx] > 0.0)) {
							tnext_kid_list[len] = t_states_sorted_sfx[_idx];
							tnext_kmer_list[len] = curr_nbhd->kmer_ptrs[_idx];
							++len;
						}
					}
				}

				kbits_t st = (j == 20) ? kbits_last_base(_repr_kid, kmer_len) : 
					((j < 4) ? j : ((j-4) >> 2));

				int yctx = (j == 20) ? QCTX_SELF_TRANS :
					((j < 4) ? QCTX_KMER_EQUAL : QCTX_KP1MER);
			
				double emit_dens = hmm_compute_emit_dens(
						d, st, obs_next_base, 
						read->qscore[tnext_kp] - qoff,
						yctx);
				
				for (int k = 0; k < len; k++) {
					/* for every kmer only at current position, add itself
					to the candidate list at the next position */
					kbits_t tnext_kid = tnext_kid_list[k];
					kbits_t tnext_kmer = tnext_kmer_list[k];

					/*
					printf("NEXT(%d): ", j);
					kmer_id_to_seq(tnext_kid, kmer_len);
					putchar('\n');
					*/

					/* transition probability to j */
					double tp2j[n_ups_trans];
					hmm_fwd_trans_prob(&tp2j[0], ups_trans_flag, 
							kmer_idx, curr_nbhd->kmer_ptrs, 
							curr_nbhd->edit_distances,
							curr_nbhd->ins_prop,
							j, tnext_kid, c_pi, c_pd, c_pn, 
							shift, dmax, kdict, 
							d->opt->unlimited_insertions & INS_MODE_EM);

					double alpha = hmm_compute_alpha(
						n_ups_trans,
						ups_trans_flag, 
						kmer_idx,
						t_states_alphas,
						tp2j,
						emit_dens	/* emit probability is 1 */
						);

					//printf(">> alpha: %lg\n", EXP(alpha));

					/* temporarily sorted by prefix */
					tnext_candidates[index_pfx].id = tnext_kid;
					tnext_candidates[index_pfx].kmer_ptr = tnext_kmer;
					tnext_candidates[index_pfx].pfx_order = index_pfx;
					tnext_candidates[index_pfx].edist = k_sfx_edist | 
						((uint64_t) (j != obs_next_base) << (kmer_len-1));
					tnext_candidates[index_pfx].is_insertion = (j == 20);
					tnext_candidates[index_pfx].emit_dens = emit_dens;

					tnext_alphas[index_pfx++] = alpha;
				}
			}

			//printf("\n");

			kmer_idx += n_ups_trans;
		}

		qsort(tnext_candidates, tnext_nbhd_size, sizeof(candidate_t), 
			kmer_cmp_suffix);

		/* terminator */
		tnext_candidates[tnext_nbhd_size].id = UINT64_MAX;

		/*
		printf("[t=%d] %.*s -> %c\nNext candidates [%d]:\n", 
				t, kmer_len, rseq+t, rseq[t+kmer_len], tnext_nbhd_size);
		*/

		/* reduce the potentially redundant neighborhood */
		int reduced_nbhd_size = 0;
		int prev_i = 0;
		kbits_t prev_cand_id = (tnext_nbhd_size == 0) ? UINT64_MAX : tnext_candidates[0].id;
		candidate_t *prev_cand = NULL;
		double _alpha;
		double _ins_prop;
		uint64_t _min_edist = UINT64_MAX;
#ifdef PMR_RABINER_SCALING
		double ct = 0.0; //NAN;
#endif

		for (int i = 0; i <= tnext_nbhd_size; i++) {
			candidate_t *cand = &tnext_candidates[i];
			kbits_t cand_id = cand->id;

			/*
			if (i < tnext_nbhd_size) {
				printf("CANDIDATE %d ", i);
				kmer_id_to_seq(cand_id, kmer_len);
				printf(" %lg\n", EXP(tnext_alphas[cand->pfx_order]));
			}
			*/

			if (cand_id != prev_cand_id) {
				//kmer_id_to_seq(prev_cand->id, kmer_len);
				if (i > (prev_i + 1)) {
					/* to reduce redundant hidden states */
					_alpha = NAN; 
					double _ins_alpha = NAN;
					_min_edist = UINT64_MAX;
					for (int j = prev_i; j < i; ++j) {
						candidate_t *cj = &tnext_candidates[j];	
						next_nbhd->suffixes_order[cj->pfx_order] = 
							reduced_nbhd_size;
						double _alp = tnext_alphas[cj->pfx_order];
						_alpha = log_add(_alpha, _alp);
						if (cj->is_insertion)
							_ins_alpha = log_add(_ins_alpha, _alp); 

						if (cj->edist < _min_edist) {
							_min_edist = cj->edist;
						}
					}

					_ins_prop = EXP(_ins_alpha - _alpha);
					//printf(" -[R]-> alpha: %lg\n", EXP(_alpha));
				}
				else {
					/* no redundant hidden states */
					_ins_prop = (double) prev_cand->is_insertion;
					_min_edist = prev_cand->edist;
					_alpha = tnext_alphas[prev_cand->pfx_order];
					next_nbhd->suffixes_order[prev_cand->pfx_order] = 
						reduced_nbhd_size;
					//printf(" -[-]-> alpha: %lg\n", EXP(_alpha));
				}

#ifdef PMR_RABINER_SCALING
				ct += EXP(_alpha);
#endif
				//ct = log_add(ct, _alpha);

				next_nbhd->alphas[reduced_nbhd_size] = _alpha;
				next_nbhd->states_sorted_sfx[reduced_nbhd_size] = prev_cand->id;
				next_nbhd->kmer_ptrs[reduced_nbhd_size] = prev_cand->kmer_ptr;
				next_nbhd->emit_dens[reduced_nbhd_size] = prev_cand->emit_dens;
				next_nbhd->edit_distances[reduced_nbhd_size] = _min_edist;
				next_nbhd->ins_prop[reduced_nbhd_size] = _ins_prop; 

				bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
					reduced_nbhd_size, prev_cand->kmer_ptr->trans_flag, 5);

				update_max_quantity(_alpha, reduced_nbhd_size, 
						tnext_max_alpha, tnext_argmax_alpha);

				/*
				printf("%d ", reduced_nbhd_size);
				kmer_id_to_seq(cand_id, kmer_len);
				printf("\n");
				*/

				++reduced_nbhd_size;
				prev_i = i;
			}

			prev_cand = cand;
			prev_cand_id = cand_id;
		}

		next_nbhd->size = reduced_nbhd_size;

#ifdef PMR_RABINER_SCALING
		rabiner_loglik += -LOG(1/ct);

		/* rabiner scaling */
		for (int i = 0; i < reduced_nbhd_size; i++) {
			//kmer_id_to_seq(next_nbhd->states_sorted_sfx[i], kmer_len);
			next_nbhd->alphas[i] = next_nbhd->alphas[i] - LOG(ct);
			//printf(" %.15lg\n", EXP(next_nbhd->alphas[i]));
		}
#endif

		obs_base = obs_next_base;
		obs_n_flag = next_obs_n_flag;
		//getchar();
	}
	//getchar();

	/* compute log-likelihood at the last position */
	double hmm_loglik = log_sum(T_nbhds[tmax].size, T_nbhds[tmax].alphas,
			tnext_argmax_alpha, tnext_max_alpha);

#ifdef PMR_RABINER_SCALING
	//printf("%d %lg %lg\n", read->id, hmm_loglik, rabiner_loglik);
	//getchar();
#endif

	if (!isnan(hmm_loglik)) {
		atomic_add_dbl(&d->loglik, hmm_loglik);
		//d->loglik += hmm_loglik;

		int tmax_kp = tmax + kmer_len - 1;
		state_nbhd_t *tmax_nbhd = &T_nbhds[tmax];
		int tmax_nbhd_size = tmax_nbhd->size;

		/*
		double *tmax_base_emit_exp = &d->base_emit_exp[tmax_kp * 
			d->base_ctx_radix];
		double *tmax_qual_emit_exp = &d->qual_emit_exp[tmax_kp * 
			d->qual_ctx_radix];
		*/

		/* --- update emission expected values --- */
		kbits_t xt = kmer_effective_base(rseq[tmax_kp]);
		kbits_t yt = read->qscore[tmax_kp] - qoff;

		for (int i = 0; i < tmax_nbhd_size; i++) {
			double norm_exp_count = EXP(PRODUCT(tmax_nbhd->alphas[i],
					-hmm_loglik));

			kbits_t kid = tmax_nbhd->states_sorted_sfx[i];
			kbits_t st = kbits_last_base(kid, kmer_len);

			atomic_add_dbl(&d->base_emit_exp[st * 5 + xt], norm_exp_count);
			// d->base_emit_exp[st * 5 + xt] += norm_exp_count;

			if (kid >> 62) {
				atomic_add_dbl(&d->qual_emit_exp[qmax * QCTX_KP1MER + yt], 
						norm_exp_count);
				//d->qual_emit_exp[qmax * QCTX_KP1MER + yt] += norm_exp_count;
			}
			else {
				double ins_prop = tmax_nbhd->ins_prop[i];
				// if the i-th kmer is purely an insertion kmer,
				// ins_prop, the contribution from the insertion copy, should
				// be 1.0. otherwise it is between 0 and 1.
				atomic_add_dbl(&d->qual_emit_exp[qmax * QCTX_SELF_TRANS + yt], 
						norm_exp_count * ins_prop);
				atomic_add_dbl(&d->qual_emit_exp[qmax * (st != xt) + yt], 
						norm_exp_count * (1.0 - ins_prop));
				/*
				d->qual_emit_exp[qmax * QCTX_SELF_TRANS + yt] += 
					norm_exp_count * ins_prop;
				d->qual_emit_exp[qmax * (st != xt) + yt] += 
					norm_exp_count * (1-ins_prop);
				*/
			}
		}
		/*
		int yt = (int) conv_q[tmax_kp];

		for (int i = 0; i < tmax_nbhd_size; i++) {
			double norm_exp_count = EXP(PRODUCT(tmax_nbhd->alphas[i],
					-hmm_loglik));

			kbits_t ref_kid = tmax_nbhd->states_sorted_sfx[i];
			kbits_t st = ref_kid >> shift;
			int bctx = (int) (ref_kid >> ctx_shift);

			tmax_base_emit_exp[bctx * 5 + xt] += norm_exp_count;
			tmax_qual_emit_exp[qual_contexts[tmax_kp] * (qmax << 1) + 
				(yt << 1) + (st == xt)] += norm_exp_count;
		}
		*/
	}
	else {
		atomic_add_u32(&d->n_reads_zero_likelihood, 1);
		//++d->n_reads_zero_likelihood;

		goto imm_return;
	}

	state_nbhd_t *next_nbhd = &T_nbhds[tmax];
	/* set BETA's at the tmax to be 1.0 */
	memset(next_nbhd->betas, 0, 
			next_nbhd->size * sizeof(*next_nbhd->betas));

	for (int t = tmax - 1; t >= tmin; t--) {
		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		obs_next_base = kmer_effective_base(rseq[t + kmer_len]);

		const int t_kp = t + kmer_len - 1;
		kbits_t xt = kmer_effective_base(rseq[t_kp]);
		int t_nbhd_size = curr_nbhd->size;

		int n_exp_trans = 0;
		for (int i = 0; i < t_nbhd_size; i++) {
			n_exp_trans += (curr_nbhd->kmer_ptrs[i]->n_trans + 1);
		}

		/* expected counts and expected transitions, the gamma's and
		 * xi's in Rabiner's notation. */
		double *t_exp_counts = mempool_alloc(mp, sizeof(*t_exp_counts) *
				t_nbhd_size);
		double *t_exp_trans = mempool_alloc(mp, sizeof(*t_exp_trans) *
				n_exp_trans);

		for (int i = 0; i < t_nbhd_size; i++) {
			t_exp_counts[i] = NAN;
		}
		for (int i = 0; i < n_exp_trans; i++) {
			t_exp_trans[i] = NAN;
		}

		double *kmer_exp_trans = t_exp_trans;

		int n_distinct_suffixes = curr_nbhd->n_uniq_sfx;
		int t_kmer_idx = 0;
		int tnext_kmer_idx = 0;
#if PMR_USE_LOG_ADD
		double sum_exp_count = 0.0;	/* for normalization */
#else
		int argmax_gamma = -1;
		double max_gamma = NAN;
#endif
		for (int i = 0; i < n_distinct_suffixes; i++) {
			kbits_t _tf = curr_nbhd->ba_pfx_sfx_flags[i];
			kbits_t ups_trans_flag = _tf & UINT32_MAX;
			kbits_t dns_trans_flag = _tf >> 32;

			/* the index of kmers that share this common suffix, not (k+1)-kmers */
			int kidx = 0;
			int n_self_trans = 0;
			for (int j = 0; j < 20; j++) {
				if ((ups_trans_flag & (1 << j)) == 0) continue;
				
				kmer_t *pk = curr_nbhd->kmer_ptrs[t_kmer_idx];

				/* reset value */
				for (int k = 0; k <= pk->n_trans; k++) {
					kmer_exp_trans[k] = NAN;
				}

				int _edist = _mm_popcnt_u64(curr_nbhd->edit_distances[t_kmer_idx]);
				double _alpha = curr_nbhd->alphas[t_kmer_idx];

				double k_tp[pk->n_trans+1];
				hmm_bwd_trans_prob(&k_tp[0], pk, 
						c_pi, c_pi_self, c_pd, c_pn, shift, kdict);

				/* only allow transitioning into observed base */
				kbits_t _tf = (_edist < dmax) ? 
					pk->trans_flag : ((1UL << obs_next_base) & pk->trans_flag);	

				if ((d->opt->unlimited_insertions & INS_MODE_EM) && 
						curr_nbhd->ins_prop[t_kmer_idx] > 0.0 && 
						((pk->id >> 62) == 0)) {
					_tf |= (1UL << 20);
				}

				/* self transition flag */
				register kbits_t _stf = (dns_trans_flag & _tf) >> 20;
				n_self_trans += _stf; 

				double _beta = hmm_compute_beta(
						dns_trans_flag, 
						tnext_kmer_idx,
						kidx,
						next_nbhd->betas,
						k_tp, 
						kmer_exp_trans,
						pk->trans_flag,
						_tf,
						next_nbhd->suffixes_order,
						next_nbhd->emit_dens,	
						_alpha,
						pk->id,
						next_nbhd->states_sorted_sfx
						);
		
				curr_nbhd->betas[t_kmer_idx] = _beta;
				kidx += _stf;

				double _gamma = PRODUCT(_alpha, _beta);
				t_exp_counts[t_kmer_idx] = _gamma;

#if PMR_USE_LOG_ADD
				sum_exp_count += EXP(_gamma);
#else
				update_max_quantity(_gamma, t_kmer_idx,
						max_gamma, argmax_gamma);
#endif
				t_kmer_idx++;
				kmer_exp_trans += (pk->n_trans + 1);
			}

			/* k-k, k-k+1, k+1-k, k+1-k+1 transitions + self transitions */
			tnext_kmer_idx += _mm_popcnt_u64((dns_trans_flag & ((1 << 20) - 1))) + 
				n_self_trans;
		}

#if PMR_USE_LOG_ADD
		double inv_sum_exp_count = -LOG(sum_exp_count);
#else
		double inv_sum_exp_count = -log_sum(curr_nbhd->size, t_exp_counts, 
				argmax_gamma, max_gamma);
#endif

		hmm_update_expected_vals(d, curr_nbhd, t, xt, 
				read->qscore[t_kp] - qoff, 
				t_exp_counts, t_exp_trans, inv_sum_exp_count, kmer_len);
	}

imm_return:
	mempool_destroy(mp);
}


void hmm_viterbi_wrapper(data *d, read_t *read, void *fdata)
{
	static char COMPLEMENTS[4] = {'T', 'G', 'A', 'C'};

	FILE *stream = d->errors_fp ? d->errors_fp : stdout;

	read_t *r1 = read_create();
	// allocate ample space just in case the read contains too many deletions
	r1->id = read->id;
	r1->sequence = calloc(read->length * 4, sizeof(char));
	r1->qscore = calloc(read->length * 4, sizeof(char));

	hmm_viterbi(d, read, (void *) r1);

	// discard the preconstructed first neighborhood since we are on the
	// reverse complement.
	d->preconstructed_nbhds[read->id].n = 0;
	d->preconstructed_nbhds[read->id].nbhd = NULL;

	// take the reverse complement of the read
	int rlen = r1->length;
	if (rlen == 0) {
		// if the read can't be decoded, take the original read
		r1->length = read->length;
		rlen = r1->length;
		memcpy(r1->sequence, read->sequence, rlen);
		memcpy(r1->qscore, read->qscore, rlen);
	}

	if (d->opt->viterbi_revcomp == 0) {
		// output the corrected read
		#pragma omp critical
		{
			fprintf(stream, "@%s\n%s\n+%s\n%s\n", 
					read->identifier, r1->sequence, 
					read->identifier, r1->qscore);
		}

		free(r1->sequence);
		free(r1->qscore);
		read_destroy(r1);

		return;
	}

	char *rseq_rc = calloc(rlen + 1, sizeof(char));
	char *qual_rc = calloc(rlen + 1, sizeof(char));

	for (int p = 0; p < rlen; ++p) {
		rseq_rc[p] = COMPLEMENTS[base_to_bits(r1->sequence[rlen - p - 1])];
		qual_rc[p] = r1->qscore[rlen - p - 1];
	}

	read_t *r2 = read_create();
	r2->id = r1->id;
	r2->sequence = r1->sequence;
	r2->qscore = r1->qscore;

	// reset memory to 0
	memset(r2->sequence, 0, read->length * 4);
	memset(r2->qscore, 0, read->length * 4);

	// use the reverse complement as the new "read"
	r1->sequence = rseq_rc;
	r1->qscore = qual_rc;

	hmm_viterbi(d, r1, (void *) r2);

	if (r2->length == 0) {
		r2->length = r1->length;
		memcpy(r2->sequence, r1->sequence, r1->length);
		memcpy(r2->qscore, r1->qscore, r1->length);
	}

	free(rseq_rc);
	free(qual_rc);
	read_destroy(r1);

	rlen = r2->length;
	rseq_rc = calloc(rlen + 1, sizeof(char));
	qual_rc = calloc(rlen + 1, sizeof(char));

	for (int p = 0; p < rlen; ++p) {
		rseq_rc[p] = COMPLEMENTS[base_to_bits(r2->sequence[rlen - p - 1])];
		qual_rc[p] = r2->qscore[rlen - p - 1];
	}

	// output the corrected read
	#pragma omp critical
	{
		fprintf(stream, "@%s\n%s\n+%s\n%s\n", 
				read->identifier, rseq_rc, read->identifier, qual_rc);
	}
	

	free(rseq_rc);
	free(qual_rc);

	free(r2->sequence);
	free(r2->qscore);

	read_destroy(r2);
}

/** 
 * hmm_viterbi
 *
 * Implementation of the Viterbi algorithm used for error decoding.
 *
 ****/
void hmm_viterbi(data *d, read_t *read, void *fdata)
{
	static char *BASES[] = {"A", "C", "T", "G",
		"AA", "CA", "TA", "GA",
		"AC", "CC", "TC", "GC",
		"AT", "CT", "TT", "GT",
		"AG", "CG", "TG", "GG",
		"self"};

	/*
	int actg_trans_count[4];
	int cum_actg_trans_count[4];
	*/
	kbits_t tnext_kid_list[4];
	kmer_t* tnext_kmer_list[4];

	const int kmer_len = d->opt->kmer_length;
	const int num_kmers = read->length - kmer_len + 1;
	const kbits_t klen_flag = (1UL << (kmer_len << 1)) - 1;
	const int shift = (kmer_len - 1) << 1;
	const int qmax = d->opt->max_qual_score + 1;
	const int qoff = d->opt->qual_score_offset;

	const double c_pi = d->ins_error_rate;
	const double c_pi_self = d->ins_error_rate_self;
	const double c_pd = d->del_error_rate;
	/* 1 - c_pi - c_pd */
	const double c_pn = d->error_free_p; 

	// corrected read
	read_t *rcor = (read_t *) fdata;

	if (read->length <= kmer_len) return;

	mempool_t *mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

	dict *kdict = d->kmer_dict;

	/* sequence */
	char *rseq = read->sequence;
	/*
	char *conv_q = convert_quality_scores(read->qscore, read->length, 
			33, mp);
	*/
	/* precompute quality score context information */
	/*
	int *qual_contexts = compute_qual_contexts(conv_q, read->length, 
			kmer_len, d->opt->qual_emit_context, qmax, mp);
	*/

	state_nbhd_t *T_nbhds = mempool_nalloc(mp,
			sizeof(*T_nbhds) * num_kmers, 16);

	/* flag for N (not determined) base, first kmer only */
	kbits_t obs_n_flag = kmer_n_base_flag(rseq, kmer_len);
	kbits_t next_obs_n_flag = kmer_n_base_flag(rseq + 1, kmer_len);

	double tnext_max_delta = NAN;
	int tnext_argmax_delta = -1;

	/* convert string literal to numeric representation */
	/*
	kbits_t observed_kmers[BITS_TO_U64(read->length)];
	read_seq_to_numeric(rseq, observed_kmers, read->length, kmer_len);
	*/

	/* set up first kmer neighborhood */
	

	int tmin = d->opt->adjust_read_start_pos & 0x2 ? d->read_start_pos[read->id] : 0;
	kbits_t obs_1st_kid = kmer_seq_to_id(rseq + tmin, kmer_len, 0);
	if (dict_find(kdict, kbits_cast_to_ptr(obs_1st_kid)) == NULL) {
		goto imm_return;
	}

	T_nbhds[tmin].size = hmm_load_preconstructed_nbhd(
			d, read->id, obs_1st_kid,
			obs_n_flag, &T_nbhds[tmin], kmer_len, true, 
			tmin > 0,
			kdict, mp);
	
	/*
	kbits_t obs_kid = observed_kmers[0] & ((1UL << (kmer_len << 1)) - 1);
	*/
	kbits_t obs_base = kmer_effective_base(rseq[tmin+kmer_len-1]);
	kbits_t obs_next_base = kmer_effective_base(rseq[tmin+kmer_len]);

	/* initial state distribution */
#ifdef PMR_RABINER_SCALING
	double rabiner_loglik = 0.0;
#endif

	//printf("[READ === %d]\n", read->id);
	int tmax = read->length - kmer_len;
	int t = tmin;
	for (; t < tmax; t++) {
		int tnext_kp = t + kmer_len;
		int tnext_nbhd_size;

		obs_next_base = kmer_effective_base(rseq[t + kmer_len]);
		/*
		obs_kid = bitarray_get(observed_kmers, t << 1, kmer_len << 1);
		*/

		register kbits_t _new_nf = ((obs_next_base >> 2) << 1) | 
			(obs_next_base >> 2);
		next_obs_n_flag = (next_obs_n_flag >> 2) | (_new_nf << shift);

		state_nbhd_t *curr_nbhd = &T_nbhds[t],
					 *next_nbhd = &T_nbhds[t+1];

		const kbits_t dmax = d->opt->max_edit_dist;

		int t_nbhd_size = curr_nbhd->size;
		/*
		if (t_nbhd_size > max_nbhd_size) {
			printf("Max nbhd size: %d\n", t_nbhd_size);
			max_nbhd_size = t_nbhd_size;
		}
		*/

		kbits_t *t_states_sorted_sfx = T_nbhds[t].states_sorted_sfx;
		double *t_states_deltas = T_nbhds[t].alphas;

		/* ----- II. compute one iteration of forward algorithm ----- */

		/* --- II a. set up kmer neighborhood at next position --- */
		//memset(actg_trans_count, 0, sizeof(int) << 2);

		/* one linear scan to determine # of distinct suffix, 
		 * compute downstream actg nucleotide counts etc. */
		int n_distinct_suffixes = hmm_count_distinct_suffixes(curr_nbhd, 
				&tnext_nbhd_size, dmax, obs_n_flag,
				obs_next_base, d->transition_p, kmer_len,
				d->opt->unlimited_insertions & INS_MODE_VITERBI);

		curr_nbhd->n_uniq_sfx = n_distinct_suffixes;

		tnext_max_delta = NAN;
		tnext_argmax_delta = -1;

		/* number of kmers in next position */
		if (tnext_nbhd_size == 0) break;

		/* allocate memory given tnext_nbhd_size */
		char *pmem = mempool_nalloc(mp,
				nbhd_alloc_size(tnext_nbhd_size), 16);
		candidate_t *tnext_candidates = mempool_alloc(mp,
			(tnext_nbhd_size + 1) * sizeof(*tnext_candidates));
		double *tnext_alphas = mempool_alloc(mp,
			tnext_nbhd_size * sizeof(double));
		double *tnext_betas = mempool_alloc(mp,
			tnext_nbhd_size * sizeof(double));

		next_nbhd->size = tnext_nbhd_size;
		hmm_setup_nbhd_ptrs(next_nbhd, tnext_nbhd_size, pmem);


		//printf("[t=%2d / size:%d]\n", t+2, tnext_nbhd_size);

		/* FIXME: is there any means to wrap following code block into a 
		 * macro template? */

		int kmer_idx = 0, index_pfx = 0;

		register kbits_t mismatch_flag = ~(1UL << obs_next_base) |
			~((obs_next_base >> 2) - 1UL);

		for (int i = 0; i < n_distinct_suffixes; i++) {
			kbits_t _repr_kid = t_states_sorted_sfx[kmer_idx];
			//kbits_t mut_flag = ((_repr_kid ^ obs_kid) | obs_n_flag) & 3UL;
			/*
			kbits_t sfx_hd = bitarray_get_pwr2(curr_nbhd->ba_hamming_dist, 
					i, 2);
			*/
			kbits_t common_sfx = kbits_trans_suffix(_repr_kid);

			/* v 4 bits                   v 4 bits (higher)
			 *
			 * A  -                       - A
			 * C  |  -- k-1 substring --  | C
			 * T  |                       | T
			 * G  -                       - G
			 */
			kbits_t _tf = curr_nbhd->ba_pfx_sfx_flags[i];
			kbits_t ups_trans_flag = _tf & UINT32_MAX;
			/* which kmers/nucleotides current suffix transition into */
			kbits_t dns_trans_flag = _tf >> 32;

			int n_dns_trans = _mm_popcnt_u64(dns_trans_flag);
			int n_ups_trans = _mm_popcnt_u64(ups_trans_flag);

			/*
			printf("[T] ");
			kmer_id_to_seq(_repr_kid, kmer_len);
			printf(" --> ");
			*/

			uint64_t k_sfx_edist = UINT64_MAX;
			int k_sfx_n_err = kmer_len + 1;
			for (int j = 0; j < n_ups_trans; j++) {
				kbits_t _id = t_states_sorted_sfx[kmer_idx + j];
				uint64_t _edist = kbits_suffix_edist(_id,
						curr_nbhd->edit_distances[kmer_idx + j]);
				int _ne = _mm_popcnt_u64(_edist);
				if (_ne < k_sfx_n_err) {
					k_sfx_n_err = _ne;
					k_sfx_edist = _edist;
				}
			}

			for (int j = 0; j < 21; j++) {
				if (((dns_trans_flag >> j) & 1) == 0) continue;

				//printf(" %s ", BASES[j]);
				kbits_t next_bases = (j < 4) ? j : j - 4;
				int len = 0;

				if (j < 20) {
					tnext_kid_list[0] = (common_sfx & KBITS_KLEN_REMOVAL_FLAG) | 
						(next_bases << shift) | ((uint64_t) (j > 3) << 62);
					tnext_kmer_list[0] = dict_find(kdict, 
						kbits_cast_to_ptr(tnext_kid_list[0]))->value;
					
					len = 1;
				}
				else {
					len = 0;
					kbits_t _ktp = trans_flag2packed(ups_trans_flag & 15);
					int _nk = _ktp & 7;
					_ktp >>= 3;
					for (int k = 0; k < _nk; k++, _ktp >>= 2) {
						kbits_t ups_base = _ktp & 3;

						int _idx = kmer_idx + k;
						uint64_t _edist = _mm_popcnt_u64(
								curr_nbhd->edit_distances[_idx]);
						if (_edist < dmax || 
								((d->opt->unlimited_insertions & INS_MODE_VITERBI) &&
								 curr_nbhd->ins_prop[_idx] > 0.0)) {
							tnext_kid_list[len] = t_states_sorted_sfx[_idx];
							tnext_kmer_list[len] = curr_nbhd->kmer_ptrs[_idx];
							++len;
						}
					}
				}

				kbits_t st = (j == 20) ? kbits_last_base(_repr_kid, kmer_len) : 
					((j < 4) ? j : ((j-4) >> 2));

				int yctx = (j == 20) ? QCTX_SELF_TRANS :
					((j < 4) ? QCTX_KMER_EQUAL : QCTX_KP1MER);

				double emit_dens = hmm_compute_emit_dens(
						d, st, obs_next_base, 
						read->qscore[tnext_kp] - qoff,
						yctx);
				
				for (int k = 0; k < len; k++) {
					/* for every kmer only at current position, add itself
					to the candidate list at the next position */
					kbits_t tnext_kid = tnext_kid_list[k];
					kbits_t tnext_kmer = tnext_kmer_list[k];

					/*
					printf("NEXT(%d): ", j);
					kmer_id_to_seq(tnext_kid, kmer_len);
					putchar('\n');
					*/

					double tp2j[n_ups_trans];
					hmm_fwd_trans_prob(&tp2j[0], ups_trans_flag, 
							kmer_idx, curr_nbhd->kmer_ptrs, 
							curr_nbhd->edit_distances,
							curr_nbhd->ins_prop,
							j, tnext_kid, c_pi, c_pd, c_pn, 
							shift, dmax, kdict, 
							d->opt->unlimited_insertions & INS_MODE_VITERBI);

					/* we will be use betas as the psi's following Rabiner's
					 * notation. */
					double delta;
					int argmax_delta = hmm_compute_delta(
							n_ups_trans,
							ups_trans_flag,
							kmer_idx,
							t_states_deltas,
							tp2j,
							emit_dens, &delta);

					//printf(">> alpha: %lg\n", EXP(alpha));

					/* temporarily sorted by prefix */
					tnext_candidates[index_pfx].id = tnext_kid;
					tnext_candidates[index_pfx].kmer_ptr = tnext_kmer;
					tnext_candidates[index_pfx].pfx_order = index_pfx;
					tnext_candidates[index_pfx].edist = k_sfx_edist | 
						((j != obs_next_base) << (kmer_len-1));
					tnext_candidates[index_pfx].is_insertion = (j == 20);
					tnext_candidates[index_pfx].emit_dens = emit_dens;

					tnext_alphas[index_pfx] = delta;
					tnext_betas[index_pfx++] = (double) argmax_delta;
				}
			}

			//printf("\n");

			kmer_idx += n_ups_trans;
		}

		qsort(tnext_candidates, tnext_nbhd_size, sizeof(candidate_t), 
			kmer_cmp_suffix);

		/* terminator */
		tnext_candidates[tnext_nbhd_size].id = UINT64_MAX;

		/*
		if (read->id == 0) {
			printf("[t=%d] %.*s -> %c\nNext candidates [%d]:\n", 
					t, kmer_len, rseq+t, rseq[t+kmer_len], tnext_nbhd_size);
		}
		*/

		/* reduce the potentially redundant neighborhood */
		int reduced_nbhd_size = 0;
		int prev_i = 0;
		kbits_t prev_cand_id = (tnext_nbhd_size == 0) ? UINT64_MAX : tnext_candidates[0].id;
		candidate_t *prev_cand = NULL;
		double _delta, _argmax_delta;
		double _ins_prop;
		uint64_t _argmax_edist, _min_edist;

		for (int i = 0; i <= tnext_nbhd_size; i++) {
			candidate_t *cand = &tnext_candidates[i];
			kbits_t cand_id = cand->id;

			/*
			if (i < tnext_nbhd_size) {
				printf("CANDIDATE %d ", i);
				kmer_id_to_seq(cand_id, kmer_len);
				printf(" %lg\n", EXP(tnext_alphas[cand->pfx_order]));
			}
			*/

			if (cand_id != prev_cand_id) {
				//kmer_id_to_seq(prev_cand->id, kmer_len);
				if (i > (prev_i + 1)) {
					/* to reduce redundant hidden states */
					_delta = NAN;
					_min_edist = UINT64_MAX;

					_argmax_delta = tnext_betas[tnext_candidates[prev_i].pfx_order];
					for (int j = prev_i; j < i; ++j) {
						candidate_t *cj = &tnext_candidates[j];	
						next_nbhd->suffixes_order[cj->pfx_order] = 
							reduced_nbhd_size;

						/* reduce by finding the maximizer */
						if (isnan(_delta) || _delta < tnext_alphas[cj->pfx_order]) {
							_delta = tnext_alphas[cj->pfx_order];
							_argmax_delta = tnext_betas[cj->pfx_order];
							_argmax_edist = cj->edist;
							_ins_prop = (double) cj->is_insertion;
						}
						if (_min_edist > cj->edist) 
							_min_edist = cj->edist;
					}
					//printf(" -[R]-> delta: %lg\n", EXP(_delta));
				}
				else {
					/* no redundant hidden states */
					_ins_prop = (double) prev_cand->is_insertion;
					_argmax_edist = prev_cand->edist;
					_min_edist = prev_cand->edist;
					_delta = tnext_alphas[prev_cand->pfx_order];
					_argmax_delta = tnext_betas[prev_cand->pfx_order];
					next_nbhd->suffixes_order[prev_cand->pfx_order] = 
						reduced_nbhd_size;
					//printf(" -[-]-> delta: %lg\n", EXP(_delta));
				}

				next_nbhd->alphas[reduced_nbhd_size] = _delta;
				next_nbhd->betas[reduced_nbhd_size] = _argmax_delta;
				next_nbhd->states_sorted_sfx[reduced_nbhd_size] = prev_cand->id;
				next_nbhd->kmer_ptrs[reduced_nbhd_size] = prev_cand->kmer_ptr;
				next_nbhd->emit_dens[reduced_nbhd_size] = prev_cand->emit_dens;
				//next_nbhd->edit_distances[reduced_nbhd_size] = _argmax_edist; 
				// strategy 2: use minimal edit distance
				next_nbhd->edit_distances[reduced_nbhd_size] = _min_edist; 
				next_nbhd->ins_prop[reduced_nbhd_size] = _ins_prop; 

				bitarray_set_pwr2(next_nbhd->ba_kmer_trans_flag,
					reduced_nbhd_size, prev_cand->kmer_ptr->trans_flag, 5);

				update_max_quantity(_delta, reduced_nbhd_size, 
						tnext_max_delta, tnext_argmax_delta);

				/*
				if (read->id == 0) {
					printf("%d ", reduced_nbhd_size);
					kmer_id_to_seq(prev_cand->id, kmer_len);
					printf(" %lg [>%d//%d]\n", EXP(_delta), (int) _argmax_delta, _argmax_edist);
				}
				*/

				++reduced_nbhd_size;
				prev_i = i;
			}

			prev_cand = cand;
			prev_cand_id = cand_id;
		}

		next_nbhd->size = reduced_nbhd_size;

		obs_base = obs_next_base;
		obs_n_flag = next_obs_n_flag;

		//if (read->id == 0) getchar();
	}

	kbits_t *optimal_pathway = mempool_alloc(mp, 
			(tmax + 1) * sizeof(*optimal_pathway));

	if (t == tmax && tnext_argmax_delta >= 0) {
		state_nbhd_t *nbhd = &T_nbhds[tmax];
		optimal_pathway[tmax] = nbhd->states_sorted_sfx[tnext_argmax_delta];

		int optimal_st_idx = (int) nbhd->betas[tnext_argmax_delta];
		for (int t = tmax - 1; t >= tmin; t--) {
			nbhd = &T_nbhds[t];

			optimal_pathway[t] = nbhd->states_sorted_sfx[optimal_st_idx];
			optimal_st_idx = (int) nbhd->betas[optimal_st_idx];
		}

		// if the read is started at a shifted position, we need to copy all
		// the bases from the original sequence to the "corrected" sequence
		// before the starting position.
		for (int t = 0; t < tmin; ++t) {
			rcor->sequence[t] = rseq[t];
		}

		/* read identifier */
		kmer_sprint_id(optimal_pathway[tmin], kmer_len, rcor->sequence + tmin);

		// corrected read length
		int cor_read_len = kmer_len + tmin;
		for (int t = tmin + 1; t <= tmax; t++) {
			kbits_t opt_id = optimal_pathway[t];
			/* skip if transition into itself */
			if ((opt_id) == (optimal_pathway[t-1])) continue;

			int l = opt_id >> 62;
			if (l == 0) {
				rcor->sequence[cor_read_len++] = BASES[opt_id >> shift][0];
			}
			else if (l == 1) {
				// (k+1)-mer
				memcpy(&rcor->sequence[cor_read_len], 
						BASES[4+((opt_id & KBITS_KLEN_REMOVAL_FLAG) >> shift)], 2);
				cor_read_len += 2;
			}
		}

		rcor->length = cor_read_len;

		/* quality score */
		if (cor_read_len <= read->length)
			memcpy(rcor->qscore, read->qscore, cor_read_len);
		else {
			// need to pad some "pseudo quality score".
			int pad_len = cor_read_len - read->length;
			memcpy(rcor->qscore, read->qscore, cor_read_len);
			memcpy(&rcor->qscore[read->length], read->qscore, pad_len);
		}
	}

#if 0
	FILE *stream = d->errors_fp ? d->errors_fp : stdout;
	if (t == tmax && tnext_argmax_delta >= 0) {
		state_nbhd_t *nbhd = &T_nbhds[tmax];
		optimal_pathway[tmax] = nbhd->states_sorted_sfx[tnext_argmax_delta];

		int optimal_st_idx = (int) nbhd->betas[tnext_argmax_delta];
		for (int t = tmax - 1; t >= 0; t--) {
			nbhd = &T_nbhds[t];

			optimal_pathway[t] = nbhd->states_sorted_sfx[optimal_st_idx];
			optimal_st_idx = (int) nbhd->betas[optimal_st_idx];
		}

		/* read identifier */
		fprintf(stream, "%s\n", read->identifier);
		kmer_print_id(optimal_pathway[0], kmer_len, stream);

		/*
		printf("[READ %d]\n  0 ", read->id+1);
		kmer_id_to_seq(optimal_pathway[0], kmer_len);
		*/

		// corrected read length
		int cor_read_len = kmer_len;
		for (int t = 1; t <= tmax; t++) {
			/*
			printf("\n%3d ", t);
			kmer_id_to_seq(optimal_pathway[t], kmer_len);
			*/

			kbits_t opt_id = optimal_pathway[t];
			/* skip if transition into itself */
			if ((opt_id) == (optimal_pathway[t-1])) continue;

			int l = opt_id >> 62;
			if (l == 0) {
				fprintf(stream, "%s", BASES[opt_id >> shift]);
				++cor_read_len;
			}
			else if (l == 1) {
				fprintf(stream, "%s", BASES[
						4+((opt_id & KBITS_KLEN_REMOVAL_FLAG) >> shift)]);

				cor_read_len += 2;
			}
		}

		/* quality score */
		fprintf(stream, "\n+%s\n", read->identifier + 1);
		if (cor_read_len <= read->length)
			fprintf(stream, "%.*s\n", cor_read_len, read->qscore);
		else {
			// need to pad some "pseudo quality score".
			fprintf(stream, "%.*s", read->length, read->qscore);
			int pad_len = cor_read_len - read->length;
			char qpad[pad_len + 1];
			qpad[pad_len] = '\0';
			memcpy(qpad, read->qscore, pad_len);
			fprintf(stream, "%s\n", qpad);
		}
	}
	else {
		++d->n_reads_cannot_decode;
		// failed to perform error correction: output the original reads.
		fprintf(stream, "%s\n%.*s\n+%s\n%.*s\n", read->identifier,
				read->length, read->sequence, read->identifier + 1,
				read->length, read->qscore);
	}
#endif

	//getchar();

imm_return:
	mempool_destroy(mp);
}

static void hmm_init_model_params(data *d)
{
	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) << 1;
	int penal_init_params = d->opt->penalize_init_params;

	/* --- emission --- */
	INFO("Initiating emission probabilities...\n");
	const double init_qual_emit_p = LOG(1.0 / 
			(d->opt->max_qual_score + 1));

	for (int i = 0; i < 4; ++i) {
		for (int obs_b = 0; obs_b < 5; ++obs_b) {
			d->base_emission_p[i * 5 + obs_b] = (i == obs_b) ?
				LOG(1 - d->emit_error_rate) : LOG(d->emit_error_rate / 3);
		}
	}

	for (int i = 0; i < d->qual_ctx_radix; ++i) {
		d->qual_emission_p[i] = init_qual_emit_p;
	}
}

static void hmm_init_e_step(data *d)
{
	/* reset the emission densities table */
	for (int i = 0; i < d->base_ctx_radix; ++i) {
		d->base_emit_exp[i] = 0.0;
	}

	for (int i = 0; i < d->qual_ctx_radix; ++i) {
		d->qual_emit_exp[i] = 0.0;
	}

	d->all_sum_exp_trans = 0.0;
	d->ins_sum_exp_trans = 0.0;
	d->ins_sum_exp_trans_self = 0.0;
	d->del_sum_exp_trans = 0.0;

	/* reset table for recording transition probabilities */
	int i;
	dictentry *de;
	for (de = d->kmer_dict->table, i = d->kmer_dict->used; i > 0; de++) {
		if (de->value != NULL) {
			i--;
			kmer_t *pk = (kmer_t *) de->value;

			if (pk->id >> 62) continue;

			pk->exp_count = 0.0;
			for (int j = 0; j < pk->n_ktrans; j++) {
				pk->exp_trans_p[j] = 0.0;
			}
		}
	}
}

static double hmm_compute_beta(kbits_t dns_trans_flag,
		int dns_kmer_idx, int pfx_idx, double *tnext_state_betas, 
		double *kmer_trans_p, double *kmer_exp_trans,
		kbits_t k_unc_tf, kbits_t k_dmax_tf,
		int *tnext_sfx_order, double *emit_dens, double alpha,
		kbits_t kid, kbits_t *tnext_states)
{
	dbl_t *tails = (dbl_t *) kmer_trans_p;
	dbl_t max_tail = {.dbl = NAN};
	int argmax_tail = 0;
	
	int dnsidx = 0;
	int transidx = 0;
	for (int j = 0; j < 20; j++) {
		/* non-admissible transition */
		if ((dns_trans_flag & (1 << j)) == 0) continue;
		
		/* check if kmer (under dmax constraint) can transition into the jth suffix */
		if ((k_dmax_tf & (1 << j))) {
			int jth_tran = _mm_popcnt_u64(k_unc_tf & ((1 << j) - 1));

			int beta_idx = tnext_sfx_order[dns_kmer_idx + dnsidx];
			dbl_t _tail = {.dbl = PRODUCT(tails[jth_tran].dbl,
					PRODUCT(emit_dens[beta_idx],
					tnext_state_betas[beta_idx]))};
			tails[transidx] = _tail;

			uint64_t ge_mask = (uint64_t) (isnan(max_tail.dbl) || 
					(max_tail.dbl < _tail.dbl)) - 1UL;
			uint64_t inv_ge_mask = ~ge_mask;

			max_tail.u64 = (ge_mask & max_tail.u64) | 
				(inv_ge_mask & _tail.u64);
			argmax_tail = (argmax_tail & ge_mask) | 
				(transidx & inv_ge_mask);

			kmer_exp_trans[jth_tran] = PRODUCT(alpha, _tail.dbl);
			transidx++;
		}

		dnsidx++;
	}

	/* check self-transition */
	if (((dns_trans_flag & k_dmax_tf) >> 20) && (pfx_idx < 4)) {
		int beta_idx = tnext_sfx_order[dns_kmer_idx + dnsidx + pfx_idx];
		int jth_tran = _mm_popcnt_u64(k_unc_tf & ((1 << 20) - 1));
		tails[transidx].dbl = PRODUCT(kmer_trans_p[jth_tran],
				PRODUCT(emit_dens[beta_idx],
					tnext_state_betas[beta_idx]));
				/*
				PRODUCT(emit_dens[beta_idx], 
					tnext_state_betas[beta_idx]));
				*/

		uint64_t ge_mask = (uint64_t) (isnan(max_tail.dbl) || 
				(max_tail.dbl < tails[transidx].dbl)) - 1UL;
		uint64_t inv_ge_mask = ~ge_mask;

		max_tail.u64 = (ge_mask & max_tail.u64) | 
			(inv_ge_mask & tails[transidx].u64);
		argmax_tail = (argmax_tail & ge_mask) | 
			(transidx & inv_ge_mask);

		kmer_exp_trans[jth_tran] = PRODUCT(alpha, tails[transidx].dbl);
		transidx++;
	}

	double _expsum = 0.0; 
	for (int j = 0; j < transidx; j++) {
		register uint64_t argmax_mask = ((j == argmax_tail) - 
				1UL);
		register dbl_t summand = {.dbl = EXP(tails[j].dbl - 
				max_tail.dbl)};
		summand.u64 &= argmax_mask;

		_expsum += summand.dbl; 
	}

	return PRODUCT(max_tail.dbl, log1p(_expsum));
}

/* hmm_compute_emit_dens
 * st: s_{[t]}, the k-th (last) base of the hidden state s_t.
 * xt: x_{[t]}, the k-th (last) base of the observed state x_t.
 * yctx: context/scenario for emitting the quality score.
 */
inline static double hmm_compute_emit_dens(
		data *d, kbits_t st, kbits_t xt, int yt, int yctx)
{
	if (yctx == QCTX_KMER_EQUAL) 
		yctx += (xt != st);

	return PRODUCT(d->qual_emission_p[
			(d->opt->max_qual_score + 1) * yctx + yt],
			d->base_emission_p[st * 5 + xt]);
}

#if 0
inline static double hmm_normalize_emission(double *emit_p, kbits_t enf, 
		kbits_t xt)
{
	if (enf == 15) return LOG(emit_p[xt]);
	else {
		double sum_emit_p = emit_p[BASE_N];
		register kbits_t en_packed = trans_flag2packed(enf);
		register int i = (int) en_packed & 7;
		for (en_packed >>= 3; i > 0; i--, en_packed >>= 2) {
			kbits_t base = en_packed & 3;
			sum_emit_p += emit_p[base];
		}

		return LOG(emit_p[xt] / sum_emit_p);
	}
}
#endif

/**
 * hmm_fwd_trans_prob
 *
 * Compute the transition probabilities in forward algorithm 
 *
 * */
static void hmm_fwd_trans_prob(double *tp, kbits_t ups_trans_flag, int ups_kmer_idx, 
		kmer_t **t_kmers, uint64_t *edists, double *ins_prop, kbits_t dns_bases, 
		kbits_t tnext_kid, double c_pi, double c_pd, double c_pn, int shift, int dmax, 
		dict *kdict, int unlimited_ins)
{
	int nt = _mm_popcnt_u64(ups_trans_flag);
	int n = 0;
	for (int i = 0; i < 20; i++) {
		if (((ups_trans_flag >> i) & 1) == 0) continue;

		register int _idx = ups_kmer_idx + n;
		kbits_t obs_next_base = tnext_kid >> shift;
		kmer_t *ref_k = t_kmers[_idx];
		int _edist = _mm_popcnt_u64(edists[_idx]);

		kbits_t _tf = _edist < dmax ? ref_k->trans_flag :
			ref_k->trans_flag & (1 << obs_next_base);
		if (unlimited_ins && ins_prop[_idx] > 0.0 && 
				(ref_k->id >> 62) == 0) _tf |= (1UL << 20);
			
		if ((_tf & (1 << dns_bases)) != 0) {
			if (dns_bases < 4) {
				int dns_trans_idx = _mm_popcnt_u64(ref_k->trans_flag 
						& ((1UL << dns_bases) - 1));
				tp[n] = PRODUCT(c_pn, 
						ref_k->transition_p[dns_trans_idx]);
			}
			else if (dns_bases == 20) {
				tp[n] = ref_k->id == tnext_kid ? c_pi : NAN;
			}
			else {
				/* -> (k+1)-mer */
				register kbits_t base1 = (dns_bases-4) & 3;
				register kbits_t base2 = (dns_bases-4) >> 2;
				kbits_t kid2 = kbits_raw_suffix(ref_k->id) | (base1 << shift);
				kmer_t *ref_k2 = dict_find(kdict, 
						kbits_cast_to_ptr(kid2))->value;

				int dns_tidx1 = _mm_popcnt_u64(ref_k->trans_flag & 
						((1UL << base1) - 1));
				int dns_tidx2 = _mm_popcnt_u64(ref_k2->trans_flag & 
						((1UL << base2) - 1));

				tp[n] = PRODUCT(c_pd, PRODUCT(
							ref_k->transition_p[dns_tidx1], 
							ref_k2->transition_p[dns_tidx2]));
			}
		}
		else {
			tp[n] = NAN;
		}

		n++;
	}
}

static void hmm_bwd_trans_prob(double *tp, kmer_t *pk, 
		double c_pi, double c_pi_self, double c_pd, double c_pn, int shift, dict *kdict)
{
	int n = 0;
	
	for (int j = 0; j <= 20; j++) {
		if ((pk->trans_flag & (1 << j)) == 0) continue;
				
		if (j < 4) {
			int dns_trans_idx = _mm_popcnt_u64(pk->trans_flag 
					& ((1UL << j) - 1));
			tp[n] = PRODUCT(c_pn, pk->transition_p[dns_trans_idx]);
		}
		else if (j == 20) {
			tp[n] = c_pi;
		}
		else {
			register kbits_t base1 = (j-4) & 3;
			register kbits_t base2 = (j-4) >> 2;
			kbits_t kid2 = kbits_raw_suffix(pk->id) | (base1 << shift);
			kmer_t *pk2 = dict_find(kdict, kbits_cast_to_ptr(kid2))->value;

			int dns_tidx1 = _mm_popcnt_u64(pk->trans_flag & 
					((1UL << base1) - 1));
			int dns_tidx2 = _mm_popcnt_u64(pk2->trans_flag & 
					((1UL << base2) - 1));

			tp[n] = PRODUCT(c_pd, PRODUCT(pk->transition_p[dns_tidx1], 
						pk2->transition_p[dns_tidx2]));
		}

		++n;
	}
}

static double hmm_compute_alpha(int n_ups_trans,
	kbits_t ups_trans_flag, int ups_kmer_idx, double *t_state_alphas, 
	double *summands, double emit_dens)
{
	dbl_t max_summand = {.dbl = NAN};
	int argmax_idx = 0;

	for (int i = 0; i < n_ups_trans; i++) {
		register int _idx = ups_kmer_idx + i;
		register dbl_t _summand = {.dbl = PRODUCT(t_state_alphas[_idx], 
					summands[i])};

		summands[i] = _summand.dbl;

		uint64_t ge_mask = (uint64_t) (
				isnan(max_summand.dbl) || 
				(max_summand.dbl < _summand.dbl)) - 1UL;
		uint64_t inv_ge_mask = ~ge_mask;

		max_summand.u64 = (ge_mask & max_summand.u64) | 
			(inv_ge_mask & _summand.u64);
		argmax_idx = (ge_mask & argmax_idx) | (inv_ge_mask & i);

	}

/* compute the sum of all summands gives us alpha */
	double _expsum = 0.0;
	for (int i = 0; i < n_ups_trans; i++) {
		register uint64_t argmax_mask = ((i == argmax_idx) - 1UL);
		register dbl_t _summand = {
			.dbl = EXP(summands[i] - max_summand.dbl)};
		_summand.u64 &= argmax_mask;

		_expsum += _summand.dbl; 
	}

	return PRODUCT(emit_dens,
			PRODUCT(max_summand.dbl, log1p(_expsum)));
}

inline static int hmm_compute_delta(int n_ups_trans,
		kbits_t ups_trans_flag, int ups_kmer_idx, double *t_state_deltas,
		double *tp, double emit_dens, double *delta)
{
	dbl_t max_tail = {.dbl = NAN};
	int argmax_tail = 0;

	int n = 0;
	for (int i = 0; i < n_ups_trans; i++) {
		register int _idx = ups_kmer_idx + i;
		register dbl_t _tail = {.dbl = PRODUCT(t_state_deltas[_idx], 
					tp[i])};

		/* (greater or equal) ge_mask is set to 0, 
		 * if max_summand < _summand */
		uint64_t ge_mask = (uint64_t) (isnan(max_tail.dbl) || 
				(max_tail.dbl < _tail.dbl)) - 1UL;
		uint64_t inv_ge_mask = ~ge_mask;

		max_tail.u64 = (ge_mask & max_tail.u64) | 
			(inv_ge_mask & _tail.u64);
		argmax_tail = (ge_mask & argmax_tail) | (inv_ge_mask & _idx);
	}

	*delta = PRODUCT(max_tail.dbl, emit_dens);
	return argmax_tail;	
}

static inline void hmm_correct_errors(data *d, read_t *read, int tmax,
		kbits_t st, kbits_t xt, kbits_t xt_n_flag, int kmer_len)
{
	static char BASES[5] = {'A', 'C', 'T', 'G', 'N'};

	for (int i = 1; i <= kmer_len; i++) {
		if (xt_n_flag & 3) {
			printf("%s %d %c %c\n", read->identifier + 1,
					tmax + i, BASES[st & 3], 'N');
		}
		else if ((xt & 3) != (st & 3)) {
			printf("%s %d %c %c\n", read->identifier + 1,
					tmax + i, BASES[st & 3], BASES[xt & 3]);
		}

		st >>= 2;
		xt >>= 2;
		xt_n_flag >>= 2;
	}
}

static inline void hmm_correct_one_error(data *d, read_t *read, int t,
		kbits_t st, kbits_t xt)
{
	static char BASES[5] = {'A', 'C', 'T', 'G', 'N'};

	if (st != xt) {
		printf("%s %d %c %c\n", read->identifier + 1, 
				t + 1, BASES[st], BASES[xt]);
	}
}

static int hmm_count_distinct_suffixes(state_nbhd_t *nbhd, 
		int *tnext_nbhd_size, kbits_t t_dmax, kbits_t obs_n_flag, 
		kbits_t obs_next_base, double *unique_trans_p, int kmer_len,
		int unlimit_ins)
{
	register int t_nbhd_size = nbhd->size;
	kbits_t *t_sts_sfx = nbhd->states_sorted_sfx;

	kbits_t prev_k_sfx = kbits_trans_suffix(t_sts_sfx[0]);

	kbits_t pfx_trans_flag = 0, sfx_trans_flag = 0;
	kbits_t _tfhd = 0;
	int n_distinct_suffixes = 0;
	int sfx_idx_lb = 0;

	int tnext_size = 0;

	/* FIXME: consider use larger step size for i to optimize cache access */
	for (int i = 0; i < t_nbhd_size; i++) {
		/* kmers, sorted by suffix */
		register kbits_t kid = t_sts_sfx[i];
		kbits_t kpfxlen = kid >> 62;
		kbits_t k_sfx = kbits_trans_suffix(kid);
		/*
		*/

		/* change in kmer suffix */
		if ((k_sfx & KBITS_KLEN_REMOVAL_FLAG) != (prev_k_sfx & KBITS_KLEN_REMOVAL_FLAG)) {
			/*
			bitarray_set_pwr2(nbhd->ba_distinct_sfx, 
					n_distinct_suffixes, 
					(i-sfx_idx_lb-1), 1);
			*/
			/* overwrite the original Hamming distance */
			/*
			bitarray_set_pwr2(nbhd->ba_hamming_dist, 
					n_distinct_suffixes,
					_tfhd >> 4, 2);
			*/
			nbhd->ba_pfx_sfx_flags[n_distinct_suffixes] = 
					(sfx_trans_flag << 32) | pfx_trans_flag;

			/*
			for (int j = 0; j < 4; j++) {
				actg_trans_count[j] += ((sfx_trans_flag >> j) & 1);
			}
			*/

			tnext_size += _mm_popcnt_u64(sfx_trans_flag & ((1UL << 20) - 1));

			++n_distinct_suffixes;

			sfx_trans_flag = 0;
			pfx_trans_flag = 0;
			sfx_idx_lb = i;
		}

		/* unconstrained transitions */
		kbits_t k_unc_tf = bitarray_get_pwr2(nbhd->ba_kmer_trans_flag,
				i, 5);
		int edist = _mm_popcnt_u64(nbhd->edit_distances[i]);
		/*
		kbits_t k_hamming_dist = bitarray_get_pwr2(nbhd->ba_hamming_dist, 
				i, 2);
		*/

		/* _tfhd: 8 bits packed integer in following format:
		 * HDHD TFTF
		 */
		
		kbits_t k_tf = trans_admissible_transitions(kid, obs_n_flag, k_unc_tf, 
			obs_next_base, kmer_len, (edist < t_dmax));

		// always allow self-transition if it cames from a self-transition
		if (unlimit_ins && nbhd->ins_prop[i] > 0.0 && kpfxlen == 0) {
			if (k_tf & (1UL << 20) == 0) {
				nbhd->ins_prop[i] = 1.0;
				nbhd->alphas[i] = PRODUCT(nbhd->alphas[i], LOG(nbhd->ins_prop[i]));
			}
			k_tf |= (1UL << 20);
		}

		if (kpfxlen == 0 && (k_tf & (1UL << 20))) {
			tnext_size++;
		}

		kbits_t discard_pfx = kbits_discarded_prefix(kid) + (kpfxlen << 2);
		pfx_trans_flag |= (1UL << discard_pfx);
		sfx_trans_flag |= k_tf;

		prev_k_sfx = k_sfx;
	}

	/*
	bitarray_set_pwr2(nbhd->ba_distinct_sfx, 
			n_distinct_suffixes, 
			(t_nbhd_size-sfx_idx_lb-1), 1);
	*/

	/* overwrite the original Hamming distance */
	/*
	bitarray_set_pwr2(nbhd->ba_hamming_dist, 
			n_distinct_suffixes,
			_tfhd >> 4, 2);
	*/

	nbhd->ba_pfx_sfx_flags[n_distinct_suffixes] = (sfx_trans_flag << 32) | pfx_trans_flag;

	tnext_size += _mm_popcnt_u64(sfx_trans_flag & ((1UL << 20) - 1));

	/*
	for (int j = 0; j < 4; j++) {
		actg_trans_count[j] += ((sfx_trans_flag >> j) & 1);
	}
	*/

	++n_distinct_suffixes;

	*tnext_nbhd_size = tnext_size;

	return n_distinct_suffixes;
}

static inline int hmm_load_preconstructed_nbhd(data *d, int rid, 
		kbits_t obs_1st_kid, kbits_t obs_n_flag, state_nbhd_t *first_nbhd, 
		int kmer_len, int compute_emit_dens, int observed_only,
		dict *kdict, mempool_t *mp)
{
	kmer_t *ref_kmer = dict_find(kdict, 
			kbits_cast_to_ptr(obs_1st_kid))->value;

	kmer_nbhd_t knbhd = d->preconstructed_nbhds[rid];
	int t_nbhd_size = knbhd.n;

	if (observed_only || t_nbhd_size == 0 || knbhd.nbhd == NULL) {
		char *pmem = mempool_nalloc(mp, nbhd_alloc_size(1), 16);
		hmm_setup_nbhd_ptrs(first_nbhd, 1, pmem);

		// no preconstructed neighborhood, simply use the observed kmer
		first_nbhd->size = 1;
		first_nbhd->alphas[0] = 0.0;
		first_nbhd->states_sorted_sfx[0] = obs_1st_kid;
		first_nbhd->kmer_ptrs[0] = ref_kmer; 
		first_nbhd->kmer_trans_prob[0] = ref_kmer->transition_p;
		first_nbhd->edit_distances[0] = 0; 

		bitarray_set_pwr2(first_nbhd->ba_kmer_trans_flag, 0, 
				ref_kmer->trans_flag, 5);

		return 1;
	}

	char *pmem = mempool_nalloc(mp, nbhd_alloc_size(t_nbhd_size), 16);
	hmm_setup_nbhd_ptrs(first_nbhd, t_nbhd_size, pmem);

	first_nbhd->size = t_nbhd_size;
	for (int i = 0; i < t_nbhd_size; i++) {
		register kmer_t *nb_k = knbhd.nbhd[i];
		kbits_t kid = nb_k->id;
		first_nbhd->alphas[i] = 0.0;
		first_nbhd->states_sorted_sfx[i] = kid;
		first_nbhd->kmer_ptrs[i] = nb_k; 
		first_nbhd->kmer_trans_prob[i] = nb_k->transition_p;
		first_nbhd->edit_distances[i] = kmer_n_hamming_dist(obs_1st_kid, kid, obs_n_flag);

		bitarray_set_pwr2(first_nbhd->ba_kmer_trans_flag, i, 
				nb_k->trans_flag, 5);
	}

	return t_nbhd_size;

	/*
	first_nbhd->size = t_nbhd_size;
	first_nbhd->alphas[0] = 0.0;//ref_kmer->init_p;
	first_nbhd->states_sorted_sfx[0] = ref_kmer->id;
	first_nbhd->kmer_ptrs[0] = ref_kmer;
	first_nbhd->kmer_trans_prob[0] = ref_kmer->transition_p;
	first_nbhd->edit_distances[0] = 0;

	bitarray_set_pwr2(first_nbhd->ba_kmer_trans_flag, 0, 
				ref_kmer->trans_flag, 5);

	return t_nbhd_size;
	*/
}

static void hmm_update_expected_vals(data *d, state_nbhd_t *curr_nbhd, 
		int t, kbits_t xt, int yt, double *t_exp_counts, double *t_exp_trans, 
		double inv_sum_exp_count, int kmer_len)
{
	int qmax = d->opt->max_qual_score + 1;
	int shift = (kmer_len - 1) << 1;
	int t_nbhd_size = curr_nbhd->size;

	double *kmer_exp_trans = t_exp_trans;
	for (int i = 0; i < t_nbhd_size; i++) {
		kmer_t *pk = curr_nbhd->kmer_ptrs[i];
		kbits_t kid = pk->id;
		register kbits_t tflag = pk->trans_flag;

		int is_triplet = 0;

		if (kid >> 62 == 0) {
			kbits_t _b1 = (kid >> shift) & 3;
			kbits_t _b2 = (kid >> (shift-2)) & 3;
			kbits_t _b3 = (kid >> (shift-4)) & 3;

			if (_b1 == _b2 && _b2 == _b3) 
				is_triplet = 1;
		}

		double norm_exp_count = EXP(PRODUCT(t_exp_counts[i], 
						inv_sum_exp_count));
		//pk->exp_count += norm_exp_count;

		kbits_t st = kbits_last_base(kid, kmer_len);
		// quality score emission
		if (kid >> 62) {
			atomic_add_dbl(&d->qual_emit_exp[qmax * QCTX_KP1MER + yt],
					norm_exp_count);
			//d->qual_emit_exp[qmax * QCTX_KP1MER + yt] += norm_exp_count;
		}
		else {
			double ins_prop = curr_nbhd->ins_prop[i];

			// is this a insertion kmer?
			atomic_add_dbl(&d->qual_emit_exp[qmax * QCTX_SELF_TRANS + yt],
					norm_exp_count * ins_prop);
			// normal transition kmer
			atomic_add_dbl(&d->qual_emit_exp[qmax * (st != xt) + yt],
					norm_exp_count * (1.0-ins_prop));
					
			/*
			d->qual_emit_exp[qmax * QCTX_SELF_TRANS + yt] += 
				norm_exp_count * ins_prop;
			d->qual_emit_exp[qmax * (st != xt) + yt] += 
				norm_exp_count * (1-ins_prop);
			*/
		}

		// base emission
		atomic_add_dbl(&d->base_emit_exp[st * 5 + xt],
				norm_exp_count);
		//d->base_emit_exp[st * 5 + xt] += norm_exp_count;

		int n = 0;
		for (int j = 0; j <= 20; j++) {
			if ((tflag & (1 << j)) == 0) continue;

			double exp_trans = EXP(PRODUCT(kmer_exp_trans[n], 
							inv_sum_exp_count));

			atomic_add_dbl(&d->all_sum_exp_trans, exp_trans);

			if (is_triplet) {
				atomic_add_dbl(&triplet_del_exp[((kid >> shift) & 3) * 2],
						exp_trans);
			}
			//d->all_sum_exp_trans += exp_trans;

			if (j < 4) {
				int dns_tidx = _mm_popcnt_u64(tflag & 
						((1UL << j) - 1));
				//selfflag=0;
				atomic_add_dbl(&pk->exp_trans_p[dns_tidx], exp_trans);
				//pk->exp_trans_p[dns_tidx] += exp_trans;
			}
			else if (j == 20) {
				atomic_add_dbl(&d->ins_sum_exp_trans, exp_trans);
				/*
				if (selfflag==0){
				d->ins_sum_exp_trans += exp_trans;
				}
				else {
				d->ins_sum_exp_trans_self += exp_trans;
				}
				selfflag=1;
				*/
				// use exp_count to record the expected # of (self-)transitions
				// for each kmer.
				// pk->exp_count += exp_trans;
			}	
			else {
				register kbits_t base1 = (j-4) & 3;
				register kbits_t base2 = (j-4) >> 2;
				kbits_t kid2 = kbits_raw_suffix(pk->id) | (base1 << shift);
				kmer_t *pk2 = dict_find(d->kmer_dict, 
						kbits_cast_to_ptr(kid2))->value;

				int dns_tidx1 = _mm_popcnt_u64(tflag & 
						((1UL << base1) - 1));
				int dns_tidx2 = _mm_popcnt_u64(pk2->trans_flag & 
						((1UL << base2) - 1));

				atomic_add_dbl(&pk->exp_trans_p[dns_tidx1], exp_trans);
				atomic_add_dbl(&pk2->exp_trans_p[dns_tidx2], exp_trans);
				atomic_add_dbl(&d->del_sum_exp_trans, exp_trans);

				if (is_triplet) {
					// triplet in the kmer suffix
					atomic_add_dbl(&triplet_del_exp[((kid >> shift) & 3) * 2 + 1],
							exp_trans);
				}

				/*
				pk->exp_trans_p[dns_tidx1] += exp_trans;
				pk2->exp_trans_p[dns_tidx2] += exp_trans;

				d->del_sum_exp_trans += exp_trans;
				*/
			}

			++n;
		}

		kmer_exp_trans += (pk->n_trans + 1);
	}
}

#if 0
static double hmm_compute_kmer_emission(data *d, kbits_t obs_kid, 
		kbits_t obs_n_flag, kmer_t *pk, int kmer_len, const char *conv_q)
{
	char emit_norm_flag[kmer_len];
	double emit_dens = 0.0;

	memset(emit_norm_flag, 0, kmer_len);

	/* compute emission renormalization flag */
	kmer_t **k_nbhd = pk->data;
	for (int n = 0; n < pk->n_neighbors; n++) {
		kbits_t nb_kid = (k_nbhd[n])->id;
		for (int i = 0; i < kmer_len; i++) {
			emit_norm_flag[i] |= (1 << ((nb_kid >> (i << 1)) & 3));
		}
	}

	kbits_t kid = pk->id;
	for (int t = 0; t < kmer_len; t++, kid >>= 2, obs_kid >>= 2) {
		register kbits_t xt = obs_kid & 3, st = kid & 3;
		double *t_base_emit_p = &d->base_emission_p[t * d->base_ctx_radix + 
			st * 5];

		emit_dens += hmm_normalize_emission(t_base_emit_p, 
				(kbits_t) emit_norm_flag[t], xt);
		emit_dens += d->qual_emission_p[t * d->qual_ctx_radix +
			((int) conv_q[t] << 1) + (xt == st)];
	}

	return emit_dens;
}
#endif

static inline void hmm_setup_nbhd_ptrs(state_nbhd_t *nbhd, int n, 
		char *pmem)
{
	memset(pmem, 0, nbhd_alloc_size(n));

	/* FIXME: consider lookup table and struct initialization */
	nbhd->kmer_ptrs = (kmer_t **) pmem;
	nbhd->kmer_trans_prob = (double **) (pmem + n * sizeof(*nbhd->kmer_ptrs));

	nbhd->ba_kmer_trans_flag = (kbits_t *) (
			(char *) nbhd->kmer_trans_prob + n * 
			sizeof(*nbhd->kmer_trans_prob));
	//nbhd->ba_distinct_sfx = nbhd->ba_kmer_trans_flag + BITS_TO_U64(n << 2);
	//nbhd->ba_hamming_dist = nbhd->ba_distinct_sfx + BITS_TO_U64(n << 1);
	nbhd->ba_pfx_sfx_flags = nbhd->ba_kmer_trans_flag + BITS_TO_U64(n << 5);

	nbhd->suffixes_order = (int *) (nbhd->ba_pfx_sfx_flags + 
			BITS_TO_U64(n << 6));
	nbhd->edit_distances = (uint64_t *) (nbhd->suffixes_order + n);
	nbhd->states_sorted_pfx = (kbits_t *) (nbhd->edit_distances + n);
	nbhd->states_sorted_sfx = nbhd->states_sorted_pfx + n;

	nbhd->alphas = (double *) (nbhd->states_sorted_sfx + n);
	nbhd->betas = nbhd->alphas + n;
	nbhd->emit_dens = nbhd->betas + n;
	nbhd->ins_prop = nbhd->emit_dens + n;
}

inline char *convert_quality_scores(const char *qscore, int len, 
		int offset, mempool_t *mp)
{
#ifdef __SSE2__
	/* FIXME: cache this result maybe */
	int i, w_m128 = (len >> 4) + ((len & 0xF) ? 1 : 0);
	char *conv_qscore = mempool_nalloc(mp, sizeof(*conv_qscore) * len, 16);
	__m128i *src_ptr = (__m128i *) qscore;
	__m128i *dst_ptr = (__m128i *) conv_qscore;
	__m128i qual_offset = _mm_set1_epi8(offset);
	__m128i xmm;

	for (i = 0; i < w_m128; i++){
		xmm = _mm_loadu_si128(src_ptr);
		xmm = _mm_sub_epi8(xmm, qual_offset);
		_mm_store_si128(dst_ptr, xmm);

		++dst_ptr;
		++src_ptr;
	}
#endif

	return conv_qscore;
}

inline static int *compute_qual_contexts(const char *conv_q, int rlen, 
		int klen, int ctxlen, int qmax, mempool_t *mp)
{
	int *qual_ctx = mempool_alloc(mp, sizeof(*qual_ctx) * rlen);
	for (int i = 0; i < klen; i++) {
		qual_ctx[i] = 0;
	}

	for (int i = klen; i < rlen; i++) {
		int qctx_idx = 0, radix = 1;
		for (int j = ctxlen - 1; j > 0; j--) {
			qctx_idx += radix * conv_q[i - j];
			radix *= qmax;
		}
		qual_ctx[i] = qctx_idx;
	}
	return qual_ctx;
}

static int kmer_cmp_suffix(const void *a, const void *b)
{
	const candidate_t *ca = a;
	const candidate_t *cb = b;
	kbits_t ra = kbits_raw_suffix(ca->id);
	kbits_t rb = kbits_raw_suffix(cb->id);

	if (ra < rb) {
		return -1;
	}
	else if (ra == rb) {
		return (ca->id < cb->id) ? -1 : 1;
	}
	else {
		return 1;
	}
}

#if 0
struct kmer_s {
	double expected_count;
	uint64_t trans_flag;
	/*
	 *      |  next step  | |kmer->A->ACGT|     |kmer->G->ACGT|
	 * aaaa bcc/bcc/bcc/bcc bcc/bcc/bcc/bcc ... bcc/bcc/bcc/bcc
	 *
	 * aaaa: 4 bit mask for valid transitions;
	 *		 set to 1 if transition is admissible, 0 otherwise
	 *
	 * bcc:  b == 1 if transition is unique, cc then is the downstream base.
	 *       bcc == 000 if no transition is allowed, 
	 *       011 if multiple transitions admissible.
	 */
}
#endif
