#ifndef __PREMIER_H__
#define __PREMIER_H__

#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "oadict.h"
#include "mempool.h"

/* 20 = 4 (true base) * 5 (observed base, ACTGN) */
#define PMR_BASE_EMIT_COMBS 20
#define PMR_UNIF_BASE_EMIT_P -1.6094379124341002818 
#define PMR_LINE_BUFFER_SIZE 1024
#define PMR_READ_BUFFER_SIZE 1073741824 //(1<<30) 

/* --- some precompiled switches --- */

/* slightly faster, but not numerically stable,
 * especially for *LONG READS* */
#define PMR_USE_LOG_ADD 0

/* NOTE: not in regular ACGT order for computation convenience */
enum bases {BASE_A, BASE_C, BASE_T, BASE_G, BASE_N};
enum qual_emit_ctxs { QCTX_KMER_EQUAL, QCTX_KMER_NEQ, QCTX_KP1MER,
	QCTX_SELF_TRANS };
enum unlimit_ins_mode { INS_MODE_OFF, INS_MODE_EM = 0x1, INS_MODE_VITERBI = 0x2,
	INS_MODE_BOTH = 0x3 };

typedef uint64_t kbits_t;
typedef struct kmer_s kmer_t;
typedef struct state_nbhd_s state_nbhd_t;
typedef struct kmer_nbhd_s kmer_nbhd_t;
typedef struct _options options;

typedef union {
	double dbl;
	uint64_t u64;
} dbl_t;

typedef struct {
	int size;
	int used;
	kmer_t **kmer_refs; /*< pointer to these reference kmers */
	/* hash_t *raw_hashes; */
} nbhd_t;

typedef struct {
	kbits_t masked_id;
	kbits_t pmask;	/* position_mask */
} nbhd_key_t;

typedef struct {
	int length;	/*< short read length */
	size_t id;	/* read identifier */

	char *identifier;
	char *sequence;	/*< literal read sequence */
	char *qscore;	/*< quality score associated with a read */
} read_t;

struct kmer_s {
	double count;
	double exp_count;
	double init_p;
	double *transition_p;
	double *exp_trans_p;

	kbits_t id;
	uint32_t trans_flag : 21; /* 4 + 16 */
	/* # of nonzero transitions to A/C/G/T (kmer) */
	uint32_t n_ktrans: 4;
	/* # of nonzero transitions out of 21 */ 
	uint32_t n_trans: 5;
	uint32_t dirty : 1;
	uint32_t label : 1; // the labeling of a particular kmer
	
	int alloc_nbhd_size;
	int n_neighbors;
	void *data;
	/*
	uint32_t unique_trans : 1;
	uint32_t alloc_nbhd_size : 4;
	uint32_t emit_norm_flag : 4;
	uint32_t n_neighbors : 19;
	*/

	//void *data;
};

/* FIXME: consider rearrange memory to improve locality / optimize
 * cache usage */
struct state_nbhd_s {
	int size;
	int n_uniq_sfx;
	kmer_t **kmer_ptrs;
	double **kmer_trans_prob;

	uint64_t *ba_kmer_trans_flag;
	//uint64_t *ba_distinct_sfx;
	//uint64_t *ba_hamming_dist;
	uint64_t *ba_pfx_sfx_flags;

	int *suffixes_order;
	uint64_t *edit_distances;
	kbits_t *states_sorted_pfx;
	kbits_t *states_sorted_sfx;

	double *alphas;
	double *betas;
	double *emit_dens;
	// since we aggregate the two version of the same nominal kmers,
	// (i.e. normal kmers and insertion copies) into one single representation,
	// we need to record the relative contribution of the forward-algorithm 
	// inductive variable, from the insertion copy. 
	double *ins_prop;
};

typedef struct {
	dict *kmer_dict;
	dict *nbhd_dict;

	options *opt;

	int max_nbhd_size;
	int max_read_length;

	int n_reads;
	int n_reads_cannot_decode;
	uint32_t n_reads_zero_likelihood;

	size_t n_qual_context;
	size_t qual_ctx_radix;
	size_t n_base_context;
	size_t base_ctx_radix;

	/* --- Baum-Welch --- */
	double loglik;
	double total_kmer_exp_counts;
	
	/* for computing pd and pi,
	 * which are probability of deletion/insertion errors */
	double ins_error_rate;
	double ins_error_rate_self;	
	double del_error_rate;
	double error_free_p;
	double emit_error_rate;

	double ins_sum_exp_trans;
	double ins_sum_exp_trans_self;
	double del_sum_exp_trans;
	double all_sum_exp_trans;

	/* - transition - */
	double *transition_p;

	/* - emission - */
	double *base_emission_p;
	double *qual_emission_p;
	double *base_emit_exp;
	double *qual_emit_exp;

	double *qual_err_rate;

	int *read_start_pos;

	/* --- Precomputed tables --- */
	char **n_choose_k;	/*< binomial coefficients (n choose k) */
	char **position_perms;	/*< permutation of all {k \choose d} erroneous positions */
	uint64_t *position_masks;	/*< masks for {k \choose d} erroneous positions. */
	uint64_t *error_perms;
	kbits_t *transition_masks;
	double *phred_to_prob; /*< table that converts Phred quality score to error probability. */ 

	/* --- output files --- */
	FILE *log_fp;
	FILE *errors_fp;
	FILE *params_fp;

	kmer_nbhd_t *preconstructed_nbhds;

	mempool_t *mp;
	mempool_t *tmpmp;
} data;

struct _options {
	int kmer_length;
	int num_threads;

	int preconstructed_nbhd_hd;
	int max_edit_dist; 

	int pagesize;

	int qual_score_offset;
	int max_qual_score;

	int discarding_threshold;

	double tol;
	double error_rate;

	uint32_t output_params : 1;
	uint32_t random_param_init : 1; 
	uint32_t disable_penalized_em : 1;
	uint32_t enable_params_output : 1;
	uint32_t filename_append_timeinfo : 1;
	uint32_t filename_append_params : 1;
	uint32_t penalize_init_params : 1;
	uint32_t viterbi_revcomp : 1;
	// as a "kludge" to deal with first neighborhood issues
	uint32_t adjust_read_start_pos : 1;
	uint32_t skip_em : 1;
	uint32_t unlimited_insertions : 2;
	uint32_t dummy : 23; 

	/* penalized EM */
	double penalty_eta;
	double penalty_gamma;

	const char *fastq_file;
	char *output_prefix;
	char *params_data;
};

typedef struct {
	double kappa;
	double eta;
	double *xi;
	int32_t signs;
	int n_xi;
} mstep_data_t ;

struct kmer_nbhd_s {
	int n;
	kmer_t **nbhd;
};

typedef struct {
	kbits_t id;
	kbits_t mut_flag;
} mut_kmer_t;

void hmm_build_1st_nbhd(data *d, read_t *read, void *fdata);

#endif
