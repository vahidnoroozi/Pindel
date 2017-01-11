#define _GNU_SOURCE

#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <getopt.h>
#include <time.h>
#include <signal.h>

#include "premier.h"
#include "hashfunc.h"
#include "read.h"
#include "kmer.h"
#include "data.h"
#include "iodata.h"
#include "em.h"
#include "numeric.h"
#include "nbhd.h"

static int make_options(options **opt, int argc, char **argv);
static int make_data(data **d, options *opt);
static int destroy_data(data *d);
static void print_usage(const char *basename);

static char *pmr_usage = ""
	"usage: %s [options] <input_fastq>\n\n"
	"Indel-aware NGS error-corrector.\n\n"
	"positional arguments:\n"
	"  input_fastq   Input FASTQ file for error-correction.\n\n"
	"optional arguments:\n"
	"  -h            show this help message and exit\n"
	"  -k K          Kmer length.\n"
	"  -d D          Maximum edit distance allowed.\n"
	"  -E ETA        Penalty threshold parameter.\n"
	"  -G GAMMA      Penalty shape parameter [default: 1e-20].\n"
	"  -o PREFIX     Prefix for the output FASTQ file.\n"
	"  -R RELTOL     Threshold to declare convergence of the EM algorithm, measured\n"
	"                by relative increment of the penalized objective function.\n"
	"  -t N_THREADS  Number of threads for parallel execution [default: 1]\n"
	"  -i [options]  Toggle 'unlimited insertions' mode. [default: off]\n"
	"                -i em:      enable this mode for the EM algorithm;\n"
	"                -i viterbi: enable this mode for the Viterbi algorithm;\n"
	"                -i both:    enable this mode for both algorithms.\n"
	"  -e            Skip the EM algorithm all together. [default: off]\n";

int main(int argc, char **argv)
{
	data *d = NULL;
	options *opt = NULL;

	struct timeval ts, te;
	time_t rt;

	if (make_options(&opt, argc, argv))
		exit(EXIT_FAILURE);

	/* record timestamp */
	gettimeofday(&ts, NULL);

	if (make_data(&d, opt)){
		/* error processing here */
		exit(EXIT_FAILURE);
	}

	/*
	fprintf(stderr, "Allocated slots in dict: %u\n",
			d->kmer_dict->sizemask - 1);
	*/

	/* data prepared, proceed to em */
	EM(d, &ts);

	/* timing info */
	time(&rt);
	gettimeofday(&te, NULL);

	/*
	if (d->errors_fp != NULL){
		fprintf(d->errors_fp, "# HMMEC finished at %s\n", ctime(&rt));
		fprintf(d->errors_fp, "# Total seconds elapsed: %.3f\n",
				(0.0 + te.tv_sec + (double) te.tv_usec/1e6) - 
				(0.0 + ts.tv_sec + (double) ts.tv_usec/1e6));
		fclose(d->errors_fp);
	}
	*/

	/*
	if (d->param_fp != NULL)
		fclose(d->param_fp);
	if (d->ll_fp != NULL)
		fclose(d->ll_fp);
	*/

	/* cleanups */
	destroy_data(d);

	return 0;
}

int make_data(data **d, options *opt)
{
	int kmer_len = opt->kmer_length;
	int dmax = opt->preconstructed_nbhd_hd;

	data *pd;

	*d = (data *) malloc(sizeof(**d));
	if (*d == NULL){
		return message(stderr, __FILE__, __LINE__,
				ERROR_MSG, MEMORY_ALLOCATION, NULL);
	}

	pd = *d;
	pd->opt = opt;

	/* initialize some default values */
	pd->max_read_length = 0;
	pd->total_kmer_exp_counts = 0.0;

	pd->ins_error_rate = LOG(0.34);
	//pd->ins_error_rate_self=LOG(0.04);
	pd->del_error_rate = LOG(0.34);
	pd->error_free_p = LOG(1-0.34-0.34);
	pd->emit_error_rate = 0.01;

	pd->errors_fp = NULL;
	pd->params_fp = NULL;

	pd->n_reads = 0;

	/* FIXME: consider an exclusive function for opening files */
	if (pd->opt->output_prefix != NULL) {

		char *out_fname = io_cat_output_filename(pd->opt, "", "fastq");	
		if (out_fname != NULL) {
			pd->errors_fp = fopen(out_fname, "w+");
			free(out_fname);
		}

		/*
		out_fname = io_cat_output_filename(pd->opt, "LOG", NULL);
		if (out_fname != NULL) {
			pd->log_fp = fopen(out_fname, "w+");
			free(out_fname);
		}
		*/

	}

#if 0
	if (pd->opt->logfile_prefix != NULL) {
		ec_out_fname = generate_log_filename(pd->opt, "EC");
		if (ec_out_fname != NULL){
			pd->errors_fp = fopen(ec_out_fname, "w+");
			free(ec_out_fname);
		}

		/* output the absolute path to the FASTQ file */
		char fastq_realpath[PATH_MAX + 1];
		realpath(pd->opt->fastq_file, fastq_realpath);

		/* generate some useful information in the output */
		time(&rt);
		fprintf(pd->errors_fp, "# HMMEC starts at %s\n", ctime(&rt));
		fprintf(pd->errors_fp, "# Sequencing FASTQ datafile: %s\n", 
				fastq_realpath);	
		fprintf(pd->errors_fp, "# Kmer length: %d\n", 
				kmer_len);
		fprintf(pd->errors_fp, "# Max substitution errors: %d\n",
				pd->opt->max_substitution_errors);

		fprintf(pd->errors_fp, "# Error rates initiated to truth: %d\n",
				INIT_ERROR_RATE_TO_TRUTH);

		if (pd->opt->enable_param_output){
			param_out_fname = generate_log_filename(pd->opt, "PARAM");
			if (param_out_fname != NULL){
				pd->param_fp = fopen(param_out_fname, "w+");
				free(param_out_fname);
			}
		}

		if (pd->opt->enable_ll_output) {
			param_out_fname = generate_log_filename(pd->opt, "LL");
			if (param_out_fname != NULL) {
				pd->ll_fp = fopen(param_out_fname, "w+");
				free(param_out_fname);
			}
		}

		if (pd->opt->enable_err_rate_output) {
			param_out_fname = generate_log_filename(pd->opt, "ER");
			if (param_out_fname != NULL) {
				pd->err_rate_fp = fopen(param_out_fname, "w+");
				free(param_out_fname);
			}
		}

		if (pd->opt->dump_viterbi) {
			param_out_fname = generate_log_filename(pd->opt, "VD");
			if (param_out_fname != NULL){
				pd->vtrace_fp = fopen(param_out_fname, "wb");
				free(param_out_fname);
			}
		}
	}
#endif

	INFO("Creating memory pools...\n");
	pd->mp = mempool_create(MEMPOOL_DEFAULT_SIZE);

#if 0
	INFO("Precomputing lookup tables...\n");
	/* --- precompute lookup tables --- */
	int err;
	if ((err = compute_binomial_coefficients(pd)) != NO_ERROR) {
		exit(err);
	}
#endif

	/* use open addressing to resolve conflict (implemented in oadict.c) */
	/* FIXME: for large k where k-spectrum is sparse, 
	 * pow(4, kmer_len-2) is way too more aggresive */
	pd->kmer_dict = dict_create(pow(4, kmer_len-(kmer_len>>1)),
			dict_uint64_hash_raw, dict_compare_identity);

#if 0
	/* some precomputed tables for error properties */
	if ((err = compute_error_positions(pd)) != NO_ERROR)
		exit(err);

	if ((err = compute_all_errors(pd)) != NO_ERROR)
		exit(err);

	/* update some constants */
	pd->max_nbhd_size = 0;
	for (int ne = 1; ne <= pd->opt->preconstructed_nbhd_hd; ne++){
		pd->max_nbhd_size += pd->n_choose_k[kmer_len][ne] * 
			(int) pow(3, ne);
	}
#endif

	INFO("Constructing k-spectrum for k=%d...\n", pd->opt->kmer_length);

	/* Iterate over all available reads provided by the FASTQ data file,
	 * with callback function 'read_parse'.
	 *
	 * After this step, we will fill the (*d)->kmer_dict buckets with
	 * observed kmers (k-spectrum).
	 */
	// iterate_reads_fastq(*d, (*d)->opt->fastq_file, dummy_callback);
	io_iterate_reads_fastq(pd, 1, false, read_parse, NULL);

	INFO("Number of total k-mers and (k+1)-mers: %zu\n", (size_t) pd->kmer_dict->used);
	size_t n_kmers = 0, n_kp1mers = 0, n_discarded = 0;
	
	INFO("Allocating memory for transition parameters...\n");

	pd->read_start_pos = calloc(pd->n_reads, sizeof(*pd->read_start_pos));

	dictentry *de = pd->kmer_dict->table;
	int i = pd->kmer_dict->used;
	for (; i > 0; de++) {
		if (de->value == NULL) continue;

		i--;
		kmer_t *pk = de->value;
		/* filtering by discarding kmers with multiplity below some threshold */
		/*
		if (pk->count < pd->opt->discarding_threshold) {
			de->value = NULL;
			de->key = NULL;
			pd->kmer_dict->used--;

			n_discarded++;
			continue;
		}
		*/
		if (pk->id >> 62) continue;

		int n_transitions = _mm_popcnt_u64(pk->trans_flag);
		int n_ktrans = _mm_popcnt_u64(pk->trans_flag & 15);
		pk->n_trans = n_transitions;
		pk->n_ktrans = n_ktrans;

		pk->exp_trans_p = mempool_alloc(pd->mp, n_ktrans * sizeof(*pk->exp_trans_p));
		memset(pk->exp_trans_p, 0, n_ktrans * sizeof(*pk->exp_trans_p));
		pk->transition_p = mempool_alloc(pd->mp, n_ktrans * sizeof(*pk->transition_p));
		memset(pk->transition_p, 0, n_ktrans * sizeof(*pk->transition_p));
	}

	INFO("Setting (k+1)-mers' transition flags...\n");
	de = pd->kmer_dict->table;
	i = pd->kmer_dict->used;
	for (; i > 0; de++) {
		if (de->value == NULL) continue;

		i--;
		kmer_t *pkp1 = de->value;

		if ((pkp1->id >> 62) == 0) continue;

		// inherit all properties from the corresponding kmer.
		kbits_t _kid = (pkp1->id & KBITS_KLEN_REMOVAL_FLAG) >> 2;
		kmer_t *pk = dict_find(pd->kmer_dict, kbits_cast_to_ptr(_kid))->value;
		// do not allow self-transition
		pkp1->trans_flag = pk->trans_flag & ((1U << 20) - 1);
		pkp1->n_trans = pk->n_trans;
		pkp1->n_ktrans = pk->n_ktrans;

		// do NOT introduce new FREE parameters.
		pkp1->exp_trans_p = pk->exp_trans_p;
		pkp1->transition_p = pk->transition_p;
	}

	INFO("Tallying transition incidence matrix for the k-spectrum...\n");
	io_iterate_reads_fastq(pd, 1, true, read_parse_trans_count, NULL);

	de = pd->kmer_dict->table;
	i = pd->kmer_dict->used;

	// label the kmer/transitions to exploit dependence inherent to the 
	// DNA double-strandedness.
	INFO("Adjusting the labeling (strandedness) of kmers...\n");
	//pd->total_kmer_exp_counts = 0.0;
	//adjust_kmer_strandedness(pd);
	adjust_kmer_labeling(pd);

	INFO("Initializing the initial state probabilities...\n");
	/* initialize the initial distribution probability. */
	for (; i > 0; de++) {
		if (de->value == NULL) continue;

		i--; 
		kmer_t *pk = de->value;
		kbits_t ktype = pk->id >> 62;

		if (ktype || pk->dirty == 0) continue;

		kmer_t *rc_pk = dict_find(pd->kmer_dict,
				kbits_cast_to_ptr(kmer_reverse_complement(pk->id, kmer_len)))->value;

		double _init_p = LOG((pk->exp_count + rc_pk->exp_count) /
				pd->total_kmer_exp_counts);

		pk->init_p = _init_p;
		rc_pk->init_p = _init_p;

		pk->dirty = 0;
		rc_pk->dirty = 0;

		// reset expected counts
		pk->exp_count = 0.0;
		rc_pk->exp_count = 0.0;

		/* for kmer only */
		pk->trans_flag |= (1UL << 20);
		rc_pk->trans_flag |= (1UL << 20);
	}

	INFO("Initializing the transition probabilities...\n");
	mstep_transition(pd);

	//pd->total_kmer_exp_counts = LOG(pd->total_kmer_exp_counts);
#if 0
	/* compute masks for all 16 possible combination of transitions */
	pd->transition_masks = mempool_alloc(pd->mp, sizeof(kbits_t) << 6);
	memset(pd->transition_masks, 0xFF, sizeof(kbits_t) << 6);
	for (int t = 0; t < 16; t++) {
		int tsl2 = t << 2;
		for (int i = 0; i < 4; i++) {
			if ((t >> i) & 1) {
				pd->transition_masks[tsl2 + i] = 0; 
			}
		}
	}

	INFO("Computing 1-neighborhood given k-spectrum...\n");
	compute_kmer_hamming_nbhds(pd);
#endif

	if (pd->opt->qual_score_offset < 64) 
		pd->opt->qual_score_offset = 33;
	else 
		pd->opt->qual_score_offset = 64;

	pd->opt->max_qual_score -= pd->opt->qual_score_offset;

	INFO("Quality score format: Phred+%d, range : (0, %d)\n", 
			pd->opt->qual_score_offset, pd->opt->max_qual_score);

	/* Allocate memory for storing and computing expected emission
	 * densities */
	INFO("Allocating memory for HMM emission distribution parameters...\n");

	pd->base_ctx_radix = 5 * 4; 
	pd->base_emission_p = mempool_alloc(pd->mp, 
			pd->base_ctx_radix * sizeof(*pd->base_emission_p));
	pd->base_emit_exp = mempool_alloc(pd->mp, 
			pd->base_ctx_radix * sizeof(*pd->base_emit_exp));

	// 4 different scenarios.
	pd->qual_ctx_radix = ((pd->opt->max_qual_score + 1) * 4);
	pd->qual_emission_p = mempool_alloc(pd->mp,
			pd->qual_ctx_radix * sizeof(*pd->qual_emission_p));
	pd->qual_emit_exp = mempool_alloc(pd->mp,
			pd->qual_ctx_radix * sizeof(*pd->qual_emit_exp));

	pd->preconstructed_nbhds = calloc(pd->n_reads, sizeof(kmer_nbhd_t));

	INFO("Creating first neighborhoods (substitution errors only)...\n");
	io_iterate_reads_fastq(pd, 1, true, hmm_build_1st_nbhd, NULL);

	return NO_ERROR;
}

int destroy_data(data *d)
{
	/* free precomputed tables */
	/*
	int n = d->opt->kmer_length + 1;
	for (int i = 0; i < n; i++){
		free(d->n_choose_k[i]);
	}

	free(d->n_choose_k);

	n = d->opt->preconstructed_nbhd_hd; 
	for (int i = 0; i < n; i++){
		free(d->position_perms[i]);
	}

	free(d->position_perms);

	free(d->error_perms);
	free(d->position_masks);
	free(d->phred_to_prob);
	*/
 
	/* destroy all dictionaries */
	dict_destroy(d->kmer_dict);

	/* release all memory allocated for memory pools */
	mempool_destroy(d->mp);

	free(d->opt);
	free(d);

	return 0;
}

int make_options(options **opt, int argc, char **argv)
{
	/*
	static struct option hmmec_options[] = {
		{"lambda", required_argument, 0, 'L'},
		{"gamma", required_argument,  0, 'G'},
		{"rel-tol", required_argument, 0, 'R'},
		{"base-emit-context", required_argument, 0, 0},
		{"qual-emit-context", required_argument, 0, 0},
		{"mph-file", required_argument, 0, 'M'},
	};
	*/

	int rand_seed;
	size_t arg_len;
	*opt = (options *) malloc(sizeof(**opt));

	options *popt = *opt;

	popt->num_threads = 1;
	popt->kmer_length = 0;
	popt->max_edit_dist = 4;
	popt->preconstructed_nbhd_hd = 1;
	popt->pagesize = getpagesize();
	popt->qual_score_offset = 255;
	popt->max_qual_score = 0;
	popt->random_param_init = 0;
	popt->filename_append_timeinfo = 1;
	popt->filename_append_params = 1;
	popt->discarding_threshold = 4;
	popt->viterbi_revcomp = 1;
	popt->adjust_read_start_pos = 1; // EM only (for higher coverage)?
	popt->unlimited_insertions = 0;
	popt->skip_em = 0;

	popt->output_prefix = NULL;
	popt->output_params = 0;

	popt->disable_penalized_em = 0;
	popt->penalize_init_params = 0;
	popt->error_rate = 0.02;
	popt->tol = 1e-5;

	popt->penalty_gamma = 1e-18;
	/* FIXME: automatic way to determine lambda ? */
	popt->penalty_eta = 100;

	/* parse user specified options */
	char *ptr_optarg;
	char c;
	while ((c = getopt(argc, argv, "k:d:D:o:E:G:R:P:s:t:i:NVpeh")) != -1) {
		switch (c) {
			/* kmer-length */
			case 'k':
				popt->kmer_length = strtol(optarg, &ptr_optarg, 0);
				if (ptr_optarg == optarg){
					message(stderr, __FILE__, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"-k %s\n", optarg);

					exit(INVALID_CMD_ARGUMENT);
				}

				/*
				if (popt->kmer_length <= 0 || popt->kmer_length > 31){

					exit(INVALID_CMD_ARGUMENT);
				}
				*/

				break;
			/* maximum substitutional errors */
			case 'd':
				popt->max_edit_dist = strtol(optarg, &ptr_optarg, 0);

				if (ptr_optarg == optarg){
					message(stderr, __FILE__, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"-d %s\n", optarg);

					exit(INVALID_CMD_ARGUMENT);
				}

				/*
				if (popt->max_substitution_errors < 1 || 
						popt->max_substitution_errors > 16){
					message(stderr, file_name, fxn_name, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"please provide valid error distance [1, 6]\n");

					exit(INVALID_CMD_ARGUMENT);
				}
				*/

				break;

			case 't':
				popt->num_threads = strtol(optarg, &ptr_optarg, 0);
				break;
			case 'D':
				popt->discarding_threshold = strtol(optarg, &ptr_optarg, 0);
				break;
			case 'V':
				popt->viterbi_revcomp = 0;
				break;

			case 'E':
				popt->penalty_eta = strtod(optarg, &ptr_optarg);
				/*
				if (ptr_optarg == optarg || popt->lambda < 0) {
					message(stderr, file_name, fxn_name, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"-L %s\n", optarg);

					exit(INVALID_CMD_ARGUMENT);
				}
				*/

				break;
			case 'G':
				popt->penalty_gamma = strtod(optarg, &ptr_optarg);
				/*
				if (ret_optarg_ptr == optarg || popt->kappa < 0) {
					message(stderr, file_name, fxn_name, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"-K %s\n", optarg);

					exit(INVALID_CMD_ARGUMENT);
				}
				*/
				break;

			case 'R':
				popt->tol = strtod(optarg, &ptr_optarg);
				break;

			/* output all trained parameters to local file PREFIX.param */
			case 'p':
				popt->output_params = 1;
				break;

			/* load parameters from file */
			case 'P':
				arg_len = strlen(optarg);
				popt->params_data = calloc(arg_len + 1, sizeof(char));
				strncpy(popt->params_data, optarg, arg_len);
				break;

			case 'N':
				popt->disable_penalized_em = 1;
				break;

			/* prefix for log files */
			case 'o':
				arg_len = strlen(optarg);
				popt->output_prefix = calloc(arg_len + 1, sizeof(char));
				strncpy(popt->output_prefix, optarg, arg_len);

				break;
			case 'e':
				popt->skip_em = 1;
				break;
			case 'i':
				if (strncmp(optarg, "em", 2) == 0) {
					popt->unlimited_insertions = 1;
				}
				else if (strncmp(optarg, "viterbi", 7) == 0) {
					popt->unlimited_insertions = 2;
				}

				else if (strncmp(optarg, "both", 4) == 0) {
					popt->unlimited_insertions = 3;
				}
				else {
					message(stderr, __FILE__, __LINE__,
							ERROR_MSG, INVALID_CMD_ARGUMENT, 
							"-i %s\n", optarg);

					exit(INVALID_CMD_ARGUMENT);
				}

				break;
			/* random seed */
			case 's':
				rand_seed = strtol(optarg, &ptr_optarg, 0);
				if (rand_seed == 0) rand_seed = time(NULL);
				srand((unsigned int) rand_seed);

				break;
			case 'h':
				print_usage(argv[0]);
				exit(0);
				break;
		}
	}

	int invalid_opts = 0;
	if (optind < argc) {
		popt->fastq_file = argv[optind];
	}
	else {
		invalid_opts |= 1;
	}

	/* validate specified options/arguments */
	if (popt->kmer_length < 4 || popt->kmer_length > 31) {
		message(stderr, __FILE__, __LINE__, ERROR_MSG, 
				INVALID_CMD_ARGUMENT, "kmer size should be within range [4, 32].\n");
		invalid_opts |= 1;
	}

	if (popt->max_edit_dist < 1  ||
			popt->max_edit_dist > popt->kmer_length) {
		message(stderr, __FILE__, __LINE__, ERROR_MSG, 
				INVALID_CMD_ARGUMENT, 
				"-d %s : should be within range [1, min(k, 15)].\n", optarg);
		invalid_opts |= 1;
	}

	if (popt->penalty_eta < 0) {
		message(stderr, __FILE__, __LINE__,
				ERROR_MSG, INVALID_CMD_ARGUMENT, 
				"-L %s : negative lambda value.\n", optarg);
		invalid_opts |= 1;
	}

	if (popt->penalty_gamma < 0) {
		message(stderr, __FILE__, __LINE__,
				ERROR_MSG, INVALID_CMD_ARGUMENT, 
				"-G %s : negative gamma value.\n", optarg);
		invalid_opts |= 1;
	}

	if (popt->tol < 0 || popt->tol > 1) {
		message(stderr, __FILE__, __LINE__,
				ERROR_MSG, INVALID_CMD_ARGUMENT, 
				"-R %s : should be within range (0, 1).\n", optarg);
		invalid_opts |= 1;
	}

	if (invalid_opts) {
		print_usage(argv[0]);
		exit(1);
	}

	return NO_ERROR;
}

void print_usage(const char *basename)
{
	fprintf(stderr, pmr_usage, basename);
}

void timer(int report)
{
	static struct timeval tv;
	if (report == 0)
		gettimeofday(&tv, NULL);
	else {
		struct timeval _ntv;
		gettimeofday(&_ntv, NULL);
		fprintf(stderr, 
				"[TIMER] Elapsed: %.3f\n",
				(0.0 + _ntv.tv_sec + (double) _ntv.tv_usec/1e6) - 
				(0.0 + tv.tv_sec + (double) tv.tv_usec/1e6));

		tv = _ntv;
	}
}

