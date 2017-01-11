#include "read.h"
#include "kmer.h"

/** Decompose a read into kmers, record kmer counts, construct
 * kmer neighborhood, and add kmer into the k-spectrum unless
 * the kmer has undetermined base(s).
 */
void read_parse(data *d, read_t *tmp_read, void *fdata)
{
	int kmer_len = d->opt->kmer_length;
	int read_len = tmp_read->length;
	int tmax = read_len - kmer_len;
	kbits_t kp1mer_reset_flag = ~(3UL << 62);
	int k_shift = (kmer_len - 1) << 1;
	int kp1_shift = (kmer_len) << 1;

	if (read_len <= kmer_len) return;

	char *rseq = tmp_read->sequence;
	
	kmer_t *pkmer;

	if (tmp_read->length > d->max_read_length)
		d->max_read_length = tmp_read->length;
		
	kbits_t kmer_id = kmer_seq_to_id(rseq, kmer_len, 0);
	kbits_t kmer_rc_id = kmer_reverse_complement(kmer_id, kmer_len);

	kbits_t kp1mer_id = kmer_seq_to_id(rseq, kmer_len+1, 1);
	kbits_t kp1_rc_id = kmer_reverse_complement(kp1mer_id, kmer_len+1);
	
	kbits_t k2k_base = kmer_effective_base(rseq[kmer_len]);
	kbits_t k2kp1_base = (read_len - kmer_len < 2) ? 0xFF :
		k2k_base + kmer_effective_base(rseq[kmer_len+1]) * 4;
	//kbits_t kp12k_base = kmer_effective_base(rseq[kmer_len+1]);
	//kbits_t kp12kp1_base = kmer_effective_base(rseq[kmer_len+2]) * 4 + kp12k_base;

	kmer_incr_count_by_id(d, kmer_id, k2k_base, k2kp1_base, &pkmer, 0);
	kmer_incr_count_by_id(d, kmer_rc_id, 0xFF, 0xFF, &pkmer, 1);

	//kmer_incr_count_by_id(d, kp1mer_id, kp12k_base, kp12kp1_base, &pkmer);
	kmer_incr_count_by_id(d, kp1mer_id, 0xFF, 0xFF, &pkmer, 0);
	kmer_incr_count_by_id(d, kp1_rc_id, 0xFF, 0xFF, &pkmer, 0);

	for (int t = 1; t < tmax; t++) {
		//kbits_t kth_n_flag = (rseq[t - 1 + kmer_len] >> 3) & 1;
		//kbits_t kp1th_n_flag = (rseq[t - 1 + kmer_len] >> 3) & 1;
		kbits_t rc_next_base = base_to_rc_bits(rseq[t - 1]);
		kbits_t rc_kp1_base = t == 1 ? 0xFF :
			base_to_rc_bits(rseq[t - 2]) * 4 + rc_next_base;

		kbits_t kth_base = base_to_bits(rseq[t - 1 + kmer_len]);
		kbits_t kp1th_base = base_to_bits(rseq[t + kmer_len]);

		//obs_next_base = kmer_effective_base(rseq[t + kmer_len]);

		//register kbits_t _nf = ((kbits_t) kth_n_flag << 1) | kth_n_flag;
		//kmer_n_flag = kbits_suffix(kmer_n_flag) | (_nf << shift);

		/* update kmer id and k+1-mer id */
		kmer_id = kbits_suffix(kmer_id) | ((kbits_t) kth_base << k_shift);
		kmer_rc_id = kmer_reverse_complement(kmer_id, kmer_len);

		/* FIXME: add a dedicated function/macro for this mess */
		kp1mer_id = kbits_suffix((kp1mer_id & kp1mer_reset_flag)) | ((kbits_t) kp1th_base << kp1_shift) | (1UL << 62);
		kp1_rc_id = kmer_reverse_complement(kp1mer_id, kmer_len+1);
		
		k2k_base = base_to_bits(rseq[t + kmer_len]);
		if (t < tmax - 1) {
				k2kp1_base = base_to_bits(rseq[t + 1 + kmer_len]) * 4 + k2k_base;
				//kp12k_base = base_to_bits(rseq[t + 1 + kmer_len]);
		}
		else {
			k2kp1_base = 0xFF;
			//kp12k_base = 0xFF;
		}

		//if (t < tmax - 2) {
		//	kp12kp1_base = base_to_bits(rseq[t + 2 + kmer_len]) * 4 + kp12k_base;
		//}
		//else {
		//	kp12kp1_base = 0xFF;
		//}
				
		/* skip any kmer with at least one 'N' base */
		//if (kmer_n_flag) continue;
		
		kmer_incr_count_by_id(d, kmer_id, k2k_base, k2kp1_base, &pkmer, 0);
		kmer_incr_count_by_id(d, kmer_rc_id, rc_next_base, 
				rc_kp1_base, &pkmer, 1);

		//kmer_incr_count_by_id(d, kp1mer_id, kp12k_base, kp12kp1_base, &pkmer);
		kmer_incr_count_by_id(d, kp1mer_id, 0xFF, 0xFF, &pkmer, 0);
		kmer_incr_count_by_id(d, kp1_rc_id, 0xFF, 0xFF, &pkmer, 0);
	}

	/* last kmer */
	//kbits_t kth_n_flag = (rseq[tmax + kmer_len - 1] >> 3) & 1;
	//register kbits_t _nf = ((kbits_t) kth_n_flag << 1) | kth_n_flag;

	//kmer_n_flag = kbits_suffix(kmer_n_flag) | (_nf << shift);

	//if (kmer_n_flag == 0) {
		kbits_t kth_base = (rseq[tmax + kmer_len - 1] >> 1) & 3;
		kmer_id = kbits_suffix(kmer_id) | ((kbits_t) kth_base << k_shift);
		kmer_rc_id = kmer_reverse_complement(kmer_id, kmer_len);

		kbits_t rc_next_base = base_to_rc_bits(rseq[tmax - 1]);
		kbits_t rc_kp1_base = tmax < 2 ? 0xFF :
			base_to_rc_bits(rseq[tmax - 2]) * 4 + rc_next_base;

		kmer_incr_count_by_id(d, kmer_id, 0xFF, 0xFF, &pkmer, 0);
		kmer_incr_count_by_id(d, kmer_rc_id, rc_next_base, rc_kp1_base, &pkmer, 1);
	//}
}	/* read_parse */

void read_parse_trans_count(data *d, read_t *tmp_read, void *fdata)
{
	int kmer_len = d->opt->kmer_length;
	int read_len = tmp_read->length;
	int tmax = read_len - kmer_len;
	kbits_t kp1mer_reset_flag = ~(3UL << 62);
	int k_shift = (kmer_len - 1) << 1;
	int kp1_shift = (kmer_len) << 1;
	int qoff = d->opt->qual_score_offset;
	int qmax = d->opt->max_qual_score;

	if (read_len <= kmer_len) return;

	char *rseq = tmp_read->sequence;
	
	kmer_t *pkmer;

	for (int p = 0; p < read_len; ++p) {
		int Q = tmp_read->qscore[p];

		if (qmax < Q) {
			qmax = Q;
		}
		if (qoff > Q) {
			qoff = Q; 
		}
	}

	d->opt->qual_score_offset = qoff;
	d->opt->max_qual_score = qmax;
		
	kbits_t kmer_id = kmer_seq_to_id(rseq, kmer_len, 0);
	//kbits_t kp1mer_id = kmer_seq_to_id(rseq, kmer_len+1, 1);
	
	/*
	kbits_t kmer_n_flag = kmer_n_base_flag(rseq, kmer_len);
	*/

	kbits_t k2k_base = kmer_effective_base(rseq[kmer_len]);
	//kbits_t kp12k_base = kmer_effective_base(rseq[kmer_len+1]);

	//if (kmer_n_flag == 0) {
		kmer_count_trans_by_id(d, kmer_id, k2k_base);
		//kmer_count_trans_by_id(d, kp1mer_id, kp12k_base);
	//}

	for (int t = 1; t < tmax; t++) {
		//kbits_t kth_n_flag = (rseq[t - 1 + kmer_len] >> 3) & 1;
		//kbits_t kp1th_n_flag = (rseq[t - 1 + kmer_len] >> 3) & 1;

		kbits_t kth_base = base_to_bits(rseq[t - 1 + kmer_len]);
		//kbits_t kp1th_base = base_to_bits(rseq[t + kmer_len]);

		//obs_next_base = kmer_effective_base(rseq[t + kmer_len]);

		//register kbits_t _nf = ((kbits_t) kth_n_flag << 1) | kth_n_flag;
		//kmer_n_flag = kbits_suffix(kmer_n_flag) | (_nf << shift);

		/* update kmer id and k+1-mer id */
		kmer_id = kbits_suffix(kmer_id) | ((kbits_t) kth_base << k_shift);
		/* FIXME: add a dedicated function/macro for this mess */
		//kp1mer_id = kbits_suffix((kp1mer_id & kp1mer_reset_flag)) | ((kbits_t) kp1th_base << kp1_shift) | (1UL << 62);
		
		k2k_base = base_to_bits(rseq[t + kmer_len]);
		/*
		if (t < tmax - 1) {
			kp12k_base = base_to_bits(rseq[t + 1 + kmer_len]);
			kmer_count_trans_by_id(d, kp1mer_id, kp12k_base);
		}
		*/

		/* skip any kmer with at least one 'N' base */
		//if (kmer_n_flag) continue;
		kmer_count_trans_by_id(d, kmer_id, k2k_base);

		//if (kp1th_n_flag) continue;
	}
}	/* read_parse */
