#ifndef __KMER_H_
#define __KMER_H_

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

#ifdef __SSE3__
#include <pmmintrin.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __POPCNT__
#include <nmmintrin.h>
#endif

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "premier.h"

#define KBITS_KLEN_REMOVAL_FLAG ~(3UL << 62)
#define KBITS_HI_MASK 0xAAAAAAAAAAAAAAAAUL
#define KBITS_LO_MASK 0x5555555555555555UL
#define KBITS_MASK_MAX (__UINT64_C(18446744073709551615)) 
#define KMER_ID_FMT "%lu"

/* get k-1 suffix/prefix of a kmer bit field */
#define kbits_suffix(b) ((b) >> 2)
#define kbits_prefix(b, k) ((b) & ((1UL << ((k) << 1)) - 1))

/* (k+1)th base */
#define kbits_first_base(b) ((b) & 3)
#define kbits_kp1_base(b, k) ((b) << (k << 1))

/* convert a nucleotide literal (char type) to 2 bit representation */
#define base_to_bits(b) (((b) >> 1) & 3)
#define base_to_rc_bits(b) (((((b) >> 1) & 3) + 2) & 3)
#define complementary_base(b) (((b) + 2) & 3)

#define complementary_trans_flag(f) (((f) >> 2) | (((f) & 3) << 2))

/* cast a kbits_t type into a void pointer (used as key for a dictentry) */
#define kbits_cast_to_ptr(k) ((void *) ((uintptr_t) (k) + 1))

#define kmer_hamming_dist(p, q) kmer_n_hamming_dist(p, q, 0UL)

int kmer_incr_count_by_id(data *d, kbits_t kid, kbits_t next_base,
	kbits_t next_2bases, kmer_t **kaddr, uint32_t label);
void kmer_count_trans_by_id(data *d, kbits_t kid, kbits_t next_base);

static inline kbits_t kbits_raw_suffix(kbits_t k)
{
	register kbits_t kmer_len = k >> 62;
	return (k & KBITS_KLEN_REMOVAL_FLAG) >> (2 + (kmer_len << 1));
}

static inline kbits_t kbits_trans_suffix(kbits_t k)
{
	register kbits_t kmer_len = k >> 62;
	register kbits_t raw_suffix = (k & KBITS_KLEN_REMOVAL_FLAG) >> (2 + (kmer_len << 1));
	return raw_suffix | (kmer_len << 62);
}

static inline uint64_t kbits_suffix_edist(kbits_t k, uint64_t ed)
{
	register kbits_t kmer_len = k >> 62;
	return (ed >> (1+kmer_len));
}

static inline kbits_t kbits_discarded_prefix(kbits_t k)
{
	register kbits_t kmer_len = k >> 62;
	return k & ((1 << (2 + (kmer_len << 1))) - 1);
}

static inline kbits_t kbits_last_base(kbits_t k, int kmer_len)
{
	kmer_len += (k >> 62);
	return (k >> ((kmer_len - 1) << 1)) & 3;
}

static inline kmer_t *kmer_create(data *d, kbits_t kid)
{
	kmer_t *kmer = (kmer_t *) mempool_nalloc(d->mp, sizeof(*kmer), 16);

	kmer->id = kid;
	kmer->trans_flag = 0;
	kmer->count = 1.0;
	kmer->exp_count = 0.0;
	//kmer->unique_trans = 0;
	//kmer->emit_norm_flag = 0;
	//kmer->gamma = 0.0;

	return kmer;
}	/* kmer_create */

static inline void kmer_destroy(kmer_t *kmer)
{
	free(kmer);
}	/* kmer_destroy */

/* Convert a sequence literal to its bit representation */
static inline kbits_t kmer_seq_to_id(const char *s, int size, kbits_t len_flag)
{
	register int size_m128 = size < 16 ? size : 16;
	kbits_t kid = 0;

#ifdef __SSE2__
	char tmp[16] __attribute__ ((aligned (16)));

	__m128i base_mask = _mm_set1_epi8(3);
	__m128i seq_epi8 = _mm_loadu_si128((__m128i *) s);

	seq_epi8 = _mm_srli_epi64(seq_epi8, 1);
	seq_epi8 = _mm_and_si128(seq_epi8, base_mask);

	_mm_store_si128((__m128i *) tmp, seq_epi8);

	for (int i = 0; i < size_m128; i++){
		kid |= ((kbits_t) tmp[i] << (i << 1));
	}

	if (size > 16) {
		size_m128 = size - 16;
		seq_epi8 = _mm_and_si128(base_mask, 
				_mm_srli_epi64(_mm_loadu_si128(((__m128i *) s) + 1), 1));
		_mm_store_si128((__m128i *) tmp, seq_epi8);
		for (int i = 0; i < size_m128; i++) {
			kid |= ((kbits_t) tmp[i] << ((i+16) << 1));
		}
	}

#else

#endif
	return kid | (len_flag << 62);
}

static inline void kmer_id_to_seq(kbits_t k, int l)
{
	static char BASES[4] = {'A', 'C', 'T', 'G'};
	l += (k >> 62);
	for (int i = 0; i < l; i++, k >>= 2) {
		putchar(BASES[k & 3]);
	}
}

static inline void kmer_print_id(kbits_t k, int l, FILE *fp)
{
	static char BASES[4] = {'A', 'C', 'T', 'G'};
	l += (k >> 62);
	for (int i = 0; i < l; i++, k >>= 2) {
		fprintf(fp, "%c", BASES[k & 3]);
	}
}

static inline void kmer_sprint_id(kbits_t k, int l, char *kstr)
{
	static char BASES[4] = {'A', 'C', 'T', 'G'};
	l += (k >> 62);
	for (int i = 0; i < l; i++, k >>= 2) {
		kstr[i] = BASES[k & 3];
	}
}

static inline kbits_t kmer_n_base_flag(const char *s, int size)
{
	kbits_t nflag = 0;
	for (int i = 0; i < size; i++) {
		register kbits_t nf = (s[i] >> 3) & 1;
		nflag |= (((nf << 1) | nf) << (i << 1));
	}

	return nflag;
}

/* compute hamming distance between two kmers p and q.
 * if user's CPU supports POPCNT instruction, this function can be fully
 * optimized. */
static inline int kmer_n_hamming_dist(kbits_t p, kbits_t q, kbits_t n_flag)
{
	kbits_t xor_pq, xor_lo, xor_hi;

	xor_pq = (p ^ q) | n_flag;
	xor_lo = xor_pq & KBITS_LO_MASK;
	xor_hi = xor_pq & KBITS_HI_MASK;
#ifdef __POPCNT__
	/* POPCNT: count # of 1's in an unsigned integer */
	register int hd = _mm_popcnt_u64(xor_hi);

	xor_lo &= ~(xor_hi >> 1);
	hd += _mm_popcnt_u64(xor_lo);
#else
	register int hd = 0;
	/* TODO: fallback code */
#endif
	return hd;
}

static inline kbits_t kmer_effective_base(uint8_t c) {
	/* if c is not 'N', not_N_mask is UINT_MAX, 
	 * otherwise, 0.  */
	kbits_t not_N_mask = (((kbits_t) c >> 3) & 1) - 1;
	kbits_t base = (c >> 1) & 3;

	return (not_N_mask & base) | (~not_N_mask & 4);
}

static inline kbits_t kmer_reverse_complement(kbits_t k, int size)
{
	//static kbits_t rev_map[4] = {2, 3, 0, 1};
	uint64_t ktype = k >> 62;

	register kbits_t rc_k = 0;
	for (int i = size - 1; i >= 0; i--) {
		register int isll1 = i << 1;
		rc_k |= (complementary_base(k & 3) << isll1);
		//rc_k |= (rev_map[k & 3] << isll1);

		k >>= 2;
	}

	rc_k |= (ktype << 62);

	return rc_k;
}

#endif
