#ifndef __PREMIER_TRANS_H__
#define __PREMIER_TRANS_H__

#include <stdint.h>

#include "kmer.h"
#include "premier.h"


extern const uint16_t trans_flag2packed_lookup_table[];
/* translate unique transition flag to corresponding index/offset in 
 * data->transition_p. */
extern const int trans_flag2index_lookup_table[];

extern const kbits_t trans_valid_trans_masks[];

#define trans_flag2packed(f) trans_flag2packed_lookup_table[(f)]
#define trans_flag2index(f) trans_flag2index_lookup_table[(f)]

/*
static inline kbits_t trans_admissible_transitions(kbits_t kid, 
		kbits_t obs_kid, kbits_t obs_n_flag, kbits_t k_trans_flag, 
		kbits_t obs_next_base, kbits_t hd, kbits_t dmax)
{
	kbits_t mut_flag = ((kid ^ obs_kid) | obs_n_flag) & 3UL;
	kbits_t sfx_hd = hd - ((mut_flag & 1UL) | (mut_flag >> 1));

	kbits_t adm_trans_mask = (((uint64_t) (sfx_hd >= dmax)) - 1UL) |
		(1UL << obs_next_base);

	return ((sfx_hd << 4) | (k_trans_flag & adm_trans_mask));
}
*/

static inline kbits_t trans_admissible_transitions(kbits_t kid,
	kbits_t obs_n_flag, kbits_t trans_flag, kbits_t obs_next_base, int kmer_len, 
	int unconstrained)
{
#if PMR_EMIT_OBSERVED_BASE
	/* check for valid/admissible transitions such that the last base of 
	the kmer/(k+1)-mer agrees with obs_next_base, unless obs_next_base is 'N',
	in which case any transition is allowed. */
	kbits_t klf = kid >> 62;
	kbits_t valid_trans_flag = trans_flag & trans_valid_trans_masks[obs_next_base];

	if (obs_next_base < 4 && klf == 0 && kbits_last_base(kid, kmer_len) == obs_next_base) {
		valid_trans_flag |= (1UL << 20);
	}

	return valid_trans_flag;
#else
	/* allow for substitution, deletion error, and insertion */
	kbits_t klf = kid >> 62;
	kbits_t valid_trans_flag = trans_flag & trans_valid_trans_masks[obs_next_base];

	if (obs_next_base < 4 && klf == 0) { // && kbits_last_base(kid, kmer_len) == obs_next_base) {
		valid_trans_flag |= (1UL << 20);
	}

	return unconstrained ? valid_trans_flag : (trans_flag & (1 << obs_next_base));
#endif
}

#endif
