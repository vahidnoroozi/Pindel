#include "trans.h"

const uint16_t trans_flag2packed_lookup_table[16] = {  
	0, 1, (1 << 3) | 1, (1 << 5) | 2, (2 << 3) | 1,
	(2 << 5) | 2, (2 << 5) | (1 << 3) | 2,
	(2 << 7) | (1 << 5) | 3, (3 << 3) | 1,
	(3 << 5) | 2, (3 << 5) | (1 << 3) | 2,
	(3 << 7) | (1 << 5) | 3, (3 << 5) | (2 << 3) | 2,
	(3 << 7) | (2 << 5) | 3,
	(3 << 7) | (2 << 5) | (1 << 3) | 3,
	(3 << 9) | (2 << 7) | (1 << 5) | 4
};

const int trans_flag2index_lookup_table[16] = {
	0, 1, 2, 0, 
	3, 0, 0, 0, 
	4, 0, 0, 0, 
	0, 0, 0, 0
};

#if PMR_EMIT_OBSERVED_BASE
const kbits_t trans_valid_trans_masks[5] = {
	1UL | (1UL << 4) | (1UL << 5)  | (1UL << 6) | (1UL << 7),
	2UL | (1UL << 8) | (1UL << 9)  | (1UL << 10) | (1UL << 11),
	4UL | (1UL << 12) | (1UL << 13) | (1UL << 14) | (1UL << 15),
	8UL | (1UL << 16) | (1UL << 17) | (1UL << 18) | (1UL << 19),
	UINT32_MAX,
};
#else
const kbits_t trans_valid_trans_masks[5] = {
	15UL | (1UL << 4) | (1UL << 5)  | (1UL << 6) | (1UL << 7),
	15UL | (1UL << 8) | (1UL << 9)  | (1UL << 10) | (1UL << 11),
	15UL | (1UL << 12) | (1UL << 13) | (1UL << 14) | (1UL << 15),
	15UL | (1UL << 16) | (1UL << 17) | (1UL << 18) | (1UL << 19),
	UINT32_MAX,
};
#endif
