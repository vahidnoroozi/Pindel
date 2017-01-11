#ifndef __PREMIER_IODATA_H__
#define __PREMIER_IODATA_H__

#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <unistd.h>
#include <ctype.h>

#include "premier.h"

#define SKIP_TO_NEWLINE(p)      \
	for(; *(p) != '\n'; (p)++); \
	(p)++;                      \

#define IO_LOCATE_NEWLINE												\
	for(pnr = p; *pnr != '\n' && pnr < buf_max; pnr++);					\
	if (*(pnr) != '\n') {												\
		partial_read_size = (pnr - pread + 1);							\
		break;															\
	}																	\

/*
int setup_worker_entries(data *d);
*/

char *io_cat_output_filename(options *opt, const char *suffix, 
		char *ext);

int io_iterate_reads_fastq(data *d, int nthreads, int progress_bar,
		void (*callback)(data *d, read_t *read, void *fdata), void *fdata);
		
#endif
