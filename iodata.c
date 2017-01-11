#include <omp.h>

#include "iodata.h"
#include "read.h"
#include "kvec.h"

#define PROG_BAR_LEN 25
#define PROG_REPORT_INTERVAL 500

#define PMR_BATCH_GROUP_SIZE 5000

void io_draw_progress_bar(int i, int n, double elapsed, int finished)
{
	i = (i > n) ? n : i;

	if (finished) {
		fprintf(stderr, "\rReads proc'ed: %8d [", n);
		for (int j = 0; j < PROG_BAR_LEN; ++j) {
			fputc('.', stderr);
		}

		fprintf(stderr, "] 100.0%% Elapsed: %02d:%02d:%6.3f\n", 
				(int) elapsed / 3600, ((int) elapsed % 3600) / 60, 
				elapsed - ((int) elapsed / 60) * 60
			);
	}
	else {
		// ETA, in seconds
		int eta = elapsed * ((double) n / i) - elapsed;
		int progress = i * PROG_BAR_LEN / n;

		fprintf(stderr, "\rReads proc'ed: %8d [", i);

		for (int j = 0; j <= progress; ++j) {
			fputc('.', stderr);
		}
		for (int j = progress + 1; j < PROG_BAR_LEN; ++j) {
			fputc(' ', stderr);
		}

		fprintf(stderr, "] %5.1f%% ETA: %02d:%02d:%02d ", 
				(double) i / n * 100.0,
				eta < 0 ? 99 : eta / 3600, 
				eta < 0 ? 99 : (eta % 3600) / 60, 
				eta < 0 ? 99 : eta % 60);
	}
}

char *io_cat_output_filename(options *opt, const char *suffix, 
		char *ext)
{
	time_t t;
	struct tm *lts;
	char *fname;
	char time_string[32];
	char eta_string[32];

	fname = calloc(PMR_LINE_BUFFER_SIZE, sizeof(char));
	if (fname == NULL) return NULL;

	/* generate time info */
	time(&t);
	lts = localtime(&t);

	if (ext == NULL) {
		ext = "fastq";
	}

	if (opt->filename_append_timeinfo) {
		strftime(time_string, 31, "_%m%d_%H%M", lts);
	}
	else {
		time_string[0] = '\0';
	}

	if (!isnan(opt->penalty_eta)) {
		snprintf(eta_string, 31, "%0.f", opt->penalty_eta);
	}
	else {
		snprintf(eta_string, 31, "%s", "auto");
	}

	if (opt->filename_append_params) {
		snprintf(fname, PMR_LINE_BUFFER_SIZE - 1, 
				"%s%s_K%dD%d_E%sG%1.0eR%1.0e%s.%s", 
				opt->output_prefix, time_string, 
				opt->kmer_length, opt->max_edit_dist, 
				eta_string, opt->penalty_gamma, opt->tol,
				suffix, ext);
	}
	else {
		snprintf(fname, PMR_LINE_BUFFER_SIZE - 1, 
				"%s%s%s.%s", 
				opt->output_prefix, time_string, 
				suffix, ext);
	}

	return fname;
}

int io_iterate_reads_fastq(data *d, int nthreads, 
		int progress_bar,
		void (*callback)(data *d, read_t *read, void *fdata), 
		void *fdata)
{
	int fd;
	struct stat fs;
	struct timeval tv;
	int64_t fastq_size;
	size_t size_ident;
	int read_idx = 0, n = 15;
	int read_processed = 0;
	char *pnr = NULL;

	gettimeofday(&tv, NULL);

	if ((fd = open(d->opt->fastq_file, O_RDONLY)) == -1){
		return message(stderr, __FILE__, __LINE__, ERROR_MSG,
				FILE_NOT_FOUND, d->opt->fastq_file);
	}
	else {
		if (fstat(fd, &fs) != -1)
			fastq_size = (int64_t) fs.st_size;
		else 
			return message(stderr, __FILE__, __LINE__, ERROR_MSG, 
				FILE_NOT_FOUND, d->opt->fastq_file);
	}

	close(fd);
	
	/* create a temporary read */
	read_t *tmp_read = read_create();
	tmp_read->identifier = NULL; // buf_ident;

	int chunk = 0;
	int chunks_to_read = (int) ceil((double) fastq_size / PMR_READ_BUFFER_SIZE);

	// pad 16 bytes to avoid invalid read with SSE intrinsics (_mm_loadu_xx).
	char *_buffer = malloc(16 + (fastq_size < PMR_READ_BUFFER_SIZE ? fastq_size :
			PMR_READ_BUFFER_SIZE));
	char *buf_loc = _buffer;

	/* process in blocks */
	FILE *fp = fopen(d->opt->fastq_file, "r");
	while (fastq_size > 0) {
		size_t chunk_size = fastq_size + (buf_loc - _buffer) < PMR_READ_BUFFER_SIZE ? 
			fastq_size : PMR_READ_BUFFER_SIZE - (buf_loc - _buffer);

		size_t buf_size = fastq_size + (buf_loc - _buffer) < PMR_READ_BUFFER_SIZE ? 
			fastq_size + (buf_loc - _buffer) :
			PMR_READ_BUFFER_SIZE;

		kvec_t(read_t) read_batch;
		kv_init(read_batch);
		kv_resize(read_t, read_batch, buf_size / (128 << 2));

		char *buf_max = _buffer + buf_size - 1;

		char *p = _buffer;
		++chunk;

		if (!progress_bar) { 
			INFO("Reading in FASTQ data in chunk [%2d/%2d] (%zu / REM: %zu)...\n", chunk, 
					chunks_to_read, chunk_size, fastq_size);
		}

		fread(buf_loc, chunk_size, 1, fp);

		size_t partial_read_size = 0;
		char *pread = p;

		/* locate the last read in buffer */
		//size_t rem_buf_size = buf_size;
		while ((p - _buffer) < buf_size) {
			pread = p;
			buf_loc = _buffer;
			tmp_read->id = read_idx++;

			IO_LOCATE_NEWLINE;

			/* process header */
			char *p_ident = p + 1;
			for (; *p_ident != ' ' && p_ident < pnr; p_ident++);

			size_ident = (size_t) (p_ident - p - 1);
			tmp_read->identifier = p + 1;

			/*
			strncpy(buf_ident, p + 1, size_ident);
			buf_ident[size_ident] = '\0';
			*/

			/* read sequence */
			p = pnr + 1;
			tmp_read->sequence = p;

			IO_LOCATE_NEWLINE;

			int read_length = pnr - p;

			tmp_read->length = read_length;

			/* skip quality header line */
			p = pnr + 1;
			if (*p == '+') {
				IO_LOCATE_NEWLINE;
			}

			/* quality score */
			p = pnr + 1;
			tmp_read->qscore = p;
			IO_LOCATE_NEWLINE;

			/* invoke callback function */
			//callback(d, tmp_read, fdata);

			// overwrite buffer raw data.
			tmp_read->identifier[size_ident] = '\0';
			kv_push(read_t, read_batch, *tmp_read);

			p = pnr + 1;
		}

		size_t batch_size = kv_size(read_batch);
		size_t batch_groups = batch_size / PMR_BATCH_GROUP_SIZE;
		if (batch_groups * PMR_BATCH_GROUP_SIZE < batch_size) ++batch_groups;

		int tid;
		double wt_start, wt_end;
		// reads processed by all threads in current batch
		int batch_processed = 0;

		omp_set_num_threads(nthreads);
		#pragma omp parallel private(tid, wt_start, wt_end)
		{ //omp

		tid = omp_get_thread_num();
		wt_start = omp_get_wtime();

		for (int b = 0; b < batch_groups; ++b) {
			int imin = b * PMR_BATCH_GROUP_SIZE;
			int imax = (b+1) * PMR_BATCH_GROUP_SIZE > batch_size ? batch_size :
				(b+1) * PMR_BATCH_GROUP_SIZE;

			#pragma omp for schedule(dynamic) reduction(+:batch_processed)
			for (int i = imin; i < imax; ++i) {
				callback(d, &kv_A(read_batch, i), fdata);
				++batch_processed;
			}

			// try resize, only performed at the master thread.
			#pragma omp master 
			{
				if (progress_bar && batch_processed % PROG_REPORT_INTERVAL == 0) {
					wt_end = omp_get_wtime();
					io_draw_progress_bar(read_processed + batch_processed, 
							d->n_reads, wt_end - wt_start, 0);
				}

				read_processed += batch_processed;
				batch_processed = 0;

				dict_try_resize(d->kmer_dict);
			}

			#pragma omp barrier
		}

		wt_end = omp_get_wtime();
		#pragma omp master
		{
			if (!progress_bar) {
				INFO("Threads finished iterating, time elapsed: %.3f\n", 
						wt_end - wt_start);
			}
		}

		} //end of omp region

		// move the partial read to the front of the buffer
		if (partial_read_size > 0) {
			memcpy(_buffer, pread, partial_read_size);
			buf_loc = _buffer + partial_read_size;
		}

		kv_destroy(read_batch);
		//read_processed += batch_size;

		fastq_size -= chunk_size;
	}


	struct timeval _ntv;
	gettimeofday(&_ntv, NULL);
	double elapsed = (0.0 + _ntv.tv_sec + (double) _ntv.tv_usec/1e6) - 
				(0.0 + tv.tv_sec + (double) tv.tv_usec/1e6);

	if (progress_bar) {
		io_draw_progress_bar(0, d->n_reads, elapsed, 1);
	}
	else {
		INFO("All (%9d) reads processed in %.3f seconds.\n",
				read_idx, elapsed);
	}

	d->n_reads = read_idx;

	read_destroy(tmp_read);

	fclose(fp);
	free(_buffer);
}

#if 0
void io_iterate_reads_fastq(data *d, 
		void (*callback)(data *d, read_t *read, void *fdata), void *fdata,
		int report_progress)
{
	static char buf_ident[PMR_LINE_BUFFER_SIZE];

	// progress indicator
	char prog_ind[PROG_BAR_LEN + 1];
	memset(prog_ind, ' ', PROG_BAR_LEN);
	prog_ind[PROG_BAR_LEN] = '\0';

	struct timeval tv, tv_s;
	int read_idx = 0, n = 15;
	char *p = d->fastq_mmap;

	gettimeofday(&tv, NULL);
	tv_s = tv;

	/* create a temporary read */
	read_t *tmp_read = read_create();
	tmp_read->identifier = buf_ident;

	int prev_prog = -1;
	/* fastq reader assuming the data is arranged in 4-line format,
	 * i.e., no sequence spans over multiple lines. */
	while ((*p != '\0') && (strchr(p, ' ') != 0)) {
		char *p_ident = strchr(p, ' ');
		size_t size_ident = (size_t) (p_ident - p);

		strncpy(buf_ident, p, size_ident);
		buf_ident[size_ident] = '\0';
		
		tmp_read->id = read_idx++;
		/* skip the title line */
		SKIP_TO_NEWLINE(p);
		tmp_read->sequence = p;

		/* skip to quality score title line, compute read length */
		SKIP_TO_NEWLINE(p);
		int read_length = p - tmp_read->sequence - 1;
		tmp_read->length = read_length;

		/* skip quality score title line */
		SKIP_TO_NEWLINE(p);

		tmp_read->qscore = p;
		SKIP_TO_NEWLINE(p);

		/* invoke callback function */
		callback(d, tmp_read, fdata);

		if (report_progress) {
			// show progress bar
			if (read_idx % PROG_REPORT_INTERVAL == 0) {
				gettimeofday(&tv, NULL);

				double sec_elapsed = (0.0 + tv.tv_sec + (double) tv.tv_usec/1e6) - 
						(0.0 + tv_s.tv_sec + (double) tv_s.tv_usec/1e6);
				// ETA, in seconds
				int eta = sec_elapsed * ((double) d->n_reads / read_idx) - sec_elapsed;

				int progress = read_idx * PROG_BAR_LEN / d->n_reads;

				if (progress > prev_prog) {
					prog_ind[progress] = '.';
					prev_prog = progress;
				}

				fprintf(stderr, "\rReads processed: %8d [%s] %5.1f%% ETA: %02d:%02d:%02d ", 
						read_idx,
						prog_ind,
						(double) read_idx / d->n_reads * 100.0,
						eta / 3600, 
						(eta % 3600) / 60, eta % 60
						);
			}
		}
		else {
			if (read_idx >> n) {
				gettimeofday(&tv, NULL);

				INFO("Processed %9d reads, time elapsed: %.3f\n", 
						read_idx,
						(0.0 + tv.tv_sec + (double) tv.tv_usec/1e6) - 
						(0.0 + tv_s.tv_sec + (double) tv_s.tv_usec/1e6));

				n++;
			}
		}

	}

	gettimeofday(&tv, NULL);
	double total_elapsed = (0.0 + tv.tv_sec + (double) tv.tv_usec/1e6) - 
				(0.0 + tv_s.tv_sec + (double) tv_s.tv_usec/1e6);

	if (report_progress) {
		fprintf(stderr, "\rReads processed: %8d [%s] %5.1f%% Elapsed: %02d:%02d:%06.3f\n", 
				d->n_reads,
				prog_ind,
				100.0,
				(int) total_elapsed / 3600, ((int) total_elapsed % 3600) / 60, 
				total_elapsed - ((int) total_elapsed / 60) * 60
			   );
	}
	else {
		INFO("All (%9d) reads processed in %.3f seconds.\n", 
				read_idx, total_elapsed);
	}

	d->n_reads = read_idx;

	read_destroy(tmp_read);
}
#endif
