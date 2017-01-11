#include "mstep.h"
#include "trans.h"

#include <signal.h>

#define EPSILON 1e-6

static double mstep_estimate_trans_p(data *d, kmer_t *pk, double *jnt_exp_trans,
		double thresh, mstep_data_t *md);

static void mstep_reciprocal_exp_trans(data *d, kmer_t *pk, double *recip_exp_trans)
{

	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) << 1;
	kbits_t k_trans_packed = trans_flag2packed(pk->trans_flag & 15);
	kbits_t n_trans = k_trans_packed & 7;
	int i = 0;

	kbits_t kid = pk->id; 

	for (k_trans_packed >>= 3; i < n_trans; i++, k_trans_packed >>= 2) {
		kbits_t dns_base = k_trans_packed & 3;
		kbits_t rc_dns_base = complementary_base(kid & 3);

		kbits_t dns_kid = kbits_suffix(kid) | (dns_base << shift);
		kbits_t rc_dns_kid = kmer_reverse_complement(dns_kid, kmer_len);
		kmer_t *rc_pk = dict_find(d->kmer_dict, kbits_cast_to_ptr(rc_dns_kid))->value;

		if ((rc_pk->trans_flag & (1 << rc_dns_base)) == 0) {
			/* no such transition on reverse complement kmer */
			recip_exp_trans[i] = 0.0;
		}
		else {
			int rc_nbidx = _mm_popcnt_u64(rc_pk->trans_flag &  
					((1 << rc_dns_base) - 1));

			recip_exp_trans[i] = rc_pk->exp_trans_p[rc_nbidx];
		}
	}
}

static double mstep_compute_reverse_trans(data *d, kmer_t *rev_pk, 
		int kmer_len, tflag_t *tf)
{
	int shift = (kmer_len - 1) * 2;
	kbits_t rev_kid = rev_pk->id; 
	kbits_t k_trans_packed = trans_flag2packed(rev_pk->trans_flag & 15);
	kbits_t n_trans = k_trans_packed & 7;
	int i = 0;

	double sum_tp = 0.0;
	for (k_trans_packed >>= 3; i < n_trans; i++, k_trans_packed >>= 2) {
		kbits_t dns_base = k_trans_packed & 3;
		kbits_t fwd_dns_base = complementary_base(rev_kid & 3);

		kbits_t rc_fwd_kid = kbits_suffix(rev_kid) | (dns_base << shift);
		kbits_t fwd_kid = kmer_reverse_complement(rc_fwd_kid, kmer_len);

		kmer_t *fwd_pk = dict_find(d->kmer_dict, 
				kbits_cast_to_ptr(fwd_kid))->value;
		double mu_dns = fwd_pk->init_p; 

		if ((fwd_pk->trans_flag & (1 << fwd_dns_base)) == 0) {
			rev_pk->transition_p[i] = 0.0;
		}
		else {
			double rtp = PRODUCT(mu_dns, -rev_pk->init_p); 
			if (fwd_pk->n_trans > 1) {
				int fwd_nbidx = _mm_popcnt_u32(fwd_pk->trans_flag & 15 &
						((1 << fwd_dns_base) - 1));
				rtp += fwd_pk->transition_p[fwd_nbidx];
			}
			rev_pk->transition_p[i] = EXP(rtp);
			sum_tp += rev_pk->transition_p[i];
		}
	}

	int transidx = 0;
	k_trans_packed = trans_flag2packed(rev_pk->trans_flag & 15) >> 3;
	for (int j = 0; j < n_trans; ++j, k_trans_packed >>= 2) {
		kbits_t base = k_trans_packed & 3;

		if (rev_pk->transition_p[j] == 0.0) {
			rev_pk->trans_flag &= ~(1U << base);
			// (k+1)-mers flags
			for (int b2 = 4; b2 <= 16; b2 += 4) {
				rev_pk->trans_flag &= ~(1U << (b2 + base));
			}
		}
		else {
			rev_pk->transition_p[transidx++] = LOG(
					rev_pk->transition_p[j] / sum_tp);
		}
	}

	rev_pk->n_trans = _mm_popcnt_u32(
			rev_pk->trans_flag & ((1 << 20) - 1));
	rev_pk->n_ktrans = _mm_popcnt_u32(rev_pk->trans_flag & 15);

	return sum_tp;
}

void mstep_transition(data *d)
{
	/*
	static char *BASES[] = {"A", "C", "T", "G",
		"AA", "CA", "TA", "GA",
		"AC", "CC", "TC", "GC",
		"AT", "CT", "TT", "GT",
		"AG", "CG", "TG", "GG",
		"self"};
	*/

	double jnt_exp_trans[4];

	int kmer_len = d->opt->kmer_length;
	int shift = (kmer_len - 1) << 1;
	int n_zero_exp_cnt_kmer = 0;
	d->total_kmer_exp_counts = 0.0;

	if (d->opt->disable_penalized_em) {
		dictentry *de = d->kmer_dict->table;
		for (int i = d->kmer_dict->used; i > 0; de++) {
			if (de->value == NULL) continue;

			i--;
			kmer_t *pk = de->value;

			d->total_kmer_exp_counts += pk->exp_count;

			if (pk->exp_count == 0.0) {
				//pk->trans_flag = 0;
				for (int j = 0; j < pk->n_ktrans; j++) {
					pk->transition_p[j] = NAN;
				}

				n_zero_exp_cnt_kmer++;
			}
			else {
				int transidx = 0;
				double sum_trans = 0.0;
				for (int j = 0; j < pk->n_ktrans; j++) {
					sum_trans += pk->exp_trans_p[j];
				}

				for (int j = 0; j < pk->n_ktrans; j++) {
					pk->transition_p[j] = LOG(
							pk->exp_trans_p[j] / sum_trans);
				}
			}
		}
	}
	else {
		d->total_kmer_exp_counts = 0.0;
		double thresh = d->opt->penalty_eta / log1p(1/d->opt->penalty_gamma);
		double loglik_penalty = 0.0;

		mstep_data_t md = {
			.kappa = d->opt->penalty_gamma,
			/* use eta to not confuse with the lambda in 
			 * Lagrange multiplier */
			.eta = d->opt->penalty_eta, 
			.xi = NULL,
			.signs = 0,
			.n_xi = 0 
		};
		
		int n_rev_kmers = 0, n_fwd_kmers = 0;
		dictentry *de = d->kmer_dict->table;
		for (int i = d->kmer_dict->used; i > 0; de++) {
			if (de->value == NULL) continue;

			i--;
			kmer_t *pk = de->value;

			if (pk->id >> 62) continue;

			//pk->exp_count = 0.0;
			for (int b = 0; b < pk->n_ktrans; ++b) {
				pk->exp_count += pk->exp_trans_p[b];
			}
			d->total_kmer_exp_counts += pk->exp_count;

			if (pk->label == 1) continue;
			
			mstep_reciprocal_exp_trans(d, pk, jnt_exp_trans);
			for (int b = 0; b < pk->n_ktrans; ++b) {
				jnt_exp_trans[b] += pk->exp_trans_p[b];
			}

			loglik_penalty += mstep_estimate_trans_p(d, pk, jnt_exp_trans, thresh, &md);

			pk->dirty = 0;
			++n_fwd_kmers;
		}

		de = d->kmer_dict->table;
		for (int i = d->kmer_dict->used; i > 0; de++) {
			if (de->value == NULL) continue;

			i--;
			kmer_t *pk = de->value;

			if (pk->id >> 62 || pk->label == 0) continue;

			tflag_t _tf = {0, 0};

			double sum_tp = mstep_compute_reverse_trans(d, pk, kmer_len, &_tf); 

			pk->dirty = 0;
			++n_rev_kmers;
		}

		INFO("Total kmers: %d, forward stranded kmers (label=1): %d\n",
				n_fwd_kmers + n_rev_kmers, n_fwd_kmers);

		INFO("Remove dangling k-(k+1) transition flags...\n");

		de = d->kmer_dict->table;
		for (int i = d->kmer_dict->used; i > 0; de++) {
			if (de->value == NULL) continue;

			i--;
			kmer_t *pk = de->value;
			if ((pk->id >> 62)) continue;

			for (kbits_t b = 4; b < 20; ++b) {
				if ((pk->trans_flag & (1 << b)) == 0) continue;

				kbits_t base1 = (b-4) & 3;
				kbits_t base2 = (b-4) >> 2;

				kbits_t kid2 = kbits_raw_suffix(pk->id) | (base1 << shift);
				kmer_t *pk2 = dict_find(d->kmer_dict, kbits_cast_to_ptr(kid2))->value;

				if ((pk2->trans_flag & (1U << base2)) == 0) {
					// dangling k-(k+1) transition, needs to be removed
					pk->trans_flag &= ~(1U << b);
					--pk->n_trans;
				}
			}
		}

		// last iteration: for (k+1)-mers: update the transition flags
		de = d->kmer_dict->table;
		for (int i = d->kmer_dict->used; i > 0; de++) {
			if (de->value == NULL) continue;

			i--;
			kmer_t *pkp1 = de->value;

			if ((pkp1->id >> 62) == 0) continue;

			kbits_t _kid = (pkp1->id & KBITS_KLEN_REMOVAL_FLAG) >> 2;
			kmer_t *pk = dict_find(d->kmer_dict, kbits_cast_to_ptr(_kid))->value;
			
			pkp1->trans_flag = pk->trans_flag & ((1U << 20) - 1);
			pkp1->n_trans = pk->n_trans;
			pkp1->n_ktrans = pk->n_ktrans;
		}

		d->loglik -= loglik_penalty;
	}
}

static inline double mstep_estimate_trans_p(data *d, kmer_t *pk, 
		double *jnt_exp_trans, double thresh, mstep_data_t *md)
{
	kbits_t n_trans = pk->n_ktrans; 
	double _k_penalty = 0.0, max_xi = 0.0;
	double sum_xi = 0.0;
	int argmax_xi = -1;

	int32_t pos_signs = (1 << n_trans) - 1;
	int32_t signs = pos_signs;
	int32_t alt_signs = signs;

	int nonzero_trans = 0;
	int transidx = 0;
	double l0_pen = 0.0;

	for (int j = 0; j < n_trans; j++) {
		_k_penalty += thresh * 
			log1p(EXP(pk->transition_p[j]) / md->kappa);

		if (jnt_exp_trans[j] > max_xi) {
			max_xi = jnt_exp_trans[j];
			argmax_xi = j;
		}

		sum_xi += jnt_exp_trans[j];

		if (jnt_exp_trans[j] > 0.0) {
			l0_pen += (jnt_exp_trans[j] > 0) * thresh;
			nonzero_trans++;
		}
		/*
		else if (pk->exp_trans_p[transidx] == 0.0) {
			pk->trans_flag &= ~(1U << j);
		}
		*/
	}

	if (nonzero_trans == 0) {
		pk->trans_flag = 0;
		pk->n_trans = 0;
		pk->n_ktrans = 0;

		return md->eta; 
	}

	if (nonzero_trans == 1) {
		for (int j = 0; j < n_trans; j++) {
			if (j == argmax_xi) {
				pk->transition_p[j] = 0.0;
			}
			else {
				pk->transition_p[j] = NAN;
			}
		}

		return _k_penalty;
	}

	/* solve for lambda in the Lagrange multiplier that 
	 * imposes constraint on transition probabilities,
	 * i.e. sum_j (p.ij = 1) for any i. */
	md->xi = jnt_exp_trans;
	md->n_xi = n_trans; 

	if (argmax_xi >= 0) 
		alt_signs = pos_signs & ~(1U << argmax_xi);

	/* compute support for lambda */
	double lamb_lb = mstep_lambda_support(md);
	double l0 = sum_xi - l0_pen;
	if (l0 < 0) {
		/* negative lambda */
		l0 = lamb_lb + 1e-6;
		signs = mstep_determine_signs(lamb_lb, md);
		md->signs = signs;
	}
	else {
		signs = pos_signs;
		md->signs = pos_signs;
	}

	double lambda = mstep_newton(mstep_lagrange_cstr_func, 
			mstep_lagrange_cstr_derv, l0, 1e-10, 100, md);
	double cstr = mstep_lagrange_cstr_func(lambda, md);

	if (isnan(lambda) || cstr < -1e-1) {
		mstep_data_t d2 = {.kappa = md->kappa, .eta = md->eta,
			.xi = jnt_exp_trans, .signs = alt_signs, .n_xi = n_trans};

		lambda = mstep_try_bisection(&d2, lamb_lb, pos_signs, &signs);

		/* still can't find solution to lambda ... ? */
		if (isnan(lambda)) {
			d2.signs = pos_signs;

			/* guessing the lower bound for bisection algorithm */
			double _ub = LOG(sum_xi), _lb = LOG(1e-12);
			double _step = (_ub - _lb) / 100;
			for(; _lb < _ub; _lb += _step) {
				double _cstr = mstep_lagrange_cstr_func(EXP(_lb), &d2);
				if (_cstr > 0) break;
			}
			/* trying bisection for positive lambda */
			lambda = mstep_bisection(mstep_lagrange_cstr_func, 
					EXP(_lb), sum_xi,
					1e-10, &d2);
		}
		
		/* okay... then it must be we need to set another sign to negative,
		 * likely the one with the largest recip_tp. */
		if (isnan(lambda)) {
			mstep_data_t d2 = {.kappa = md->kappa, .eta = md->eta,
				.xi = jnt_exp_trans, .signs = signs, .n_xi = n_trans};

			for (int j = 0; j < n_trans && isnan(lambda); j++) {
				if (j == argmax_xi) continue;

				d2.signs = pos_signs & ~(1U << j);
				lambda = mstep_try_bisection(&d2, lamb_lb, pos_signs, &signs);
			}
		}

	}

	if (!isnan(lambda)) {
		double sum_tp = 0.0;
		for (int base = 0; base < n_trans; base++) {
			/* compute the MPLE for transition prob. */
			double _xi = jnt_exp_trans[base];
			double b = md->kappa * lambda - _xi + thresh; 
			double sign = (double) ((((signs >> base) & 1) << 1) - 1);

			pk->transition_p[base] = (-b + sign * 
				sqrt(b*b + 4 * lambda * md->kappa * _xi)) /
				(2 * lambda);

			sum_tp += pk->transition_p[base];
		}

		/* renormalize the estimated transition probability to avoid
		 * numerical precision issues. */
		int transidx = 0;
		kbits_t tpacked = trans_flag2packed(pk->trans_flag & 15);
		tpacked >>= 3;
		for (int j = 0; j < n_trans; j++, tpacked >>= 2) {
			kbits_t base = tpacked & 3;

			/* if certain transition has zero probability (zeroed out by
			* the penalty), we can set the transition flag to 0 to
			* speed up later E-steps. */
			if (pk->transition_p[j] == 0.0) {
				// set corresponding transition flags to 0
				pk->trans_flag &= ~(1U << base); 
				// take care of the (k+1)-kmer transitions as well
				for (int b2 = 4; b2 <= 16; b2 += 4) {
					pk->trans_flag &= ~(1U << (b2 + base));
				}
			}
			else {
				pk->transition_p[transidx++] = LOG(pk->transition_p[j] / 
						sum_tp);
			}
		}

		pk->n_trans = _mm_popcnt_u32(
				pk->trans_flag & ((1 << 20) - 1));
		pk->n_ktrans = _mm_popcnt_u32(pk->trans_flag & 15);
	}
	else {
		INFO("BAD NUMERICS!!! Kmer: %lu\n", pk->id);
	}

	return _k_penalty;
}

void mstep_emission(data *d)
{
	const size_t qmax = d->opt->max_qual_score + 1;

	for (int t = 0; t < d->max_read_length; t++) {
		size_t base_toff = t * d->base_ctx_radix;
		size_t qual_toff = t * d->qual_ctx_radix;

		for (int i = 0; i < d->n_base_context; i++) {
			double *ctx_base_emit_exp = &d->base_emit_exp[base_toff + i * 5];
			double *ctx_base_emit_p = &d->base_emission_p[base_toff + i * 5];

			double sum_emission_p = 0.0;
			for (int j = 0; j <= BASE_N; j++) {
				sum_emission_p += ctx_base_emit_exp[j];
			}

			for (int j = 0; j <= BASE_N; j++) {
				ctx_base_emit_p[j] = ctx_base_emit_exp[j] / sum_emission_p;	
			}
		}

		for (int i = 0; i < d->n_qual_context; i++) {
			double *ctx_qual_emit_exp = &d->qual_emit_exp[qual_toff + i * (qmax << 1)];
			double *ctx_qual_emit_p = &d->qual_emission_p[qual_toff + i * (qmax << 1)];

			double sum_qual_exp_c = 0.0;
			double sum_qual_exp_e = 0.0;
			for (int j = 0; j < qmax; j++) {
				register int jsll1 = j << 1;
				sum_qual_exp_e += ctx_qual_emit_exp[jsll1];
				sum_qual_exp_c += ctx_qual_emit_exp[jsll1 + 1];
			}

			for (int j = 0; j < qmax; j++) {
				register int jsll1 = j << 1;
				ctx_qual_emit_p[jsll1] = LOG(ctx_qual_emit_exp[jsll1] / 
						sum_qual_exp_e);
				ctx_qual_emit_p[jsll1 + 1] = LOG(ctx_qual_emit_exp[jsll1 + 1] / 
						sum_qual_exp_c);
			}
		}
	}
}

void mstep_initial(data *d)
{
	dictentry *de = d->kmer_dict->table;
	for (int i = d->kmer_dict->used; i > 0; de++) {
		if (de->value == NULL) continue;

		i--;
		kmer_t *pk = de->value;
		if ((pk->id >> 62) || pk->dirty) continue;

		kmer_t *rc_pk = dict_find(d->kmer_dict,
				kbits_cast_to_ptr(kmer_reverse_complement(
						pk->id, d->opt->kmer_length)))->value;

		double _init_p = LOG(
				(pk->exp_count + rc_pk->exp_count) / d->total_kmer_exp_counts);

		pk->init_p = _init_p;
		rc_pk->init_p = _init_p;

		pk->dirty = 1;
		rc_pk->dirty = 1;
	}
}

double mstep_lambda_support(void *fdata)
{
	mstep_data_t *md = fdata;
	register double _kappa = md->kappa;
	register double *_xi = md->xi;
	register double _elk = md->eta / log(1 + 1/_kappa);
	register int n_xi = md->n_xi;

	register double max_lamb_lb = -INFINITY;
	for (int j = 0; j < n_xi; j++){
		register double xi_j = _xi[j];
		register double _lb = (-(xi_j + _elk) + 2 * sqrt(xi_j * _elk)) / 
			_kappa;
		if (_lb > max_lamb_lb) 
			max_lamb_lb = _lb;
	}

	return max_lamb_lb;
}

int32_t mstep_determine_signs(double lamb_lb, void *fdata)
{
	register mstep_data_t *md = fdata;
	register double *pxi = md->xi;
	register int n_xi = md->n_xi;
	register int32_t signs_v = (1 << n_xi) - 1; 

	mstep_data_t param = {.kappa = md->kappa, .eta = md->eta,
		.xi = md->xi, .signs = signs_v};

	int argmax_idx = -1;
	double max_xi = 0.0;
	for (int j = 0; j < n_xi; j++) {
		if (pxi[j] > max_xi) {
			max_xi = pxi[j];
			argmax_idx = j;
		}
	}

	if (argmax_idx < 0) return signs_v;

	register int32_t alt_signs_v = signs_v & ~(1 << argmax_idx);

	/* FIXME: not thread-safe */
	double l0 = lamb_lb + EPSILON;
	double l1 = l0 - mstep_lagrange_cstr_func(l0, &param) / 
		mstep_lagrange_cstr_derv(l0, &param);

	//param.signs = alt_signs_v;

	double l1_alt = l0 - mstep_lagrange_cstr_func(l0, &param) / 
		mstep_lagrange_cstr_derv(l0, &param);

	if (lamb_lb < 0 && l1 > 0) 
		return signs_v;

	if (l1_alt > lamb_lb) {
		return alt_signs_v;
	}
	else {
		return signs_v;
	}
}

double mstep_lagrange_cstr_func(double lambda, void *fdata)
{
	register mstep_data_t *md = fdata;
	register double kap = md->kappa;
	register double elk = md->eta / log(1 + 1/kap);
	register double *pxi = md->xi;
	register int32_t sgns = md->signs;
	register int n_xi = md->n_xi;
	int sgn_off = n_xi > 4;

	register double sum_root = 0.0;
	for (int j = 0; j < n_xi; j++) {
		register double xij = pxi[j];

		double b = kap * lambda - xij + elk;
		double mc = kap * xij;

		double delta = b * b + 4 * lambda * mc;
		/* sanity check */
		if (delta < 0) return NAN;

		delta = sqrt(delta);

		double sign = (double) (((((sgns >> j) & 1) + sgn_off) << 1) - 1);
		sum_root += (-b + sign * delta)/(2 * lambda);
	}

	return sum_root - 1;
}

double mstep_lagrange_cstr_derv(double lambda, void *fdata)
{
	register mstep_data_t *md = fdata;
	register double kap = md->kappa;
	register double elk = md->eta / log(1 + 1/kap);
	register double *pxi = md->xi;
	register int32_t sgns = md->signs;
	register int n_xi = md->n_xi;
	int sgn_off = n_xi > 4;

	register double derv = 0.0;
	for (int j = 0; j < n_xi; j++) {
		register double xij = pxi[j];

		double b = kap * lambda - xij + elk;
		double mc = kap * xij;

		double delta = b * b + 4 * lambda * mc;
		/* sanity check */
		if (delta < 0) return NAN;

		delta = sqrt(delta);
		double sign = (double) (((((sgns >> j) & 1) + sgn_off) << 1) - 1);
		derv += (lambda * (-kap + sign / delta * kap * (b + 2 * xij))) - 
			(-b + sign * delta);
	}

	return derv / (2 * lambda * lambda);
}

/* ----- Root finding algorithms ----- */
double mstep_newton(double (*fx)(double x, void *data), 
		double (*fderv)(double x, void *data), double x0, double eps,
		int maxiter, void *fdata)
{
	int iter = 1;
	double x1 = x0 - fx(x0, fdata) / fderv(x0, fdata);

	if (isnan(x1)) return NAN;

	while (fabs(x1 - x0) > eps && iter < maxiter) {
		x0 = x1;
		x1 = x0 - fx(x0, fdata) / fderv(x0, fdata);

		if (isnan(x1)) return NAN;

		iter++;
	}

	return x1;
}

/**
 * Finding root of f(x) using bisection method.
 *
 * @param fx function f(x)
 * @param a lower bound
 * @param b upper bound
 * @param eps accepted error
 * @param fdata extra data pass to the function pointer
 */
double mstep_bisection(double (*fx)(double x, void *data), double a, double b, 
		double eps, void *fdata)
{
	double fa = fx(a, fdata),
		   fb = fx(b, fdata);

	double mid, fmid;

	if (fa == 0) return a;
	if (fb == 0) return b;

	if (a > b || (fa * fb) > 0)
		return NAN;

	mid = (a+b) / 2;

	if (b - a < eps)
		return mid;

	fmid = fx(mid, fdata);

	if ((fa < 0 && fmid >= 0) || (fa > 0 && fmid <= 0))
		return mstep_bisection(fx, a, mid, eps, fdata);
	else
		return mstep_bisection(fx, mid, b, eps, fdata);
}	/* bisection */

double mstep_try_bisection(mstep_data_t *d2, double lamb_lb,
		int32_t pos_signs, int *signs)
{
	mstep_data_t _d1;
	mstep_data_t *d1 = &_d1;
	memcpy(d1, d2, sizeof(_d1));
	d1->signs = pos_signs; 

	double lo_delta = 1e-12;
	double lambda = NAN;

	/* (\sum_j p_{ij}) - 1 at lower bound */
	double lin_cstr_lo = mstep_lagrange_cstr_func(lamb_lb - lamb_lb * lo_delta, d2);

	int iter = 0;
	while (isnan(lin_cstr_lo) && (iter++) < 60) {
		lo_delta = pow(lo_delta, 0.5);
		lin_cstr_lo = mstep_lagrange_cstr_func(lamb_lb - lamb_lb * lo_delta, d2);
	}
	
	if (isnan(lin_cstr_lo)) return NAN;

	double bound_sign = lin_cstr_lo * mstep_lagrange_cstr_func(-1e-10, d2);

	if (!isnan(bound_sign) && bound_sign < 0) {
		lambda = mstep_bisection(mstep_lagrange_cstr_func,
				lamb_lb - lamb_lb * lo_delta, -1e-10, 
				1e-10, d2);
		*signs = d2->signs;
	}
	else if (mstep_lagrange_cstr_func(lamb_lb - lamb_lb * lo_delta, d1) *
			mstep_lagrange_cstr_func(-1e-10, d1) < 0) {
		lambda = mstep_bisection(mstep_lagrange_cstr_func,
				lamb_lb - lamb_lb * lo_delta, -1e-10, 
				1e-10, d1);
		*signs = d1->signs;
	}

	return lambda;
}
