#ifndef __PREMIER_MSTEP_H__
#define __PREMIER_MSTEP_H__

#include <stdint.h>
#include <math.h>
#include "premier.h"
#include "numeric.h"

#ifdef __POPCNT__
#include <nmmintrin.h>
#endif

typedef struct tflag_s tflag_t;

struct tflag_s {
	uint32_t trans_flag;
	int n_nonzero_trans;
};

void mstep_transition(data *d);
void mstep_emission(data *d);
void mstep_initial(data *d);

double mstep_lambda_support(void *fdata);
int32_t mstep_determine_signs(double lamb_lb, void *fdata);
double mstep_lagrange_cstr_func(double lambda, void *fdata);
double mstep_lagrange_cstr_derv(double lambda, void *fdata);

double mstep_newton(double (*fx)(double x, void *data), 
		double (*fderv)(double x, void *data), double x0, double eps,
		int maxiter, void *fdata);
double mstep_bisection(double (*fx)(double x, void *data), double a, double b, 
		double eps, void *fdata);

double mstep_try_bisection(mstep_data_t *d2, double lamb_lb,
		int32_t pos_signs, int *signs);

#endif
