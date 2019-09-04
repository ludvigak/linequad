#ifndef LINEQUAD_H
#define LINEQUAD_H
#include "complex.h"


#define MAX_EXPANSION_ORDER 64

void rsqrt_pow_weights(const double* restrict tj,
		       double complex troot,
		       int n,
		       int pmax,
		       double* restrict w1,
		       double* restrict w3,
		       double* restrict w5);

void rsqrt_pow_integrals(double complex z,
			 int n,
			 int pmax,
			 double* restrict I1,
			 double* restrict I3,
			 double* restrict I5);

int rootfinder(const double* restrict xhat,
	       const double* restrict yhat,
	       const double* restrict zhat,
	       int n,
	       double x0,
	       double y0,
	       double z0,
	       double complex tinit,
	       double complex* troot);

double complex rootfinder_initial_guess(const double* restrict tj,
					const double* restrict xj,
					const double* restrict yj,
					const double* restrict zj,
					int n,
					double x0,
					double y0,
					double z0);


#endif /* LINEQUAD_H */
