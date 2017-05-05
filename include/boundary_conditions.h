// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__boundary_conditions_h__INCLUDED
#define DPG__boundary_conditions_h__INCLUDED

#include <stdbool.h>
#include <complex.h>

struct S_BC {
	bool         ComputeQ;
	unsigned int d, Nn, Nel, BC;

	double const *XYZ, *nL;

	double const *WL,
	             *const *QL;
	double       *WB,
	             *const *QB;

	// Used for complex step verification
	double complex const *WL_c,
	                     *const *QL_c;
	double complex       *WB_c,
	                     *const *QB_c;

	// Used for linearization
	double *dWBdWL,
	       *dQBdWL;
};

// Make these DB parameters if used regularly.
#define EXACT_SLIPWALL 0
#define EXACT_NORMAL   0

extern void   compute_exact_boundary_solution (struct S_BC *const BCdata);
extern double *compute_exact_boundary_normal  (struct S_BC *const BCdata);
extern double *compute_XYZ_boundary           (struct S_BC *const BCdata);

extern void compute_boundary_values (struct S_BC *const BCdata);
extern void get_boundary_values     (const double X, const double Y, double *const rho, double *const u,
                                     double *const v, double *const w, double *const p);

#endif // DPG__boundary_conditions_h__INCLUDED
