// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__boundary_conditions_h__INCLUDED
#define DPG__boundary_conditions_h__INCLUDED

#include <complex.h>

struct S_BC {
	unsigned int d, Nn, Nel;

	double const *XYZ, *nL;

	double const *WL,
	             *const *GradWL;
	double       *WB,
	             *const *GradWB;

	// Used for complex step verification
	double complex const *WL_c,
	                     *const *GradWL_c;
	double complex       *WB_c,
	                     *const *GradWB_c;

	// Used for linearization
	double *dWBdWL,
	       *dGradWBdWL;
};

extern void get_boundary_values        (const double X, const double Y, double *const rho, double *const u,
                                        double *const v, double *const w, double *const p);
extern void boundary_Riemann           (const unsigned int Nn, const unsigned int Nel, const double *const XYZ,
                                        const double *const WL, double *const WOut, double *const WB,
                                        const double *const nL, const unsigned int d);
extern void boundary_SlipWall          (const unsigned int Nn, const unsigned int Nel, const double *const WL,
                                        double *const WB, const double *const nL, const unsigned int d);
extern void boundary_BackPressure      (const unsigned int Nn, const unsigned int Nel, const double *const WL,
                                        double *const WB, const double *const nL, const unsigned int d,
                                        const unsigned int Neq);
extern void boundary_Total_TP          (const unsigned int Nn, const unsigned int Nel, const double *const XYZ,
                                        const double *const WL, double *const WB, const double *const nL,
                                        const unsigned int d, const unsigned int Nvar);
extern void boundary_SupersonicInflow  (const unsigned int Nn, const unsigned int Nel, const double *const XYZ,
                                        const double *const WL, double *const WB, const double *const nL,
                                        const unsigned int d, const unsigned int Nvar);
extern void boundary_SupersonicOutflow (const unsigned int Nn, const unsigned int Nel, const double *const XYZ,
                                        const double *const WL, double *const WB, const double *const nL,
                                        const unsigned int d, const unsigned int Nvar);
extern void boundary_NoSlip_Dirichlet  (struct S_BC *const BCdata);
extern void boundary_NoSlip_Adiabatic  (struct S_BC *const BCdata);

#endif // DPG__boundary_conditions_h__INCLUDED
