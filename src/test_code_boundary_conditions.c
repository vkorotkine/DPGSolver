// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_code_boundary_conditions.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "Test.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions for testing relating to boundary conditions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void set_BTypes(unsigned int *NBTypesOut, char ***BTypeOut)
{
	unsigned int const NBTypes = 12;

	char **BType = malloc(NBTypes * sizeof *BType); // keep

	for (size_t i = 0; i < NBTypes; i++)
		BType[i] = malloc(STRLEN_MAX * sizeof *BType[i]); // keep

	strcpy(BType[0],"SlipWall         ");
	strcpy(BType[1],"Riemann          ");
	strcpy(BType[2],"BackPressure     ");
	strcpy(BType[3],"Total_TP         ");
	strcpy(BType[4],"SupersonicIn     ");
	strcpy(BType[5],"SupersonicOut    ");
	strcpy(BType[6],"NoSlip_Dirichlet ");
	strcpy(BType[7],"NoSlip_Adiabatic ");
	strcpy(BType[8],"Poisson_Dirichlet");
	strcpy(BType[9],"Poisson_Neumann  ");
	strcpy(BType[10],"Advection_Inflow ");
	strcpy(BType[11],"Advection_Outflow");

	*NBTypesOut = NBTypes;
	*BTypeOut   = BType;
}

void set_parameters_test_boundary_conditions(char const *const BType, unsigned int const d)
{
	DB.d = d;
	if (strstr(BType,"SlipWall") || strstr(BType,"Riemann")) {
		strcpy(DB.PDE,"Euler");
		strcpy(DB.TestCase,"SupersonicVortex");
		strcpy(DB.PDESpecifier,"Internal");
		strcpy(DB.Geometry,"n-Cylinder_HollowSection");
	} else if (strstr(BType,"BackPressure") || strstr(BType,"Total_TP")) {
		strcpy(DB.PDE,"Euler");
		strcpy(DB.TestCase,"InviscidChannel");
		strcpy(DB.PDESpecifier,"Internal/Subsonic");
		strcpy(DB.Geometry,"EllipsoidalSection");
		strcpy(DB.GeomSpecifier,"Annular/3/");
	} else if (strstr(BType,"SupersonicIn") || strstr(BType,"SupersonicOut")) {
		strcpy(DB.PDE,"Euler");
		strcpy(DB.TestCase,"InviscidChannel");
		strcpy(DB.PDESpecifier,"Internal/Supersonic");
		strcpy(DB.Geometry,"EllipsoidalSection");
		strcpy(DB.GeomSpecifier,"Annular/3/");
	} else if (strstr(BType,"NoSlip")) {
		strcpy(DB.PDE,"NavierStokes");
		strcpy(DB.TestCase,"TaylorCouette");
		strcpy(DB.Geometry,"n-Cylinder_Hollow");
	} else if (strstr(BType,"Poisson")) {
		strcpy(DB.PDE,"Poisson");
		strcpy(DB.TestCase,"Poisson");
		strcpy(DB.Geometry,"n-Ellipsoid");
		strcpy(DB.GeomSpecifier,"AR_3");
	} else if (strstr(BType,"Advection")) {
		strcpy(DB.PDE,"Advection");
		strcpy(DB.PDESpecifier,"Steady/Default");
		strcpy(DB.Geometry,"n-Cube");
		strcpy(DB.GeomSpecifier,"YL");
	} else {
		EXIT_UNSUPPORTED;
	}
}

void reset_entered_test_boundary_conditions(char const *const BType)
{
	if (strstr(BType,"Riemann")) {
		for (size_t j = 0; j < TEST_N_RIEMANN; j++)
			TestDB.EnteredRiemann[j] = 0;
	} else if (strstr(BType,"BackPressure")) {
		for (size_t j = 0; j < TEST_N_BACKPRESSURE; j++)
			TestDB.EnteredBackPressure[j] = 0;
	}
}

void update_values_BackPressure(unsigned int const Nn, unsigned int const Nel, double *const W, double *const nL,
                                unsigned int const d)
{
	unsigned int const NnTotal = Nn*Nel;

	if (NnTotal != 6)
		EXIT_UNSUPPORTED;

	if (d == 3) {
		unsigned int FinalIndices[6] = {0,2,2,3,0,6};
		for (size_t n = 0; n < NnTotal; n++) {
			if (FinalIndices[n] != n) {
				for (size_t dim = 0; dim < d; dim++)
					nL[n*d+dim] *= -1.0;
				if (FinalIndices[n] == n+1) {
					for (size_t var = 1; var < d+1; var++)
						W[n+var*NnTotal] *= 0.5;
				}
			}
		}
	} else if (d == 2) {
		unsigned int FinalIndices[6] = {0,2,0,4,0,6};
		for (size_t n = 0; n < NnTotal; n++) {
			if (FinalIndices[n] != n) {
				for (size_t dim = 0; dim < d; dim++)
					nL[n*d+dim] *= -1.0;
				if (FinalIndices[n] == n+1) {
					for (size_t var = 1; var < d+1; var++)
						W[n+var*NnTotal] *= 0.5;
				}
			}
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void check_entered_test_boundary_conditions(bool *CheckedAll, char const *const BType)
{
	*CheckedAll = 1;
	if (strstr(BType,"Riemann")) {
		for (size_t i = 0; i < TEST_N_RIEMANN; i++) {
			if (!TestDB.EnteredRiemann[i]) {
				*CheckedAll = 0;
				break;
			}
		}
		if (!(*CheckedAll))
			array_print_ui(1,4,TestDB.EnteredRiemann,'R');
	} else if (strstr(BType,"BackPressure")) {
		for (size_t i = 0; i < TEST_N_BACKPRESSURE; i++) {
			if (!TestDB.EnteredBackPressure[i]) {
				*CheckedAll = 0;
				break;
			}
		}
		if (!(*CheckedAll))
			array_print_ui(1,2,TestDB.EnteredBackPressure,'R');
	}
}
