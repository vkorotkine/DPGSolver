// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_code_fluxes.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "Test.h"

#include "test_support.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions for testing flux related functions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

double *initialize_W(unsigned int *Nn, unsigned int *Nel, const unsigned int d)
{
	unsigned int i, iMax, NnTotal, Nvar;
	double       *W;

	*Nn  = 3;
	*Nel = 2;

	NnTotal = (*Nn)*(*Nel);
	Nvar = d+2;

	W = malloc(NnTotal*Nvar * sizeof *W); // keep

	switch (d) {
		case 1: {
			double W1[18] = {  1.01,  1.02,  1.03,  1.11,  1.12,  1.13,
			                   2.01, -2.02,  2.03, -2.11,  2.12,  2.53,
			                   9.01,  9.02,  9.03,  9.11,  9.12,  9.13};

			for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
				W[i] = W1[i];

			break;
		} case 2: {
			double W2[24] = {  1.01,  1.02,  1.03,  1.11,  1.12,  1.13,
			                   2.01, -2.22,  2.03, -2.11,  2.12,  2.53,
			                   2.04, -2.25,  2.06, -2.14,  2.15,  2.56,
			                   9.01,  9.02,  9.03,  9.11,  9.12,  9.13};

			for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
				W[i] = W2[i];

			break;
		} case 3: {
			double W3[30] = {  1.01,  1.02,  1.03,  1.11,  1.12,  1.13,
			                   2.41, -2.02,  1.92, -2.11,  2.12,  2.13,
			                   2.44, -2.05,  1.95, -2.14,  2.15,  2.16,
			                   2.47, -2.08,  1.98, -2.17,  2.18,  2.19,
			                   9.01,  9.02,  9.02,  9.11,  9.12,  9.13};

			for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
				W[i] = W3[i];

			break;
		} default: {
			EXIT_UNSUPPORTED;
			break;
		}
	}
	return W;
}

double **initialize_Q(unsigned int const Nn, unsigned int const Nel, unsigned int const d)
{
	// Numbers used here are random.

	unsigned int const NnTotal = Nn*Nel,
	                   Nvar    = d+2;

	double **Q = malloc(d * sizeof *Q); // keep
	for (size_t dim = 0; dim < d; dim++)
		Q[dim] = malloc(NnTotal*Nvar * sizeof *Q[dim]); // keep

	switch (d) {
		case 2: {
			double Q2[2][24] = {{ -0.0855,  -0.9289,   0.2373,  -0.5211,   0.6791,  -0.0377,
			                      -0.2625,   0.7303,   0.4588,  -0.2316,   0.3955,  -0.8852,
			                      -0.8010,  -0.4886,  -0.9631,   0.4889,  -0.3674,  -0.9133,
			                      -0.0292,  -0.5785,  -0.5468,  -0.6241,  -0.9880,   0.7962, },
			                    { -0.4018,   0.1839,  -0.9027,  -0.3377,   0.7803,   0.0965,
			                       0.0760,   0.2400,   0.9448,  -0.9001,   0.3897,   0.1320,
			                      -0.2399,  -0.4173,  -0.4909,  -0.3692,   0.2417,  -0.9421,
			                       0.1233,  -0.0497,   0.4893,  -0.1112,  -0.4039,   0.9561, }, };

			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0; i < NnTotal*Nvar; i++)
					Q[dim][i] = Q2[dim][i];
			}
			break;
		} case 3: {
			double Q3[3][30] = {{ -0.0911,  -0.6444,  -0.2089,  -0.4501,  -0.6620,   0.6135,
			                       0.5762,   0.6476,   0.7093,  -0.4587,   0.4162,  -0.5822,
			                       0.6834,  -0.6790,  -0.2362,   0.6619,  -0.8419,   0.5407,
			                      -0.5466,  -0.6358,   0.1194,  -0.7703,   0.8329,  -0.8699,
			                       0.4257,  -0.9452,  -0.6073,  -0.3502,  -0.2564,  -0.2648, },
			                    {  0.7184,  -0.6110,   0.1537,  -0.8754,   0.2407,  -0.0680,
			                       0.9686,   0.7788,  -0.2810,  -0.5181,   0.6761,   0.2548,
			                       0.5313,   0.4235,  -0.4401,  -0.9436,   0.2891,   0.2240,
			                      -0.3251,   0.0908,  -0.5271,   0.6377,  -0.6718,  -0.6678,
			                       0.1056,  -0.2665,  -0.4574,  -0.9577,  -0.6951,   0.8444, },
			                    { -0.4899,   0.4711,   0.5216,  -0.1499,   0.8003,   0.1332,
			                      -0.1679,  -0.0596,   0.0967,  -0.6596,   0.4538,  -0.1734,
			                      -0.9787,   0.6820,   0.8181,   0.5186,   0.4324,  -0.3909,
			                       0.7127,  -0.0424,  -0.8175,   0.9730,   0.8253,   0.8314,
			                      -0.5005,   0.0714,   0.7224,  -0.6490,   0.0835,  -0.8034, },};

			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0; i < NnTotal*Nvar; i++)
					Q[dim][i] = Q3[dim][i];
			}
			break;
		} default: {
			EXIT_UNSUPPORTED;
			break;
		}
	}
	return Q;
}

double *initialize_n(const unsigned int Nn, const unsigned int Nel, const unsigned int d)
{
	unsigned int i, iMax, NnTotal;
	double       *n;

	NnTotal = Nn*Nel;
	if (NnTotal != 6)
		printf("Error: Unsupported NnTotal.\n"), EXIT_MSG;

	n = malloc(NnTotal*d * sizeof *n); // keep

	switch (d) {
		case 1: {
			double n1[6 ] = {  0.6651,
			                   0.8190,
			                   0.2449,
			                   0.7004,
			                  -0.7150,
			                  -0.1394};

			for (i = 0, iMax = NnTotal*d; i < iMax; i++)
				n[i] = n1[i];

			break;
		} case 2: {
			double n2[12] = {  0.6651,  0.7395,
			                   0.8190,  0.5670,
			                   0.2449, -0.4809,
			                   0.7004, -0.1144,
			                  -0.7150,  0.3626,
			                  -0.1394, -0.4143};

			for (i = 0, iMax = NnTotal*d; i < iMax; i++)
				n[i] = n2[i];

			break;
		} case 3: {
			double n3[18] = {  0.6651,  0.7395,  0.1037,
			                   0.8190,  0.5670, -0.0875,
			                   0.8190,  0.5670, -0.0875,
			                   0.7004, -0.1144, -0.7045,
			                  -0.7150,  0.3626, -0.5978,
			                  -0.1394, -0.4143, -0.8994};

			for (i = 0, iMax = NnTotal*d; i < iMax; i++)
				n[i] = n3[i];

			break;
		} default: {
			printf("Error: Unsupported dimension (%d).\n",d), EXIT_MSG;
			break;
		}
	}
	return n;
}

double *initialize_XYZ(const unsigned int Nn, const unsigned int Nel, const unsigned int d)
{
	unsigned int i, iMax, NnTotal;
	double       *XYZ;

	NnTotal = Nn*Nel;
	if (NnTotal != 6)
		printf("Error: Unsupported NnTotal.\n"), EXIT_MSG;

	XYZ = malloc(NnTotal*d * sizeof *XYZ); // keep


	switch (d) {
		case 1: {
			double XYZ1[6 ] = {  1.8147,  1.8568,  0.9360,  1.5012,  0.9414,  0.0975};

			for (i = 0, iMax = NnTotal*d; i < iMax; i++)
				XYZ[i] = XYZ1[i];

			break;
		} case 2: {
			double XYZ2[12] = {  1.8147,  1.8568,  0.9360,  1.5012,  0.9414,  0.0975,
			                     0.2785,  0.8559,  1.5453,  1.7739,  1.1087,  1.9706};

			for (i = 0, iMax = NnTotal*d; i < iMax; i++)
				XYZ[i] = XYZ2[i];

			break;
		} case 3: {
			double XYZ3[18] = {  1.8147,  1.8568,  0.9360,  1.5012,  0.9414,  0.0975,
			                     0.2785,  0.8559,  1.5453,  1.7739,  1.1087,  1.9706,
			                     0.1143, -0.0037,  0.0751, -0.0895, -0.0196,  0.1039};

			for (i = 0, iMax = NnTotal*d; i < iMax; i++)
				XYZ[i] = XYZ3[i];

			break;
		} default: {
			printf("Error: Unsupported dimension (%d).\n",d), EXIT_MSG;
			break;
		}
	}
	return XYZ;
}

void set_FiTypes(unsigned int *NFiTypesOut, char ***FiTypeOut)
{
	unsigned int const NFiTypes = 2;

	char **FiType = malloc(NFiTypes * sizeof *FiType); // keep
	for (size_t i = 0; i < NFiTypes; i++)
		FiType[i] = malloc(STRLEN_MAX * sizeof *FiType[i]); // keep

	strcpy(FiType[0],"Advection");
	strcpy(FiType[1],"Euler    ");

	*NFiTypesOut = NFiTypes;
	*FiTypeOut   = FiType;
}

void set_FNumTypes(unsigned int *NFNumTypesOut, char ***FNumTypeOut)
{
	unsigned int const NFNumTypes = 3;

	char **FNumType = malloc(NFNumTypes * sizeof *FNumType); // keep
	for (size_t i = 0; i < NFNumTypes; i++)
		FNumType[i] = malloc(STRLEN_MAX * sizeof *FNumType[i]); // keep

	strcpy(FNumType[0],"LF     ");
	strcpy(FNumType[1],"RoePike");
	strcpy(FNumType[2],"Upwind ");

	*NFNumTypesOut = NFNumTypes;
	*FNumTypeOut   = FNumType;
}

void set_parameters_test_flux_inviscid(char const *const FiType, unsigned int const d)
{
	DB.d = d;
	if (strstr(FiType,"Advection")) {
		strcpy(DB.PDE,"Advection");
		strcpy(DB.PDESpecifier,"Steady/Default");
		strcpy(DB.Geometry,"n-Cube");
		strcpy(DB.GeomSpecifier,"YL");
	} else if (strstr(FiType,"Euler")) {
		strcpy(DB.PDE,"Euler");
		strcpy(DB.TestCase,"SupersonicVortex");
		strcpy(DB.PDESpecifier,"Internal");
		strcpy(DB.Geometry,"n-Cylinder_HollowSection");
	} else {
		printf("Unsupported FiType: %s\n",FiType);
		EXIT_UNSUPPORTED;
	}
}

void set_parameters_test_flux_Num(char const *const FNumType, unsigned int const d)
{
	DB.d = d;
	if (strstr(FNumType,"Upwind")) {
		strcpy(DB.PDE,"Advection");
		strcpy(DB.PDESpecifier,"Steady/Default");
		strcpy(DB.Geometry,"n-Cube");
		strcpy(DB.GeomSpecifier,"YL");

		DB.InviscidFluxType = FLUX_UPWIND;
	} else if (strstr(FNumType,"LF") || strstr(FNumType,"RoePike")) {
		strcpy(DB.PDE,"Euler");
		strcpy(DB.TestCase,"SupersonicVortex");
		strcpy(DB.PDESpecifier,"Internal");
		strcpy(DB.Geometry,"n-Cylinder_HollowSection");

		if (strstr(FNumType,"LF"))
			DB.InviscidFluxType = FLUX_LF;
		else if (strstr(FNumType,"RoePike"))
			DB.InviscidFluxType = FLUX_ROE;
		else
			EXIT_UNSUPPORTED;
	} else {
		printf("Unsupported FNumType: %s\n",FNumType);
		EXIT_UNSUPPORTED;
	}
}

void reset_entered_test_num_flux(char const *const FNumType)
{
	if (strstr(FNumType,"LF")) {
		for (size_t i = 0; i < TEST_N_LF; i++)
			TestDB.EnteredLF[i] = 0;
	} else if (strstr(FNumType,"RoePike")) {
		for (size_t i = 0; i < TEST_N_ROE; i++)
			TestDB.EnteredRoe[i] = 0;
	}
}
