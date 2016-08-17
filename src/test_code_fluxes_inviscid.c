// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_code_fluxes_inviscid.h"

#include <stdlib.h>
#include <stdio.h>

#include "Macros.h"

/*
 *	Purpose:
 *		Provide functions for testing fluxes_inviscid related functions.
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
			printf("Error: Unsupported dimension (%d).\n",d), EXIT_MSG;
			break;
		}
	}
	return W;
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
