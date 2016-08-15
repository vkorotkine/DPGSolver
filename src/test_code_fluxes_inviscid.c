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
			                   2.01, -2.02,  2.03, -2.11,  2.12,  2.13,
			                   5.01,  5.02,  5.03,  5.11,  5.12,  5.13};

			for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
				W[i] = W1[i];

			break;
		} case 2: {
			double W2[24] = {  1.01,  1.02,  1.03,  1.11,  1.12,  1.13,
			                   2.01, -2.02,  2.03, -2.11,  2.12,  2.13,
			                   3.04, -3.05,  3.06, -3.14,  3.15,  3.16,
			                   5.01,  5.02,  5.03,  5.11,  5.12,  5.13};

			for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
				W[i] = W2[i];

			break;
		} case 3: {
			double W3[30] = {  1.01,  1.02,  1.03,  1.11,  1.12,  1.13,
			                   2.01, -2.02,  2.03, -2.11,  2.12,  2.13,
			                   3.04, -3.05,  3.06, -3.14,  3.15,  3.16,
			                   4.07, -4.08,  4.09, -4.17,  4.18,  4.19,
			                   5.01,  5.02,  5.03,  5.11,  5.12,  5.13};

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
			double n1[6 ] = {  0.6651,  0.8190,  0.2449,  0.7004, -0.7150, -0.1394};

			for (i = 0, iMax = NnTotal*d; i < iMax; i++)
				n[i] = n1[i];

			break;
		} case 2: {
			double n2[12] = {  0.6651,  0.8190,  0.2449,  0.7004, -0.7150, -0.1394,
			                   0.7395,  0.5670, -0.4809, -0.1144,  0.3626, -0.4143};

			for (i = 0, iMax = NnTotal*d; i < iMax; i++)
				n[i] = n2[i];

			break;
		} case 3: {
			double n3[18] = {  0.6651,  0.8190,  0.2449,  0.7004, -0.7150, -0.1394,
			                   0.7395,  0.5670, -0.4809, -0.1144,  0.3626, -0.4143,
			                   0.1037, -0.0875,  0.8419, -0.7045, -0.5978, -0.8994};

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
