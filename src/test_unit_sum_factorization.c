// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_sum_factorization.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"
#include "Test.h"

#include "test_support.h"
#include "array_norm.h"
#include "sum_factorization.h"

/*
 *	Purpose:
 *		Test correctness of implementation of sum_factorization functions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_unit_sum_factorization(void)
{
	unsigned int pass;

	/*
	 *	sf_apply_d:
	 *
	 *	Input:
	 *
	 *		NIn  = [ 7 5 3 ]
	 *		NOut = [ 6 4 2 ]
	 *		Input  = [ Equispaced nodes on HEX of order NIn-1 in each dimension ]
	 *
	 *	Expected Output:
	 *
	 *		Output = [ Equispaced nodes on HEX of order NOut-1 in each dimension ]
	 *
	 */

	unsigned int i, j, k, iMax, jMax, kMax, dim, d, Ind,
	             NIn_Total[3] = {105, 90, 72}, NOut_Total[3] = {90, 72, 48}, NCols, Diag[3];
	unsigned int NIn[3] = { 7, 5, 3 }, NOut[3] = { 6, 4, 2 };
	double *Output1_c, *Output2_c, *Output3_c, *OP[3];
	double Input[315], Output1[270], Output2[216], Output3[144];

	d = 3;

	Ind = 0;
	for (k = 0, kMax = NIn[2]; k < kMax; k++) {
	for (j = 0, jMax = NIn[1]; j < jMax; j++) {
	for (i = 0, iMax = NIn[0]; i < iMax; i++) {
		Input[NIn_Total[0]*0+Ind] = -1.0+(2.0*i)/(NIn[0]-1);
		Input[NIn_Total[0]*1+Ind] = -1.0+(2.0*j)/(NIn[1]-1);
		Input[NIn_Total[0]*2+Ind] = -1.0+(2.0*k)/(NIn[2]-1);
		Ind++;
	}}}

	Ind = 0;
	for (k = 0, kMax = NIn[2]; k < kMax; k++) {
	for (j = 0, jMax = NIn[1]; j < jMax; j++) {
	for (i = 0, iMax = NOut[0]; i < iMax; i++) {
		Output1[NOut_Total[0]*0+Ind] = -1.0+(2.0*i)/(NOut[0]-1);
		Output1[NOut_Total[0]*1+Ind] = -1.0+(2.0*j)/(NIn[1]-1);
		Output1[NOut_Total[0]*2+Ind] = -1.0+(2.0*k)/(NIn[2]-1);
		Ind++;
	}}}

	Ind = 0;
	for (k = 0, kMax = NIn[2]; k < kMax; k++) {
	for (j = 0, jMax = NOut[1]; j < jMax; j++) {
	for (i = 0, iMax = NOut[0]; i < iMax; i++) {
		Output2[NOut_Total[1]*0+Ind] = -1.0+(2.0*i)/(NOut[0]-1);
		Output2[NOut_Total[1]*1+Ind] = -1.0+(2.0*j)/(NOut[1]-1);
		Output2[NOut_Total[1]*2+Ind] = -1.0+(2.0*k)/(NIn[2]-1);
		Ind++;
	}}}

	Ind = 0;
	for (k = 0, kMax = NOut[2]; k < kMax; k++) {
	for (j = 0, jMax = NOut[1]; j < jMax; j++) {
	for (i = 0, iMax = NOut[0]; i < iMax; i++) {
		Output3[NOut_Total[2]*0+Ind] = -1.0+(2.0*i)/(NOut[0]-1);
		Output3[NOut_Total[2]*1+Ind] = -1.0+(2.0*j)/(NOut[1]-1);
		Output3[NOut_Total[2]*2+Ind] = -1.0+(2.0*k)/(NOut[2]-1);
		Ind++;
	}}}

	for (dim = 0; dim < d; dim++) {
		OP[dim] = malloc(NIn[dim]*NOut[dim] * sizeof *OP[dim]); // free

		Ind = 0;
		for (i = 0, iMax = NOut[dim]; i < iMax; i++) {
		for (j = 0, jMax = NIn[dim]; j < jMax; j++) {
			if (iMax-1 != 0) {
				if (j == 0)
					OP[dim][i*jMax+j] = (iMax-i-1)/((double) iMax-1);
				else if (j == jMax-1)
					OP[dim][i*jMax+j] = (i)/((double) iMax-1);
				else
					OP[dim][i*jMax+j] = 0.0;
			} else {
				OP[dim][i*jMax+j] = 0.5;
			}
		}}
	}

	for (dim = 0; dim < d; dim++)
		Diag[dim] = 0;

	Output1_c = malloc(NOut_Total[0]*d * sizeof *Output1_c); // free
	Output2_c = malloc(NOut_Total[1]*d * sizeof *Output2_c); // free
	Output3_c = malloc(NOut_Total[2]*d * sizeof *Output3_c); // free

	NCols = d*1;


	// r component
	sf_apply_d(Input,Output1_c,NIn,NOut,NCols,OP,Diag,1);

	pass = 0;
	if (array_norm_diff_d(NOut_Total[0]*d,Output1_c,Output1,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"sum_factorization (d = 1):");

	// r, s components
	sf_apply_d(Input,Output2_c,NIn,NOut,NCols,OP,Diag,2);

	pass = 0;
	if (array_norm_diff_d(NOut_Total[1]*d,Output2_c,Output2,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"                  (d = 2):");

	// r, s, t components
	sf_apply_d(Input,Output3_c,NIn,NOut,NCols,OP,Diag,3);

	pass = 0;
	if (array_norm_diff_d(NOut_Total[2]*d,Output3_c,Output3,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"                  (d = 3):");

	for (dim = 0; dim < d; dim++)
		free(OP[dim]);
	free(Output1_c);
	free(Output2_c);
	free(Output3_c);
}
