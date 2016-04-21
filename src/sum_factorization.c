#include <stdlib.h>
#include <stdio.h>

#include "functions.h"

/*
 *	Purpose:
 *		Perform operations using sum factorization for Tensor-Product elements.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *		ToBeModified.
 */

void sf_operate_d(const int NOut, const int NCols, const int NIn, const int BRowMaxIn, double *OP, double *Input,
                  double *Output)
{
	/*	Purpose:
	 *		Perform the matrix-matrix operation required by sf_apply_*.
	 *
	 *	Comments:
	 *		See the '*** IMPORTANT ***' comment in sf_apply_*.
	 *		The 'register' prefix is likely unused here as a lot of work is done in the mm_* call.
	 */

	register unsigned IndIn, IndOut, stepIndIn, stepIndOut, BRowMax;

	stepIndIn  = NIn*NCols;
	stepIndOut = NOut*NCols;

	IndIn  = 0;
	IndOut = 0;
	for (BRowMax = BRowMaxIn; BRowMax--; ) {
		mm_CTN_d(NOut,NCols,NIn,OP,&Input[IndIn],&Output[IndOut]);

		IndIn  += stepIndIn;
		IndOut += stepIndOut;
	}
}

void sf_swap_d(double *Input, const unsigned int dim, const unsigned int NIn, const unsigned int step,
               const unsigned int NRows, const unsigned int NCols,
               const unsigned int iBound, const unsigned int jBound, const unsigned int kBound,
               const unsigned int iStep, const unsigned int jStep1, const unsigned int jStep2, const unsigned int kStep)
{
	/*	Purpose:
	 *		Perform the row swapping operation required by sf_apply_*.
	 *
	 *	Comments:
	 *		See the '*** IMPORTANT ***' comment in sf_apply_*.
	 *		The switch based on the (dim)ension is necessary as the swapping is done in a specific looping order in the
	 *		difference cases (Note the ordering of the i/j/k loops).
	 */

	register unsigned int i, iMax, j, jMax, k, kMax,
	                      RowInd, RowSub, ReOrder;
	unsigned int lInput[NIn], RowlInput[NRows];

	for (i = 0, iMax = NIn; iMax--; i++)
		lInput[i] = i*step;
	for (i = 0, iMax = NRows; iMax--; i++)
		RowlInput[i] = i;

	switch (dim) {
	case 2:
		for (k = 0, kMax = kBound; kMax--; k++) {
		for (i = 0, iMax = iBound; iMax--; i++) {
		for (j = 0, jMax = jBound; jMax--; j++) {

			RowInd  = i*iStep+j+k*kStep;
			ReOrder = i+lInput[j]+k*kStep;

			for (RowSub = ReOrder; RowlInput[RowSub] != ReOrder; RowSub = RowlInput[RowSub])
				;

			if (RowInd != RowSub) {
				array_swap_d(&Input[RowInd],&Input[RowSub],NCols,NRows);
				array_swap_ui(&RowlInput[RowInd],&RowlInput[RowSub],1,1);
			}
		}}}
		break;
	case 3:
		for (i = 0, iMax = iBound; iMax--; i++) {
		for (j = 0, jMax = jBound; jMax--; j++) {
		for (k = 0, kMax = kBound; kMax--; k++) {
			RowInd  = i*iStep+j*jStep1+k;
			ReOrder = i+j*jStep2+lInput[k];

			for (RowSub = ReOrder; RowlInput[RowSub] != ReOrder; RowSub = RowlInput[RowSub])
				;

			if (RowInd != RowSub) {
				array_swap_d(&Input[RowInd],&Input[RowSub],NCols,NRows);
				array_swap_ui(&RowlInput[RowInd],&RowlInput[RowSub],1,1);
			}
		}}}
		break;
	default:
		printf("Error: Invalid dimension entered in sf_swap_*.\n"), exit(1);
		break;
	}
}

void sf_apply_d(double *Input, double *Output, const unsigned int NIn[3], const unsigned int NOut[3],
                const unsigned int NCols, double *OP[3], const unsigned int Diag[3], const unsigned int d)
{
	/*
	 *	Purpose:
	 *		Use TP sum factorized operators to speed up calculations.
	 *		LIKELY HAVE THIS FUNCTION ACT ON BOTH SIMPLEX AND TP ELEMENTS IN FUTURE. (ToBeDeleted)
	 *
	 *	Comments:
	 *		*** IMPORTANT ***
	 *		It is assumed that the operators are stored in row major layout, while the input is stored in a column major
	 *		layout. This has several advantages despite the obvious disadvantage of being atypical and potentially
	 *		confusing:
	 *			1) Extremely efficient matrix-matrix multiplication in the vectorized version of the code (minimizing
	 *			   memory stride) when noting that the operator matrices generally fit completely in the cache
	 *			   (ToBeModified after confirming this with testing).
	 *			2) In vectorized version of the code, when multiple elements are operated on at once, performing blas
	 *			   calls using the CblasColMajor layout results in the output being stored in continuous memory for each
	 *			   element. This results in greatly reduced memory stride when distributing the output back to the
	 *			   VOLUME/FACET structures.
	 *		*** IMPORTANT ***
	 *
	 *		For the moment, the routine is only implemented using the new non-redundant approach. (ToBeModified)
	 *		Operating in the eta/zeta directions requires re-ordering of the matrices before and after operation. To
	 *		minimize memory usage, re-ordering is done in place, requiring only a single additional row of storage.
	 *		Add support for NonRedundant == 2 (Diag == 1). (ToBeDeleted)
	 *			Likely implementation: extract diagonal from OP, then loop over rows performing BLAS 1 row scaling. Note
	 *			                       that this really does not require re-arranging. If it is found that this option
	 *			                       is important in the future, profile both implementations. (ToBeDeleted)
	 *		Add the implementation for the further 2x reduction through fourier transform of operators (ToBeDeleted).
	 *			This will only be made available for major operators as it requires storage of the decomposed operators.
	 *		Make sure that appropriate variables are declared with 'register' and 'unsigned' (ToBeDeleted).
	 *		If found to be slow during profiling, change all array indexing operations to pointer operations where
	 *		possible. Also change loops so that exit check decrements to 0 instead of using comparison (ToBeDeleted).
	 *
	 *	Notation:
	 *		Input  : Input array, number of entries prod(NIn) x NCols
	 *		Output : Output array, number of entries prod(NOut) x NCols
	 *		N()[]  : (N)umber of (In/Out)put entries in each of the coordinate directions []
	 *		OP[]   : 1D operators in each of the coordinate directions []
	 *		Diag   : Indication of whether the OPs are diagonal
	 *		         Options: 0 (Not diagonal)
	 *		                  1 (Diagonal but not identity)
	 *		                  2 (Diagonal identity)
	 *
	 *	References:
	 *		Add in Sherwin's book or perhaps my thesis as the procedure implemented is slighly modified (ToBeModified).
	 */

	register unsigned int i, iMax, dim;
	unsigned int          NonRedundant[d], BRows[d], NRows_Out[d], Indd;
	double                **Output_Inter;

	for (dim = 0; dim < d; dim++) {
		if      (Diag[dim] == 0) NonRedundant[dim] = 1;
		else if (Diag[dim] == 1) NonRedundant[dim] = 2;
		else if (Diag[dim] == 2) NonRedundant[dim] = 0;
		else
			printf("Error: Invalid entry in Diag.\n"), exit(1);
	}

	BRows[0] = NIn[1]*NIn[2];
	if (d > 1)
		BRows[1] = NOut[0]*NIn[2];
	if (d > 2)
		BRows[2] = NOut[0]*NOut[1];

	for (dim = 0; dim < d; dim++)
		NRows_Out[dim] = NOut[dim]*BRows[dim];

	Output_Inter = malloc(d * sizeof *Output); // free

	// xi
	Indd = 0;

	if (d == 1)
		Output_Inter[Indd] = Output;
	else
		Output_Inter[Indd] = malloc(NRows_Out[Indd]*NCols * sizeof *Output_Inter[Indd]); // free

	if (NonRedundant[Indd]) {
		sf_operate_d(NOut[Indd],NCols,NIn[Indd],BRows[Indd],OP[Indd],Input,Output_Inter[Indd]);
	} else {
		for (i = 0, iMax = NRows_Out[Indd]*NCols; i < iMax; i++)
			Output_Inter[Indd][i] = Input[i];
	}
//array_print_d(NRows_Out[Indd],NCols,Output_Inter[Indd],'C');

	if (d == 1) {
		free(Output_Inter);
		return;
	}

	// eta
	Indd = 1;

	if (d == 2)
		Output_Inter[Indd] = Output;
	else
		Output_Inter[Indd] = malloc(NRows_Out[Indd]*NCols * sizeof *Output_Inter[Indd]); // free

	if (NonRedundant[Indd]) {
		sf_swap_d(Output_Inter[Indd-1],Indd+1,NIn[Indd],NOut[0],NRows_Out[Indd-1],NCols,
		          NOut[0],NIn[1],NIn[2],NIn[1],0,0,NIn[1]*NOut[0]);

		sf_operate_d(NOut[Indd],NCols,NIn[Indd],BRows[Indd],OP[Indd],Output_Inter[Indd-1],Output_Inter[Indd]);

		sf_swap_d(Output_Inter[Indd],Indd+1,NOut[Indd],NOut[0],NRows_Out[Indd],NCols,
		          NOut[0],NOut[1],NIn[2],NOut[1],0,0,NOut[1]*NOut[0]);
	} else {
		for (i = 0, iMax = NRows_Out[Indd]*NCols; i < iMax; i++)
			Output_Inter[Indd][i] = Output_Inter[Indd-1][i];
	}
	free(Output_Inter[Indd-1]);
//array_print_d(NRows_Out[Indd],NCols,Output_Inter[Indd],'C');

	if (d == 2) {
		free(Output_Inter);
		return;
	}

	// zeta
	Indd = 2;

	Output_Inter[Indd] = Output;

	if (NonRedundant[Indd]) {
		sf_swap_d(Output_Inter[Indd-1],Indd+1,NIn[Indd],NOut[0]*NOut[1],NRows_Out[Indd-1],NCols,
		          NOut[0],NOut[1],NIn[2],NOut[1]*NIn[2],NIn[2],NOut[0],0);

		sf_operate_d(NOut[Indd],NCols,NIn[Indd],BRows[Indd],OP[Indd],Output_Inter[Indd-1],Output_Inter[Indd]);

		sf_swap_d(Output_Inter[Indd],Indd+1,NOut[Indd],NOut[0]*NOut[1],NRows_Out[Indd],NCols,
		          NOut[0],NOut[1],NOut[2],NOut[1]*NOut[2],NOut[2],NOut[0],0);
	} else {
		for (i = 0, iMax = NRows_Out[Indd]*NCols; i < iMax; i++)
			Output_Inter[Indd][i] = Output_Inter[Indd-1][i];
	}
	free(Output_Inter[Indd-1]);
//array_print_d(NRows_Out[Indd],NCols,Output_Inter[Indd],'C');

	free(Output_Inter);
}

