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

void sf_swap_d(double *Input, const unsigned int NRows, const unsigned int NCols,
               const unsigned int iBound, const unsigned int jBound, const unsigned int kBound,
               const unsigned int iStep, const unsigned int jStep, const unsigned int kStep)
{
	/*	Purpose:
	 *		Perform the row swapping operation required by sf_apply_*.
	 *
	 *	Comments:
	 *		See the '*** IMPORTANT ***' comment in sf_apply_*.
	 */

	register unsigned int i, iMax, j, jMax, k, kMax,
	                      RowInd, RowSub, ReOrder;
	unsigned int RowlInput[NRows];

	for (i = 0, iMax = NRows; iMax--; i++)
		RowlInput[i] = i;

	RowInd = 0;
	for (i = 0, iMax = iBound; iMax--; i++) {
	for (j = 0, jMax = jBound; jMax--; j++) {
	for (k = 0, kMax = kBound; kMax--; k++) {
		ReOrder = i*iStep+j*jStep+k*kStep;
//printf("%d %d %d %d %d\n",i,j,k,RowInd,ReOrder);

		for (RowSub = ReOrder; RowlInput[RowSub] != ReOrder; RowSub = RowlInput[RowSub])
			;

		if (RowInd != RowSub) {
			array_swap_d(&Input[RowInd],&Input[RowSub],NCols,NRows);
			array_swap_ui(&RowlInput[RowInd],&RowlInput[RowSub],1,1);
		}
		RowInd++;
	}}}
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
	 *
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
	 *
	 *		After the swapping is performed, note that applying the operator in a loop over the blocks of Input is the
	 *		same as interpretting the results as applying the operator to a matrix where each column is a block.
	 *
	 *		*** IMPORTANT ***
	 *
	 *		Operating in the s/t directions requires re-ordering of the matrices before and after operation. To minimize
	 *		memory usage, re-ordering is done in place, requiring only a single additional row of storage.
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

	// r
	Indd = 0;

	if (d == 1)
		Output_Inter[Indd] = Output;
	else
		Output_Inter[Indd] = malloc(NRows_Out[Indd]*NCols * sizeof *Output_Inter[Indd]); // free

	if (NonRedundant[Indd]) {
		mm_CTN_d(NOut[Indd],NCols*BRows[Indd],NIn[Indd],OP[Indd],Input,Output_Inter[Indd]);
	} else {
		for (i = 0, iMax = NRows_Out[Indd]*NCols; i < iMax; i++)
			Output_Inter[Indd][i] = Input[i];
	}
//array_print_d(NRows_Out[Indd],NCols,Output_Inter[Indd],'C');

	if (d == 1) {
		free(Output_Inter);
		return;
	}

	// s
	Indd = 1;

	if (d == 2)
		Output_Inter[Indd] = Output;
	else
		Output_Inter[Indd] = malloc(NRows_Out[Indd]*NCols * sizeof *Output_Inter[Indd]); // free

	if (NonRedundant[Indd]) {
		sf_swap_d(Output_Inter[Indd-1],NRows_Out[Indd-1],NCols,
		          NIn[2],NOut[0],NIn[1],NOut[0]*NIn[1],1,NOut[0]);
//array_print_d(NRows_Out[Indd-1],NCols,Output_Inter[Indd-1],'C');

		mm_CTN_d(NOut[Indd],NCols*BRows[Indd],NIn[Indd],OP[Indd],Output_Inter[Indd-1],Output_Inter[Indd]);
//array_print_d(NRows_Out[Indd],NCols,Output_Inter[Indd],'C');

		sf_swap_d(Output_Inter[Indd],NRows_Out[Indd],NCols,
		          NIn[2],NOut[1],NOut[0],NOut[1]*NOut[0],1,NOut[1]);
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

	// t
	Indd = 2;

	Output_Inter[Indd] = Output;

	if (NonRedundant[Indd]) {
		sf_swap_d(Output_Inter[Indd-1],NRows_Out[Indd-1],NCols,
		          NOut[0],NOut[1],NIn[2],NOut[1],1,NOut[0]*NOut[1]);
//array_print_d(NRows_Out[Indd-1],NCols,Output_Inter[Indd-1],'C');

		mm_CTN_d(NOut[Indd],NCols*BRows[Indd],NIn[Indd],OP[Indd],Output_Inter[Indd-1],Output_Inter[Indd]);
//array_print_d(NRows_Out[Indd],NCols,Output_Inter[Indd],'C');

		sf_swap_d(Output_Inter[Indd],NRows_Out[Indd],NCols,
		          NOut[2],NOut[1],NOut[0],1,NOut[2]*NOut[0],NOut[2]);
	} else {
		for (i = 0, iMax = NRows_Out[Indd]*NCols; i < iMax; i++)
			Output_Inter[Indd][i] = Output_Inter[Indd-1][i];
	}
	free(Output_Inter[Indd-1]);
//array_print_d(NRows_Out[Indd],NCols,Output_Inter[Indd],'C');

	free(Output_Inter);
}
