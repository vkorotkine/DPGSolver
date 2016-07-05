// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "functions.h"
#include "parameters.h"

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

void get_sf_parameters(const unsigned int NIn0, const unsigned int NOut0, double *OP0,
                       const unsigned int NIn1, const unsigned int NOut1, double *OP1,
                       unsigned int NIn_SF[3], unsigned int NOut_SF[3], double *OP_SF[3],
                       const unsigned int d, const unsigned int dim1, const unsigned int Eclass)
{
	/*
	 *	Purpose:
	 *		Set up (s)um (f)actorization parameters.
	 *
	 *	Comments:
	 *		Clearly the input parameters are not very elegant, however, this results in the least overhead when calling
	 *		this function (as compared to using a structure for example).
	 */

	unsigned int dim;

	if (Eclass == C_TP) {
		for (dim = 0; dim < 3; dim++) {
			if (dim == dim1) {
				NIn_SF[dim]  = NIn1;
				NOut_SF[dim] = NOut1;
				OP_SF[dim]   = OP1;
			} else if (dim < d) {
				NIn_SF[dim]  = NIn0;
				NOut_SF[dim] = NOut0;
				OP_SF[dim]   = OP0;
			} else {
				NIn_SF[dim]  = 1;
				NOut_SF[dim] = 1;
				OP_SF[dim]   = NULL;
			}
		}
	} else if (Eclass == C_WEDGE) {
		for (dim = 0; dim < 3; dim++) {
			if (dim == 0) {
				NIn_SF[dim]  = NIn0;
				NOut_SF[dim] = NOut0;
				OP_SF[dim]   = OP0;
			} else if (dim == 2) {
				NIn_SF[dim]  = NIn1;
				NOut_SF[dim] = NOut1;
				OP_SF[dim]   = OP1;
			} else {
				NIn_SF[dim]  = 1;
				NOut_SF[dim] = 1;
				OP_SF[dim]   = NULL;
			}
		}
	} else {
		printf("Error: Unsupported Eclass in get_sf_parameters.\n"), exit(1);
	}
}

void get_sf_parametersF(const unsigned int NIn0, const unsigned int NOut0, double **OP0,
                        const unsigned int NIn1, const unsigned int NOut1, double **OP1,
                        unsigned int NIn_SF[3], unsigned int NOut_SF[3], double *OP_SF[3],
                        const unsigned int d, const unsigned int Vf, const unsigned int Eclass)
{
	/*
	 *	Purpose:
	 *		Set up (s)um (f)actorization parameters for (F)ACET operators.
	 *
	 *	Comments:
	 *		While get_sf_parameters is sufficient for QUADs, HEX ELEMENTs may require three different FACET operators if
	 *		h-refinement is employed.
	 */

	unsigned int f, fh, dimF, dimV1, dimV2;

	f  = Vf/NFREFMAX;
	fh = Vf % NFREFMAX;

	dimF = f/2;

	switch (d) {
	default: // 3D
		switch (Eclass) {
		default: // C_TP
			// FACET term (standard treatment)
			NIn_SF[dimF]  = NIn1;
			NOut_SF[dimF] = NOut1;
			OP_SF[dimF]   = OP1[(f%2)*NFREFMAX];

			// VOLUME term NIn/NOut (standard treatment)
			dimV1 = (dimF+1)%3;
			dimV2 = (dimV1+1)%3;

			NIn_SF[dimV1]  = NIn0;
			NOut_SF[dimV1] = NOut0;
			NIn_SF[dimV2]  = NIn0;
			NOut_SF[dimV2] = NOut0;

			// VOLUME term OP_SF
			switch (fh) {
			case 0: OP_SF[dimV1] = OP0[0]; OP_SF[dimV2] = OP0[0]; break; // Conforming
			case 1: OP_SF[dimV1] = OP0[1]; OP_SF[dimV2] = OP0[0]; break;
			case 2: OP_SF[dimV1] = OP0[2]; OP_SF[dimV2] = OP0[0]; break;
			case 3: OP_SF[dimV1] = OP0[0]; OP_SF[dimV2] = OP0[1]; break;
			case 4: OP_SF[dimV1] = OP0[0]; OP_SF[dimV2] = OP0[2]; break;
			case 5: OP_SF[dimV1] = OP0[1]; OP_SF[dimV2] = OP0[1]; break;
			case 6: OP_SF[dimV1] = OP0[2]; OP_SF[dimV2] = OP0[1]; break;
			case 7: OP_SF[dimV1] = OP0[1]; OP_SF[dimV2] = OP0[2]; break;
			case 8: OP_SF[dimV1] = OP0[2]; OP_SF[dimV2] = OP0[2]; break;
			default:
				printf("Error: Unsupported fh in get_sf_parametersF.\n"), exit(1);
				break;
			}
			break;
		case C_WEDGE:
			// OP_SF[1] (standard)
			NIn_SF[1]  = 1;
			NOut_SF[1] = 1;
			OP_SF[1]   = NULL;

			// dim = 0 and dim = 2
			NIn_SF[0]  = NIn0;
			NOut_SF[0] = NOut0;
			NIn_SF[2]  = NIn1;
			NOut_SF[2] = NOut1;

			if (f < 3) { // QUAD FACETs
				switch (fh) {
				case 0: OP_SF[0] = OP0[f*NFREFMAX+0]; OP_SF[2] = OP1[0]; break;
				case 1: OP_SF[0] = OP0[f*NFREFMAX+1]; OP_SF[2] = OP1[0]; break;
				case 2: OP_SF[0] = OP0[f*NFREFMAX+2]; OP_SF[2] = OP1[0]; break;
				case 3: OP_SF[0] = OP0[f*NFREFMAX+0]; OP_SF[2] = OP1[1]; break;
				case 4: OP_SF[0] = OP0[f*NFREFMAX+0]; OP_SF[2] = OP1[2]; break;
				case 5: OP_SF[0] = OP0[f*NFREFMAX+1]; OP_SF[2] = OP1[1]; break;
				case 6: OP_SF[0] = OP0[f*NFREFMAX+2]; OP_SF[2] = OP1[1]; break;
				case 7: OP_SF[0] = OP0[f*NFREFMAX+1]; OP_SF[2] = OP1[2]; break;
				case 8: OP_SF[0] = OP0[f*NFREFMAX+2]; OP_SF[2] = OP1[2]; break;
				default:
					printf("Error: Unsupported fh for f < 3 in get_sf_parametersF (C_WEDGE).\n"), exit(1);
					break;
				}
			} else if (f < 5) { // TRI FACETs
				if (fh > 4)
					printf("Error: Unsupported fh for f < 5 in get_sf_parametersF (C_WEDGE).\n"), exit(1);

				OP_SF[0] = OP0[fh];
				OP_SF[2] = OP1[((f+1)%2)*NFREFMAX];
			} else {
				printf("Error: Unsupported f in get_sf_parametersF (C_WEDGE).\n"), exit(1);
			}
			break;
		}
		break;
	case 2:
		get_sf_parameters(NIn0,NOut0,OP0[fh],NIn1,NOut1,OP1[(f%2)*NFREFMAX],NIn_SF,NOut_SF,OP_SF,d,dimF,C_TP);
		break;
	}
}

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

double *sf_assemble_d(const unsigned int NIn[3], const unsigned int NOut[3], const unsigned int d, double *BOP[3])
{
	/*
	 *	Purpose:
	 *		Assemble ST(andard) OP(erators) for TP elements using the lower dimensional operators.
	 *
	 *	Comments:
	 *		Standard operators are assembled using sparse BLAS multiplications of lower dimensional operators.
	 *		Sparse matrices are stored in compressed sparse row (CSR) format.
	 *
	 *	Notation:
	 *		BOP : (B)lock (OP)erator in each TP direction.
	 *
	 *	References:
	 *		Intel MKL Sparse BLAS CSR Matrix Storage Format: https://software.intel.com/en-us/node/599835
	 */

	char         transa, matdescra[6];
	unsigned int dim, i, j, k, Bcol, iMax, jMax, kMax, BcolMax, Gcol,
	             Indd, IndG, IndGrow;
	double *OP_ST;

	unsigned int NNZ_BOP[3] = {NIn[0]*NOut[0], NIn[1]*NOut[1], NIn[2]*NOut[2]};
	MKL_INT      BRows[3] = {NIn[1]*NIn[2], NOut[0]*NIn[2], NOut[0]*NOut[1]},
	             dims_OP_ST[2] = {NOut[0]*NOut[1]*NOut[2], NIn[0]*NIn[1]*NIn[2]},
	             dims_DOPr[2]  = {NOut[0]*NIn[1]*NIn[2],   NIn[0]*NIn[1]*NIn[2]},
	             dims_DOPs[2]  = {NOut[0]*NOut[1]*NIn[2],  NOut[0]*NIn[1]*NIn[2]},
	             dims_DOPt[2]  = {NOut[0]*NOut[1]*NOut[2], NOut[0]*NOut[1]*NIn[2]};

	MKL_INT      OPr_rowIndex[BRows[0]*NOut[0]+1], OPs_rowIndex[BRows[1]*NOut[1]+1], OPt_rowIndex[BRows[2]*NOut[2]+1],
	             OPr_cols[BRows[0]*NNZ_BOP[0]],    OPs_cols[BRows[1]*NNZ_BOP[1]],    OPt_cols[BRows[2]*NNZ_BOP[2]];
	double       OPr_vals[BRows[0]*NNZ_BOP[0]],    OPs_vals[BRows[1]*NNZ_BOP[1]],    OPt_vals[BRows[2]*NNZ_BOP[2]],
	             OPr_ST[dims_DOPr[0]*dims_DOPr[1]], OPInter_ST[dims_DOPs[0]*dims_DOPr[1]],
	             alpha, beta, one_d[1] = {1.0};

//	OP_ST = malloc(dims_OP_ST[0]*dims_OP_ST[1] * sizeof *OP_ST); // keep (requires external free)
	OP_ST = calloc(dims_OP_ST[0]*dims_OP_ST[1] , sizeof *OP_ST); // keep (requires external free)

	if (d == 1)
		printf("Error: d must be greater than 1 in sf_assemble_d.\n"), exit(1);

	for (dim = 0; dim < 3; dim++) {
		if (BOP[dim] == NULL)
			BOP[dim] = one_d;
	}

	// r
	Indd = 0;

	IndG = 0; IndGrow = 0;
	for (k = 0, kMax = NIn[2];  k < kMax; k++) {
	for (j = 0, jMax = NIn[1];  j < jMax; j++) {
	for (i = 0, iMax = NOut[0]; i < iMax; i++) {
		OPr_rowIndex[IndGrow] = IndG;
		IndGrow++;

		for (Bcol = 0, BcolMax = NIn[Indd]; Bcol < BcolMax; Bcol++) {
			Gcol = Bcol + j*BcolMax + k*BcolMax*jMax;

			OPr_cols[IndG] = Gcol;
			OPr_vals[IndG] = BOP[Indd][i*BcolMax+Bcol];

			IndG++;
		}
	}}}
	OPr_rowIndex[IndGrow] = BRows[Indd]*NNZ_BOP[Indd];

	// s
	Indd = 1;

	IndG = 0; IndGrow = 0;
	for (k = 0, kMax = NIn[2];  k < kMax; k++) {
	for (j = 0, jMax = NOut[1]; j < jMax; j++) {
	for (i = 0, iMax = NOut[0]; i < iMax; i++) {
		OPs_rowIndex[IndGrow] = IndG;
		IndGrow++;

		for (Bcol = 0, BcolMax = NIn[Indd]; Bcol < BcolMax; Bcol++) {
			Gcol = i + Bcol*iMax + k*iMax*BcolMax;

			OPs_cols[IndG] = Gcol;
			OPs_vals[IndG] = BOP[Indd][j*BcolMax+Bcol];

			IndG++;
		}
	}}}
	OPs_rowIndex[IndGrow] = BRows[Indd]*NNZ_BOP[Indd];

/*
array_print_i(1,BRows[0]*NOut[0]+1,OPr_rowIndex,'R');
array_print_i(1,BRows[0]*NNZ_BOP[0],OPr_cols,'R');
array_print_d(1,BRows[0]*NNZ_BOP[0],OPr_vals,'R');
*/

	MKL_INT job[8] = {1,0,0,2,BRows[0]*NNZ_BOP[0],1}, info = 0;
	mkl_ddnscsr(job,&dims_DOPr[0],&dims_DOPr[1],OPr_ST,&dims_DOPr[1],OPr_vals,OPr_cols,OPr_rowIndex,&info);
	if (info != 0)
		printf("Problem converting from sparse to dense in mkl_ddnscsr.\n"), exit(1);

	alpha = 1.0;
	beta  = 0.0;
	transa = 'N';
	matdescra[0] = 'G';
	matdescra[1] = 'L'; // not used
	matdescra[2] = 'N'; // not used
	matdescra[3] = 'C';

	if (d == 2) {
		mkl_dcsrmm(&transa,&dims_DOPs[0],&dims_DOPr[1],&dims_DOPs[1],&alpha,matdescra,OPs_vals,OPs_cols,
		           OPs_rowIndex,&OPs_rowIndex[1],OPr_ST,&dims_DOPr[1],&beta,OP_ST,&dims_DOPr[1]);
	} else {
		mkl_dcsrmm(&transa,&dims_DOPs[0],&dims_DOPr[1],&dims_DOPs[1],&alpha,matdescra,OPs_vals,OPs_cols,
		           OPs_rowIndex,&OPs_rowIndex[1],OPr_ST,&dims_DOPr[1],&beta,OPInter_ST,&dims_DOPr[1]);

		Indd = 2;
		IndG = 0; IndGrow = 0;
		for (k = 0, kMax = NOut[2]; k < kMax; k++) {
		for (j = 0, jMax = NOut[1]; j < jMax; j++) {
		for (i = 0, iMax = NOut[0]; i < iMax; i++) {
			OPt_rowIndex[IndGrow] = IndG;
			IndGrow++;

			for (Bcol = 0, BcolMax = NIn[Indd]; Bcol < BcolMax; Bcol++) {
				Gcol = i + j*iMax + Bcol*iMax*jMax;

				OPt_cols[IndG] = Gcol;
				OPt_vals[IndG] = BOP[Indd][k*BcolMax+Bcol];

				IndG++;
			}
		}}}
		OPt_rowIndex[IndGrow] = BRows[Indd]*NNZ_BOP[Indd];

		mkl_dcsrmm(&transa,&dims_DOPt[0],&dims_DOPr[1],&dims_DOPt[1],&alpha,matdescra,OPt_vals,OPt_cols,
		           OPt_rowIndex,&OPt_rowIndex[1],OPInter_ST,&dims_DOPr[1],&beta,OP_ST,&dims_DOPr[1]);
	}

//array_print_d(dims_OP_ST[0],dims_OP_ST[1],OP_ST,'R');
//exit(1);

	return OP_ST;
}
