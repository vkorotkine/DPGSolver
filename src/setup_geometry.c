#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

//#include "petscsys.h"

/*
 *	Purpose:
 *		Set up geometry.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
*/

void sf_operate_d(const int NOut, const int NCols, const int NIn, const int BRowMaxIn, const double *OP,
                  const double *Input, const double *Output)
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
		mm_d(CblasColMajor,CblasTrans,CblasNoTrans,NOut,NCols,NIn,1.0,OP,&Input[IndIn],&Output[IndOut]);

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

			RowInd   = i*iStep+j+k*kStep;
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
			RowInd   = i*iStep+j*jStep1+k;
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

void sf_apply_d(const double *Input, double *Output, const int NIn[3], const int NOut[3], const int NCols,
                double *OP[3], const int Diag[3], const int d)
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
	 *		Make sure that appropriate variables are made declared with 'register' and 'unsigned' (ToBeDeleted).
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



void setup_geometry()
{
	// Initialize DB Parameters
	char *MeshType = DB.MeshType;
	int  ExactGeom = DB.ExactGeom,
	     d         = DB.d,
	     NV        = DB.NV,
		 *NE       = DB.NE,
		 *EToVe    = DB.EToVe,

	     Testing   = DB.Testing;

	double *VeXYZ  = DB.VeXYZ;

	int  PrintTesting = 0, MPIrank = DB.MPIrank;

	// Standard datatypes
	int i, ve, dim, v, P, vn,
	    Nve, Vs, PMax, NvnGs, NvnGc,
		NIn, NOut, NIn_SF[3], NOut_SF[3], NCols, Diag[3], NOut_Total,
		*VeC;
	double *XYZc, *XYZs,
	       *I_vGs_vGc, *Input_SF, *OP_SF[3];

	struct S_ELEMENT *ELEMENT, *ELEMENT_class[2];
	struct S_VOLUME  *VOLUME;

	Vs = 0; for (i = 0; i < d; i++) Vs += NE[i];

	// Modify vertex locations if exact geometry is known
	if (ExactGeom) {
		if(!MPIrank) printf("    Modify vertex nodes if exact geometry is known\n");
		printf("Did not yet verify the implementation.\n");
		vertices_to_exact_geom();
	}

	// Set up global XYZ VOLUME coordinates at (s)tart (i.e. before curving)
	VOLUME = DB.VOLUME;
	v = 0;
	while (VOLUME != NULL) {
		P = VOLUME->P;

		ELEMENT          = get_ELEMENT_type(VOLUME->type);

		VeC   = ELEMENT->VeC;
		NvnGs = ELEMENT->NvnGs[0];

		XYZc = malloc (NvnGs*d * sizeof *XYZc); // keep
		VOLUME->XYZc = XYZc;

		// Note: XYZc may be interpreted as [X Y Z] where each of X, Y, Z are column vectors
		for (ve = 0; ve < NvnGs; ve++) {
		for (dim = 0; dim < d; dim++) {
			XYZc[dim*NvnGs+ve] = VeXYZ[EToVe[(Vs+v)*8+VeC[ve]]*d+dim];
		}}

		if (!VOLUME->curved) {
			// If not curved, the P1 geometry representation sufficies to fully specify the element.
			XYZs = malloc(NvnGs*d * sizeof *XYZs); // keep
			VOLUME->XYZs = XYZs;

			for (dim = 0; dim < d; dim++) {
			for (vn = 0; vn < NvnGs; vn++) {
				XYZs[dim*NvnGs+vn] = XYZc[dim*NvnGs+vn];
			}}
		} else {
			if (VOLUME->Eclass == C_TP) {
				ELEMENT_class[0] = get_ELEMENT_Eclass(VOLUME->Eclass,C_TP);

				NvnGs = ELEMENT_class[0]->NvnGs[0];
				NvnGc = ELEMENT_class[0]->NvnGc[P];

				I_vGs_vGc = ELEMENT_class[0]->I_vGs_vGc[P];

				Input_SF = XYZc; // note multi column input

				NIn = NvnGs;
				for (dim = 0; dim < 3; dim++) {
					if (dim < d) NIn_SF[dim] = NIn;
					else         NIn_SF[dim] = 1;
				}

				NOut = NvnGc;
				NOut_Total = 1;
				for (dim = 0; dim < 3; dim++) {
					if (dim < d) NOut_SF[dim] = NOut;
					else         NOut_SF[dim] = 1;
					NOut_Total *= NOut_SF[dim];
				}

				NCols    = d*1; // d coordinates * 1 element
				OP_SF[0] = I_vGs_vGc;
				OP_SF[1] = OP_SF[0];
				OP_SF[2] = OP_SF[0];
				for (i = 0; i < 3; i++) Diag[i] = 0;

				XYZs = malloc(NOut_Total*NCols * sizeof *XYZs);
				sf_apply_d(Input_SF,XYZs,NIn_SF,NOut_SF,NCols,OP_SF,Diag,d);
			} else if (VOLUME->Eclass == C_SI) {
				NvnGs = ELEMENT->NvnGs[0];
				NvnGc = ELEMENT->NvnGc[P];
				I_vGs_vGc = ELEMENT->I_vGs_vGc[P];

				NCols = d*1; // d coordinates * 1 element
				XYZs = malloc(NvnGc*NCols * sizeof *XYZs);
				mm_d(CblasColMajor,CblasTrans,CblasNoTrans,NvnGc,NCols,NvnGs,1.0,I_vGs_vGc,XYZc,XYZs);
			}
		}
		VOLUME->XYZs = XYZs;

array_print_d(NOut_Total,d,XYZs,'C');

test_speed_mm_d();

exit(1);



		v++;
		VOLUME = VOLUME->next;
	}

	// Find node index ordering on each FACET
	/* WRITE A TEST ROUTINE TO MAKE SURE THAT THIS IS WORKING? (ToBeDeleted)
	 * Note: Ideally, this information can be given without relying on matching of physical points as there will be a
	 *       potential for non-conforming elements.
	 * For the "sum-factorization" on simplex elements, the nodes must be ordered according to their symmetries, which
	 * is unrelated to the vertex positions? Perhaps define the nodes only in one section (1/3) of the reference
	 * triangle with the associated multiplicity.
	 * Note: After having read through Hesthaven(2000) on the sum-factorization on triangles, it seemed intuitive to
	 *       attempt a similar extension for the TP case, using the symmetry about 0. This succeeded, resulting in an
	 *       asymptotic complexity reduction of 2. This does require re-ordering of both the nodes and basis functions
	 *       however => Do this. Also, while I thought I might be onto a new result, a very similar demonstration seems
	 *       to have been made in Solomonoff(1992)-A_Fast_Algorithm_for_Spectral_Differentiation, although he does not
	 *       seem to have considered the interpolation.
	 */



	// Performing analytical mesh curving if MeshType == ToBeCurved
	if (strstr(MeshType,"ToBeCurved") != NULL) {
		printf("    Modify Vertex and VOLUME Nodes of ToBeCurved Mesh\n");
	}


}
