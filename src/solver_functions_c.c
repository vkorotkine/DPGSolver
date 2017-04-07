// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_functions_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "S_OpCSR.h"

#include "solver_functions.h"
#include "matrix_functions.h"
#include "sum_factorization.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Provide solver related functions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void coef_to_values_vI_c(struct S_VDATA *VDATA, const char coef_type)
{
	/*
	 *	Purpose:
	 *		Identical to coef_to_values_vI using complex variables (for complex step verification).
	 */

	if (DB.Collocated) {
		EXIT_UNSUPPORTED; // Set A_vI = VOLUME->Ahat in calling function.
	} else {
		if (!(coef_type == 'W' || coef_type == 'Q'))
			EXIT_UNSUPPORTED;

		unsigned int d    = DB.d,
		             Nvar = DB.Nvar;

		struct S_OPERATORS_V **OPS   = VDATA->OPS;
		struct S_VOLUME      *VOLUME = VDATA->VOLUME;

		if (coef_type == 'W') {
			mm_dcc(CBCM,CBT,CBNT,OPS[0]->NvnI,Nvar,OPS[0]->NvnS,1.0,0.0,OPS[0]->ChiS_vI,VOLUME->What_c,VDATA->W_vI_c);
		} else if (coef_type == 'Q') {
			for (size_t dim = 0; dim < d; dim++)
				mm_dcc(CBCM,CBT,CBNT,OPS[0]->NvnI,Nvar,OPS[0]->NvnS,1.0,0.0,OPS[0]->ChiS_vI,VOLUME->Qhat_c[dim],VDATA->Q_vI_c[dim]);
		}
	}
}

void convert_between_rp_c(const unsigned int Nn, const unsigned int Nrc, const double *C, double complex *Ap,
                          double complex *Ar, const char *conv_type)
{
	/*
	 *	Purpose:
	 *		Identical to convert_between_rp using complex variables (for complex step verification).
	 */

	unsigned int d = DB.d;

	if (strstr(conv_type,"FluxToRef")) {
		memset(Ar,0.0,Nn*Nrc*d * sizeof *Ar);
		unsigned int IndAr, IndAp, IndC;

		for (size_t col = 0; col < Nrc; col++) {
		for (size_t dim1 = 0; dim1 < d; dim1++) {
			IndAr = (Nrc*dim1+col)*Nn;
			for (size_t dim2 = 0; dim2 < d; dim2++) {
				IndAp = (col*d+dim2)*Nn;
				IndC  = (dim1*d+dim2)*Nn;
				for (size_t n = 0; n < Nn; n++)
					Ar[IndAr+n] += Ap[IndAp+n]*C[IndC+n];
			}
		}}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_VOLUME_Inviscid_Weak_c(const unsigned int Nrc, const double complex *Ar_vI, double complex *RLHS,
                                     const char imex_type, struct S_VDATA *VDATA)
{
	/*
	 *	Purpose:
	 *		Identical to finalize_VOLUME_Inviscid_Weak_c using complex variables (for complex step verification).
	 *
	 *	Comments:
	 *		Only supports term_type = "RHS" as the complex solver functions are not used to compute linearizations.
	 */

	unsigned int d = DB.d;

	struct S_OPERATORS_V **OPS = VDATA->OPS;

	unsigned int NvnS = OPS[0]->NvnS,
	             NvnI = OPS[0]->NvnI;

	if (imex_type == 'E') {
		memset(RLHS,0.0,NvnS*Nrc * sizeof *RLHS);
		for (size_t dim = 0; dim < d; dim++)
			mm_dcc(CBCM,CBT,CBNT,NvnS,Nrc,NvnI,1.0,1.0,OPS[0]->D_Weak[dim],&Ar_vI[NvnI*Nrc*dim],RLHS);
	} else {
		EXIT_UNSUPPORTED;
	}
}
