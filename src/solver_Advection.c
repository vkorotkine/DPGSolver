// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_Advection.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "solver_functions.h"
#include "solver_Advection_functions.h"

#include "output_to_paraview.h"
#include "test_code_output_to_paraview.h"

#include "matrix_functions.h"
/*
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "S_FACE.h"
#include "Test.h"

#include "update_VOLUMEs.h"

#include "finalize_LHS.h"
#include "solver_implicit.h"

#include "array_print.h"
*/
/*
 *	Purpose:
 *		Provide functions for the Advection solver.
 *
 *	Comments:
 *		Many of the RHS terms computed are 0; they are included as they are used to check the linearization. Further,
 *		the computational cost is dominated by the global system solve making this additional cost negligible.
 *
 *		It is currently assumed that div (dot) b = 0.
 *
 *	Notation:
 *
 *	References:
 */

static void compute_implicit_VOLUME(void)
{
// Get this function to used implicit_VOLUME_info (ToBeDeleted)
// Will need to add XYZ to FLUXDATA and compute f = b*u
	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *const VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	struct S_Dxyz *const DxyzInfo = malloc(sizeof *DxyzInfo); // free

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		init_VDATA(VDATA,VOLUME);

		unsigned int const NvnS = VDATA->OPS[0]->NvnS,
		                   NvnI = VDATA->OPS[0]->NvnI;

		// compute velocity field
		double *const XYZ_vI = malloc(NvnI*d * sizeof *XYZ_vI); // free
		mm_CTN_d(NvnI,d,VOLUME->NvnG,VDATA->OPS[0]->I_vG_vI,VOLUME->XYZ,XYZ_vI);

		double const *const b_vI = compute_b_Advection(NvnI,XYZ_vI); // free
		free(XYZ_vI);

		DxyzInfo->Nbf = VDATA->OPS[0]->NvnS;
		DxyzInfo->Nn  = VDATA->OPS[0]->NvnI;
		DxyzInfo->D   = (double const *const *const) VDATA->OPS[0]->D_Weak;
		DxyzInfo->C   = VOLUME->C_vI;

		double const *const ChiS_vI = VDATA->OPS[0]->ChiS_vI;

		double *const Dxyz_b = calloc(NvnS*NvnI , sizeof *Dxyz_b); // free
		for (size_t dim = 0; dim < d; dim++) {
			DxyzInfo->dim = dim;
			double *const Dxyz = compute_Dxyz(DxyzInfo,d); // free

			// Add in contribution from velocity field
			mm_diag_d(NvnS,NvnI,&b_vI[dim*NvnI],Dxyz,Dxyz_b,1.0,1.0,'R','R');
			free(Dxyz);
		}
		free((double *) b_vI);

		memset(VOLUME->LHS,0.0,NvnS*NvnS*Neq*Nvar * sizeof *(VOLUME->LHS));
		if (strstr(DB.Form,"Weak")) {
			mm_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,-1.0,0.0,Dxyz_b,ChiS_vI,VOLUME->LHS);
		} else if (strstr(DB.Form,"Strong")) {
			mm_d(CBRM,CBT,CBT,NvnS,NvnS,NvnI,1.0,0.0,ChiS_vI,Dxyz_b,VOLUME->LHS);
		} else {
			EXIT_UNSUPPORTED;
		}
		free(Dxyz_b);

		mm_d(CBCM,CBT,CBNT,NvnS,Nvar,NvnS,1.0,0.0,VOLUME->LHS,VOLUME->What,VOLUME->RHS);
	}

	free(VDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
	free(DxyzInfo);
}

static void compute_implicit_FACE(void)
{
	struct S_OPERATORS_F *OPSL[2], *OPSR[2];
	struct S_FDATA       *const FDATAL = malloc(sizeof *FDATAL), // free
	                     *const FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NumericalFlux *const NFluxData = malloc(sizeof *NFluxData); // free
	FDATAL->NFluxData = NFluxData;
	FDATAR->NFluxData = NFluxData;

	struct S_DATA *const DATA = malloc(sizeof *DATA); // free
	DATA->FDATAL    = FDATAL;
	DATA->FDATAR    = FDATAR;
	DATA->NFluxData = NFluxData;
	DATA->feature   = 'F';
	DATA->imex_type = 'I';

	for (size_t i = 0; i < 2; i++) {
		OPSL[i]  = malloc(sizeof *OPSL[i]);  // free
		OPSR[i]  = malloc(sizeof *OPSR[i]);  // free
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		init_FDATA(FDATAL,FACE,'L');
		init_FDATA(FDATAR,FACE,'R');

	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}
	free(NFluxData);
	free(FDATAL);
	free(FDATAR);
	free(DATA);
}

void solver_Advection(bool const PrintEnabled)
{
	if (strstr(DB.SolverType,"Explicit")) {
printf("%d\n",PrintEnabled);
		EXIT_UNSUPPORTED;
	} else if (strstr(DB.SolverType,"Implicit")) {
		compute_implicit_VOLUME();
		compute_implicit_FACE();

		char *const fNameOut = get_fNameOut("SolFinal_");
		output_to_paraview(fNameOut);
		free(fNameOut);

		EXIT_UNSUPPORTED;
	}
}
