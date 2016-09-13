// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "solver_poisson.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

/*
#include "Parameters.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "adaptation.h"
#include "update_VOLUMEs.h"
#include "implicit_VOLUME_info.h"
#include "implicit_FACET_info.h"
#include "finalize_LHS.h"
#include "output_to_paraview.h"

#include "Macros.h"
#include "array_print.h"
*/

/*
 *	Purpose:
 *		Perform the implicit solve for the Poisson equation using Petsc's KSP object.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME)
{
	// Standard datatypes
	unsigned int P, type, curved;

	struct S_ELEMENT *ELEMENT;

	P      = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);

	OPS->NvnS = ELEMENT->NvnS[P];
	if (!curved) {
		OPS->NvnI = ELEMENT->NvnIs[P];

		OPS->ChiS_vI = ELEMENT->ChiS_vIs[P][P][0];
		OPS->D_Weak  = ELEMENT->Ds_Weak_VV[P][P][0];
	} else {
		OPS->NvnI = ELEMENT->NvnIc[P];

		OPS->ChiS_vI = ELEMENT->ChiS_vIc[P][P][0];
		OPS->D_Weak  = ELEMENT->Dc_Weak_VV[P][P][0];
	}
}

static void compute_qhat_u_VOLUME(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	double *MInv, **D, *C_vI, **Dxyz, *Sxyz;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		init_ops(OPS,VOLUME);

		NvnI = OPS->NvnI;
		NvnS = OPS->NvnS;

		MInv = VOLUME->MInv;
		C_vI = VOLUME->C_vI;

		// Construct physical derivative operator matrices
		D = OPS->D_Weak;

		Dxyz = malloc(d * sizeof *Dxyz); // free
		for (dim1 = 0; dim1 < d; dim1++)
			Dxyz[dim1] = malloc(NvnS*NvnI * sizeof *Dxyz[dim1]); // free

		for (dim1 = 0; dim1 < d; dim1++) {
		for (dim2 = 0; dim2 < d; dim2++) {
			IndC = (dim1+dim2*d)*NvnI;
			for (i = 0; i < NvnS; i++) {
				for (j = 0; j < NvnI; j++) {
					Dxyz[dim1][IndD+j] = D[dim2][IndD+j]*C_vI[IndC+j];
				}
				IndD += NvnI;
			}
		}}

		// Construct VOLUME q_u term (Removed "-ve" sign)
		ChiS_vI = OPS->ChiS_vI;

		for (dim1 = 0; dim1 < d; dim1++) {
			Sxyz   = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnI,NvnS,1.0,MInv,Dxyz[dim1]); // free
			qhat_u = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,Sxyz,ChiS_vI);    // keep

			free(Sxyz);
			VOLUME->qhat_u[dim1] = qhat_u;
		}

		for (dim1 = 0; dim1 < d; dim1++)
			free(Dxyz[dim1]);
		free(Dxyz);

	}
	free(OPS);
}

static void compute_qhat_u_FACET(void)
{
	// Standard datatypes

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACET     *FACET;

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn = FACET->VIn;

		VOut = FACET->VOut;



		switch (PoissonTraceFluxType) {
		case IP:
			jacobian_tflux_IP(NfnI,1,duNumduIn_fI,'L');
			jacobian_tflux_IP(NfnI,1,duNumduOut_fI,'R');
			break;
		default:
			printf("Error: Unsupported PoissonTraceFluxType.\n"), EXIT_MSG;
			break;
		}

		// RHS
// Include terms related to BCs if relevant. (ToBeDeleted)

		// LHS
		I_FF = OPSIn->I_Weak_FF[VfIn];

		MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NfnI,NvnSIn,1.0,VIn->MInv,I_FF); // free

		// InIn
		for (i = 0; i < NvnSIn; i++) {
		for (j = 0; j < NfnI; j++) {
			MInvI_FF[IndMInvI+j] *= detJF_fI[j];

			for (dim = 0; dim < d; dim++)
				FACET->qhat_uInIn[dim][IndMInvI+j]  = MInvI_FF[IndMInvI+j]*n_fI[j*d+dim]*duNumduIn_fI[j];
		}}

		if (!Boundary) {
			// OutIn
			for (i = 0; i < NvnSIn; i++) {
			for (j = 0; j < NfnI; j++) {
				for (dim = 0; dim < d; dim++)
					FACET->qhat_uOutIn[dim][IndMInvI+j] = MInvI_FF[IndMInvI+j]*n_fI[j*d+dim]*duNumduOut_fI[j];
			}}

			// InOut

		}

		free(MInvI_FF);

	}
}

void solver_poisson(void)
{
	// Initialize DB Parameters
	unsigned int OutputInterval = DB.OutputInterval,
	             Nvar           = DB.Nvar;

	unsigned int PrintTesting = 1;

	// Standard datatypes

	compute_qhat_u_VOLUME();
	compute_qhat_u_FACET();


}
