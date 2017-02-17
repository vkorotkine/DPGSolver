// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"

#include "element_functions.h"
#include "matrix_functions.h"

/*
 *	Purpose:
 *		Compute weak gradients required for the computation of viscous fluxes.
 *
 *	Comments:
 *		The weak form of the equation is used (i.e. integrated by parts once).
 *
 *	Notation:
 *
 *	References:
 *		Toro(2009)-Riemann_ Solvers_and_Numerical_Methods_for_Fluid_Dynamics
 */

static void explicit_GradW_VOLUME   (void);
static void explicit_GradW_FACE     (void);
static void explicit_GradW_finalize (void);

void explicit_GradW(void)
{
	explicit_GradW_VOLUME();
	explicit_GradW_FACE();
	explicit_GradW_finalize();
}


struct S_OPERATORS {
	unsigned int NvnS, NvnI;
	double       *ChiS_vI, **D_Weak;
};

struct S_Dxyz {
	unsigned int dim, Nbf, Nn;
	double       **D, *C;
};

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


static double *compute_Dxyz(struct S_Dxyz *DxyzInfo, unsigned int d)
{
	/*
	 *	Purpose:
	 *		Compute physical derivative operator matrices using the chain rule.
	 *
	 *	Comments:
	 *		Note the ordering of C specified in the comments of setup_geom_factors.
	 *
	 *	References:
	 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_
	 *		                 Schemes (eq. B.2)
	 */

	unsigned int Nbf, Nn, dim1, IndD, IndC;
	double       **D, *C, *Dxyz;

	Nbf  = DxyzInfo->Nbf;
	Nn   = DxyzInfo->Nn;
	dim1 = DxyzInfo->dim;

	D    = DxyzInfo->D;
	C    = DxyzInfo->C;

	Dxyz = calloc(Nbf*Nn, sizeof *Dxyz); // keep
	for (size_t dim2 = 0; dim2 < d; dim2++) {
		IndD = 0;
		IndC = (dim1+dim2*d)*Nn;
		for (size_t i = 0; i < Nbf; i++) {
			for (size_t j = 0; j < Nn; j++)
				Dxyz[IndD+j] += D[dim2][IndD+j]*C[IndC+j];
			IndD += Nn;
		}
	}

	return Dxyz;
}

static void explicit_GradW_VOLUME(void)
{
	/*
	 *	Purpose:
	 *		Compute intermediate VOLUME contribution to Qhat.
	 *
	 *	Comments:
	 *		This is an intermediate contribution because the multiplication by MInv is not included.
	 *		It is currently hard-coded that GradW is of the same order as the solution.
	 *		Note, if Collocation is enable, that D_Weak includes the inverse cubature weights.
	 */

	// Initialize DB Parameters
	unsigned int d          = DB.d,
	             Neq        = DB.Neq,
	             Collocated = DB.Collocated;

	// Standard datatypes
	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	struct S_Dxyz *DxyzInfo;

	OPS      = malloc(sizeof *OPS);      // free
	DxyzInfo = malloc(sizeof *DxyzInfo); // free

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int NvnS, NvnI;
		double       *ChiS_vI, **DxyzChiS;

		init_ops(OPS,VOLUME);

		NvnS = OPS->NvnS;
		NvnI = OPS->NvnI;

		DxyzInfo->Nbf = OPS->NvnS;
		DxyzInfo->Nn  = OPS->NvnI;
		DxyzInfo->D   = OPS->D_Weak;
		DxyzInfo->C   = VOLUME->C_vI;

		ChiS_vI = OPS->ChiS_vI;

		DxyzChiS = VOLUME->DxyzChiS;
		for (size_t dim1 = 0; dim1 < d; dim1++) {
			double *Dxyz;

			DxyzInfo->dim = dim1;
			Dxyz = compute_Dxyz(DxyzInfo,d); // free/keep

			if (Collocated) {
				DxyzChiS[dim1] = Dxyz;
			} else {
				DxyzChiS[dim1] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,Dxyz,ChiS_vI); // keep
				free(Dxyz);
			}
// Might not need to store DxyzChiS for explicit (ToBeDeleted)

			// Compute intermediate Qhat
			VOLUME->Qhat[dim1] = mm_Alloc_d(CBCM,CBT,CBNT,NvnS,Neq,NvnS,1.0,DxyzChiS[dim1],VOLUME->What); // keep
		}
	}

	free(OPS);
	free(DxyzInfo);
}

static void explicit_GradW_FACE(void)
{
	/*
	 *	Purpose:
	 *		Compute intermediate FACE contribution to Qhat.
	 *
	 *	Comments:
	 *		This is an intermediate contribution because the multiplication by MInv is not included.
	 *		It is currently hard-coded that GradW is of the same order as the solution and that a central numerical flux
	 *		is used.
	 *		Note, if Collocation is enable, that I_Weak includes the inverse cubature weights.
	 */
}

static void explicit_GradW_finalize(void)
{
}
