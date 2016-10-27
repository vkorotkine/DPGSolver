// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "memory_destructors.h"

#include <stdlib.h>

#include "Parameters.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

#include "array_free.h"
#include "element_functions.h"
#include "adaptation.h"

/*
 *	Purpose:
 *		Free memory of various structures.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void memory_destructor_E(struct S_ELEMENT *ELEMENT)
{
	// Initialize DB Parameters
	unsigned int d  = DB.d,
	             NP = DB.NP;

	// Mesh
	free(ELEMENT->Nfve);
	free(ELEMENT->VeCGmsh);
	free(ELEMENT->VeFcon);
	free(ELEMENT->NrefV);

	// h-refinement related
	free(ELEMENT->type_h);
	free(ELEMENT->Nfref);
	free(ELEMENT->NfMixed);
	array_free2_d(NFREFMAX*NFMAX,ELEMENT->VeF);
	free(ELEMENT->Nvve);
	array_free2_d(NVREFMAX,ELEMENT->VeV);

	// Normals
	free(ELEMENT->nr);

	// Plotting
	free(ELEMENT->NvnP);
	free(ELEMENT->connect_NE);
	array_free2_ui(NP,ELEMENT->connectivity);
	array_free2_ui(NP,ELEMENT->connect_types);

	// Operators
	free(ELEMENT->NvnGs);
	free(ELEMENT->NvnGc);
	free(ELEMENT->NvnCs);
	free(ELEMENT->NvnCc);
	free(ELEMENT->NvnIs);
	free(ELEMENT->NvnIc);
	free(ELEMENT->NvnS);

	array_free3_d(NP,NESUBCMAX,ELEMENT->w_fIs);
	array_free3_d(NP,NESUBCMAX,ELEMENT->w_fIc);

	array_free2_d(NP,ELEMENT->w_vIs);
	array_free2_d(NP,ELEMENT->w_vIc);

	array_free2_ui(NP,ELEMENT->NfnS);
	array_free2_ui(NP,ELEMENT->NfnIs);
	array_free2_ui(NP,ELEMENT->NfnIc);

	array_free4_d(NP,NP,1,         ELEMENT->ChiS_vP);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->ChiS_vS);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->ChiS_vIs);
	array_free4_d(NP,NP,NVREFMAX,  ELEMENT->ChiS_vIc);
	array_free4_d(NP,NP,1,         ELEMENT->ChiInvS_vS);

	array_free4_d(NP,NP,1,ELEMENT->ICs);
	array_free4_d(NP,NP,1,ELEMENT->ICc);

	array_free5_d(NP,NP,1,d,ELEMENT->GradChiS_vIs);
	array_free5_d(NP,NP,1,d,ELEMENT->GradChiS_vIc);

	array_free4_d(NP,NP,1,         ELEMENT->I_vGs_vP);
	array_free4_d(NP,NP,NVREFMAX,  ELEMENT->I_vGs_vGs);
	array_free4_d(NP,NP,1,         ELEMENT->I_vGs_vGc);
	array_free4_d(NP,NP,1,         ELEMENT->I_vGs_vCs);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vGs_vIs);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vGs_vIc);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vGs_vS);
	array_free4_d(NP,NP,1,         ELEMENT->I_vGc_vP);
	array_free4_d(NP,NP,1,         ELEMENT->I_vGc_vCc);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vGc_vIs);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vGc_vIc);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vGc_vS);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vCs_vS);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vCs_vIs);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vCs_vIc);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vCc_vS);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vCc_vIs);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vCc_vIc);

	array_free4_d(NP,NP,NVREFMAX,ELEMENT->Ihat_vS_vS);

	array_free5_d(NP,NP,1,d,ELEMENT->D_vGs_vCs);
	array_free5_d(NP,NP,1,d,ELEMENT->D_vGs_vIs);
	array_free5_d(NP,NP,1,d,ELEMENT->D_vCs_vCs);
	array_free5_d(NP,NP,1,d,ELEMENT->D_vGc_vCc);
	array_free5_d(NP,NP,1,d,ELEMENT->D_vGc_vIc);
	array_free5_d(NP,NP,1,d,ELEMENT->D_vCc_vCc);

	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->ChiS_fS);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->ChiS_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->ChiS_fIc);
	array_free4_CSR_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->ChiS_fIs_sp);
	array_free4_CSR_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->ChiS_fIc_sp);

	array_free5_d(NP,NP,NFREFMAX*NFMAX,d,ELEMENT->GradChiS_fIs);
	array_free5_d(NP,NP,NFREFMAX*NFMAX,d,ELEMENT->GradChiS_fIc);

	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGs_fS);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGs_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGs_fIc);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGc_fS);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGc_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGc_fIc);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCs_fS);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCs_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCs_fIc);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCc_fS);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCc_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCc_fIc);

	array_free5_d(NP,NP,NFREFMAX*NFMAX,d,ELEMENT->D_vGs_fIs);
	array_free5_d(NP,NP,NFREFMAX*NFMAX,d,ELEMENT->D_vGs_fIc);
	array_free5_d(NP,NP,NFREFMAX*NFMAX,d,ELEMENT->D_vGc_fIs);
	array_free5_d(NP,NP,NFREFMAX*NFMAX,d,ELEMENT->D_vGc_fIc);

	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->Is_Weak_VV);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->Ic_Weak_VV);

	array_free5_d(NP,NP,1,d,ELEMENT->Ds_Weak_VV);
	array_free5_d(NP,NP,1,d,ELEMENT->Dc_Weak_VV);
	array_free5_CSR_d(NP,NP,1,d,ELEMENT->Ds_Weak_VV_sp);
	array_free5_CSR_d(NP,NP,1,d,ELEMENT->Dc_Weak_VV_sp);

	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->Is_Weak_FF);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->Ic_Weak_FF);
	array_free4_CSR_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->Is_Weak_FF_sp);
	array_free4_CSR_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->Ic_Weak_FF_sp);

	array_free4_d(NP,NP,NVREFMAX      ,ELEMENT->L2hat_vS_vS);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->GfS_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->GfS_fIc);

	array_free3_ui(NP,NFORDMAX,ELEMENT->nOrd_fS);
	array_free3_ui(NP,NFORDMAX,ELEMENT->nOrd_fIs);
	array_free3_ui(NP,NFORDMAX,ELEMENT->nOrd_fIc);


	free(ELEMENT->ELEMENTclass);
	free(ELEMENT->ELEMENT_FACET);

	free(ELEMENT);
}

void memory_destructor_L2_projection(const unsigned int EType)
{
	unsigned int NP = DB.NP;

	unsigned int i, P, Pb, PbMin, PbMax;
	struct S_ELEMENT *ELEMENT = get_ELEMENT_type(EType);

	for (P = 0; P < NP; P++) {
		if (ELEMENT->w_vIc[P]) {
			free(ELEMENT->w_vIc[P]);
			ELEMENT->w_vIc[P] = NULL;
		}

		get_Pb_range(P,&PbMin,&PbMax);
		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			for (i = 0; i < NVREFMAX; i++) {
				if (ELEMENT->ChiS_vIc[P][Pb][i]) {
					free(ELEMENT->ChiS_vIc[P][Pb][i]);
					ELEMENT->ChiS_vIc[P][Pb][i] = NULL;
				}
			}
		}
	}
}

void memory_destructor_V(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Structures
	free(VOLUME->XYZ_vC);
	free(VOLUME->NsubF);
	free(VOLUME->neigh);
	free(VOLUME->neigh_f);

	// Geometry
	free(VOLUME->XYZ_S);
	free(VOLUME->XYZ);

	free(VOLUME->detJV_vI);
	free(VOLUME->C_vC);
	free(VOLUME->C_vI);
	array_free2_d(NFMAX,VOLUME->C_vf);

	// Initialization
	free(VOLUME->What);
	free(VOLUME->RES);

	// Solving
	free(VOLUME->RHS);
	free(VOLUME->LHS);
	free(VOLUME->MInv);

	// Linearization testing
	free(VOLUME->What_c);
	free(VOLUME->RHS_c);

	free(VOLUME->uhat_c);
	array_free2_cmplx(d,VOLUME->qhat_c);

	// Poisson
	free(VOLUME->uhat);
	array_free2_d(d,VOLUME->qhat);
	array_free2_d(d,VOLUME->qhat_uhat);
	array_free2_d(d,VOLUME->DxyzChiS);

	// structs
	free(VOLUME->FACET);

	free(VOLUME);
}

void memory_destructor_F(struct S_FACET *FACET)
{
	/*
	 *	Comments:
	 *		LHS** terms are not necessarily freed while testing functions as finalize_LHS need not be called in all
	 *		cases.
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Geometry
	free(FACET->XYZ_fI);
	free(FACET->XYZ_fS);
	free(FACET->n_fI);
	free(FACET->n_fS);
	free(FACET->detJF_fI);
	free(FACET->detJF_fS);
	free(FACET->detJVIn_fI);
	free(FACET->detJVOut_fI);

	// Solving
	if (FACET->LHSInIn)
		free(FACET->LHSInIn);
	if (FACET->LHSOutIn)
		free(FACET->LHSOutIn);
	if (FACET->LHSInOut)
		free(FACET->LHSInOut);
	if (FACET->LHSOutOut)
		free(FACET->LHSOutOut);

	// Poisson
	array_free2_d(d,FACET->qhatIn);
	array_free2_d(d,FACET->qhatOut);
	array_free2_d(d,FACET->qhat_uhatInIn);
	array_free2_d(d,FACET->qhat_uhatOutIn);
	array_free2_d(d,FACET->qhat_uhatInOut);
	array_free2_d(d,FACET->qhat_uhatOutOut);

	// Linearization testing
	free(FACET->RHSIn_c);
	free(FACET->RHSOut_c);

	array_free2_cmplx(d,FACET->qhatIn_c);
	array_free2_cmplx(d,FACET->qhatOut_c);

	free(FACET);
}
