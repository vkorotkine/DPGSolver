// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>

#include "database.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Free memory of various structures.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	Memory freed:
 *
 *	References:
 *
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

	// h-refinement related
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

	array_free2_d(NP,ELEMENT->w_vIs);
	array_free2_d(NP,ELEMENT->w_vIc);

	array_free2_ui(NP,ELEMENT->NfnIs);
	array_free2_ui(NP,ELEMENT->NfnIc);

	array_free4_d(NP,NP,1,         ELEMENT->ChiS_vP);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->ChiS_vIs);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->ChiS_vIc);
	array_free4_d(NP,NP,1,         ELEMENT->ChiInvS_vS);

	array_free4_d(NP,NP,1,ELEMENT->ICs);
	array_free4_d(NP,NP,1,ELEMENT->ICc);

	array_free4_d(NP,NP,1,         ELEMENT->I_vGs_vP);
	array_free4_d(NP,NP,1,         ELEMENT->I_vGs_vGc);
	array_free4_d(NP,NP,1,         ELEMENT->I_vGs_vCs);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vGs_vIs);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vGs_vIc);
	array_free4_d(NP,NP,1,         ELEMENT->I_vGs_vS);
	array_free4_d(NP,NP,1,         ELEMENT->I_vGc_vP);
	array_free4_d(NP,NP,1,         ELEMENT->I_vGc_vCc);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vGc_vIs);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vGc_vIc);
	array_free4_d(NP,NP,1,         ELEMENT->I_vGc_vS);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vCs_vIs);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vCs_vIc);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vCc_vIs);
	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->I_vCc_vIc);

	array_free4_d(NP,NP,NVREFSFMAX,ELEMENT->Ihat_vS_vS);

	array_free5_d(NP,NP,1,d,ELEMENT->D_vGs_vCs);
	array_free5_d(NP,NP,1,d,ELEMENT->D_vGs_vIs);
	array_free5_d(NP,NP,1,d,ELEMENT->D_vCs_vCs);
	array_free5_d(NP,NP,1,d,ELEMENT->D_vGc_vCc);
	array_free5_d(NP,NP,1,d,ELEMENT->D_vGc_vIc);
	array_free5_d(NP,NP,1,d,ELEMENT->D_vCc_vCc);

	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->ChiS_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->ChiS_fIc);
	array_free4_CSR_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->ChiS_fIs_sp);
	array_free4_CSR_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->ChiS_fIc_sp);

	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGs_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGs_fIc);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGc_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGc_fIc);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCs_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCs_fIc);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCc_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCc_fIc);

	array_free4_d(NP,NP,1,ELEMENT->Is_Weak_VV);
	array_free4_d(NP,NP,1,ELEMENT->Ic_Weak_VV);

	array_free5_d(NP,NP,1,d,ELEMENT->Ds_Weak_VV);
	array_free5_d(NP,NP,1,d,ELEMENT->Dc_Weak_VV);
	array_free5_CSR_d(NP,NP,1,d,ELEMENT->Ds_Weak_VV_sp);
	array_free5_CSR_d(NP,NP,1,d,ELEMENT->Dc_Weak_VV_sp);

	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->Is_Weak_FF);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->Ic_Weak_FF);
	array_free4_CSR_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->Is_Weak_FF_sp);
	array_free4_CSR_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->Ic_Weak_FF_sp);

	array_free3_ui(NP,NFORDMAX,ELEMENT->nOrd_fIs);
	array_free3_ui(NP,NFORDMAX,ELEMENT->nOrd_fIc);


	free(ELEMENT->ELEMENTclass);
	free(ELEMENT->ELEMENT_FACET);

	free(ELEMENT);
}

void memory_destructor_V(struct S_VOLUME *VOLUME)
{
//	int NP = DB.NP;

	// Structures
	free(VOLUME->XYZ_vC);

	// Geometry
	free(VOLUME->XYZ_S);
	free(VOLUME->XYZ);

	free(VOLUME->detJV_vI);
	free(VOLUME->C_vI);
	array_free2_d(NFMAX,VOLUME->C_vf);

	// Initialization
	free(VOLUME->What);
	free(VOLUME->RES);

	// Solving
	free(VOLUME->wdetJV_vI);
	free(VOLUME->MInv);

	free(VOLUME);
}

void memory_destructor_F(struct S_FACET *FACET)
{
	free(FACET->XYZ_fI);
	free(FACET->n_fI);
	free(FACET->detJF_fI);

	free(FACET);
}
