// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "memory_constructors.h"

#include <stdlib.h>
#include <limits.h>

#include "Parameters.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

/*
 *	Purpose:
 *		Allocate memory for and initialize new structures.
 *
 *	Comments:
 *		Split this function into memory_constructor for each struct type to minimize dependencies. (ToBeDeleted)
 *		Change all initializations from 0 to UINT_MAX. (ToBeDeleted)
 *		GfS_fIs/c may not be needed based on FACET_info testing. See comments at the start of explicit_FACET_Info.
 *		(ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

struct S_ELEMENT *New_ELEMENT(void)
{
	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             PMax = DB.PMax,
	             NP   = DB.NP;

	// Standard datatypes
	unsigned int P, Pb, Vf, PbMin, PbMax;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = malloc(sizeof *ELEMENT); // free

	// Mesh
	ELEMENT->present = 0;
	ELEMENT->type    = UINT_MAX;
	ELEMENT->Eclass  = UINT_MAX;
	ELEMENT->d       = UINT_MAX;
	ELEMENT->Nve     = UINT_MAX;
	ELEMENT->Nf      = UINT_MAX;
	ELEMENT->Nvref   = UINT_MAX;
	ELEMENT->NvrefSF = UINT_MAX;

	ELEMENT->Nfve    = calloc(NFMAX         , sizeof *(ELEMENT->Nfve));    // free
	ELEMENT->VeCGmsh = calloc(NVEMAX        , sizeof *(ELEMENT->VeCGmsh)); // free
	ELEMENT->VeFcon  = calloc(NFMAX*NFVEMAX , sizeof *(ELEMENT->VeFcon));  // free
	ELEMENT->NrefV   = calloc(NREFVVARMAX   , sizeof *(ELEMENT->NrefV));   // free

	// Operators

	// h-refinement related
	ELEMENT->NEhref  = UINT_MAX;
	ELEMENT->type_h  = calloc(NEHREFMAX      , sizeof *(ELEMENT->type_h));  // free
	ELEMENT->Nfref   = calloc(NFMAX          , sizeof *(ELEMENT->Nfref));   // free
	ELEMENT->NfMixed = calloc(NFMIXEDMAX     , sizeof *(ELEMENT->NfMixed)); // free
	ELEMENT->VeF     = calloc(NFREFMAX*NFMAX , sizeof *(ELEMENT->VeF));     // free
	ELEMENT->Nvve    = calloc(NVREFMAX       , sizeof *(ELEMENT->Nvve));    // free
	ELEMENT->VeV     = calloc(NVREFMAX       , sizeof *(ELEMENT->VeV));     // free

	// Normals
	ELEMENT->nr = calloc(NFMAX*DMAX , sizeof *(ELEMENT->nr)); // free

	// Plotting
	ELEMENT->connectivity  = calloc(NP , sizeof *(ELEMENT->connectivity));  // free
	ELEMENT->connect_types = calloc(NP , sizeof *(ELEMENT->connect_types)); // free
	ELEMENT->NvnP          = calloc(NP , sizeof *(ELEMENT->NvnP));          // free
	ELEMENT->connect_NE    = calloc(NP , sizeof *(ELEMENT->connect_NE));    // free

	// Operators
	ELEMENT->NvnGs = calloc(NP , sizeof *(ELEMENT->NvnGs)); // free
	ELEMENT->NvnGc = calloc(NP , sizeof *(ELEMENT->NvnGc)); // free
	ELEMENT->NvnCs = calloc(NP , sizeof *(ELEMENT->NvnCs)); // free
	ELEMENT->NvnCc = calloc(NP , sizeof *(ELEMENT->NvnCc)); // free
	ELEMENT->NvnIs = calloc(NP , sizeof *(ELEMENT->NvnIs)); // free
	ELEMENT->NvnIc = calloc(NP , sizeof *(ELEMENT->NvnIc)); // free
	ELEMENT->NvnS  = calloc(NP , sizeof *(ELEMENT->NvnS));  // free
	ELEMENT->NfnS  = calloc(NP , sizeof *(ELEMENT->NfnS));  // free
	ELEMENT->NfnIs = calloc(NP , sizeof *(ELEMENT->NfnIs)); // free
	ELEMENT->NfnIc = calloc(NP , sizeof *(ELEMENT->NfnIc)); // free

	ELEMENT->w_vIs = calloc(NP , sizeof *(ELEMENT->w_vIs)); // free
	ELEMENT->w_vIc = calloc(NP , sizeof *(ELEMENT->w_vIc)); // free

	ELEMENT->ChiS_vP    = calloc(NP , sizeof *(ELEMENT->ChiS_vP));    // free
	ELEMENT->ChiS_vS    = calloc(NP , sizeof *(ELEMENT->ChiS_vS));    // free
	ELEMENT->ChiS_vIs   = calloc(NP , sizeof *(ELEMENT->ChiS_vIs));   // free
	ELEMENT->ChiS_vIc   = calloc(NP , sizeof *(ELEMENT->ChiS_vIc));   // free
	ELEMENT->ChiInvS_vS = calloc(NP , sizeof *(ELEMENT->ChiInvS_vS)); // free

	ELEMENT->ICs = calloc(NP , sizeof *(ELEMENT->ICs)); // free
	ELEMENT->ICc = calloc(NP , sizeof *(ELEMENT->ICc)); // free

	ELEMENT->I_vGs_vP  = calloc(NP , sizeof *(ELEMENT->I_vGs_vP));  // free
	ELEMENT->I_vGs_vGs = calloc(NP , sizeof *(ELEMENT->I_vGs_vGs)); // free
	ELEMENT->I_vGs_vGc = calloc(NP , sizeof *(ELEMENT->I_vGs_vGc)); // free
	ELEMENT->I_vGs_vCs = calloc(NP , sizeof *(ELEMENT->I_vGs_vCs)); // free
	ELEMENT->I_vGs_vIs = calloc(NP , sizeof *(ELEMENT->I_vGs_vIs)); // free
	ELEMENT->I_vGs_vIc = calloc(NP , sizeof *(ELEMENT->I_vGs_vIc)); // free
	ELEMENT->I_vGs_vS  = calloc(NP , sizeof *(ELEMENT->I_vGs_vS));  // free
	ELEMENT->I_vGc_vP  = calloc(NP , sizeof *(ELEMENT->I_vGc_vP));  // free
	ELEMENT->I_vGc_vCc = calloc(NP , sizeof *(ELEMENT->I_vGc_vCc)); // free
	ELEMENT->I_vGc_vIs = calloc(NP , sizeof *(ELEMENT->I_vGc_vIs)); // free
	ELEMENT->I_vGc_vIc = calloc(NP , sizeof *(ELEMENT->I_vGc_vIc)); // free
	ELEMENT->I_vGc_vS  = calloc(NP , sizeof *(ELEMENT->I_vGc_vS));  // free
	ELEMENT->I_vCs_vS  = calloc(NP , sizeof *(ELEMENT->I_vCs_vS));  // free
	ELEMENT->I_vCs_vIs = calloc(NP , sizeof *(ELEMENT->I_vCs_vIs)); // free
	ELEMENT->I_vCs_vIc = calloc(NP , sizeof *(ELEMENT->I_vCs_vIc)); // free
	ELEMENT->I_vCc_vS  = calloc(NP , sizeof *(ELEMENT->I_vCc_vS));  // free
	ELEMENT->I_vCc_vIs = calloc(NP , sizeof *(ELEMENT->I_vCc_vIs)); // free
	ELEMENT->I_vCc_vIc = calloc(NP , sizeof *(ELEMENT->I_vCc_vIc)); // free

	ELEMENT->Ihat_vS_vS = calloc(NP , sizeof *(ELEMENT->Ihat_vS_vS)); // free

	ELEMENT->D_vGs_vCs = calloc(NP , sizeof *(ELEMENT->D_vGs_vCs)); // free
	ELEMENT->D_vGs_vIs = calloc(NP , sizeof *(ELEMENT->D_vGs_vIs)); // free
	ELEMENT->D_vCs_vCs = calloc(NP , sizeof *(ELEMENT->D_vCs_vCs)); // free
	ELEMENT->D_vGc_vCc = calloc(NP , sizeof *(ELEMENT->D_vGc_vCc)); // free
	ELEMENT->D_vGc_vIc = calloc(NP , sizeof *(ELEMENT->D_vGc_vIc)); // free
	ELEMENT->D_vCc_vCc = calloc(NP , sizeof *(ELEMENT->D_vCc_vCc)); // free

	ELEMENT->ChiS_fS     = calloc(NP , sizeof *(ELEMENT->ChiS_fS));     // free
	ELEMENT->ChiS_fIs    = calloc(NP , sizeof *(ELEMENT->ChiS_fIs));    // free
	ELEMENT->ChiS_fIc    = calloc(NP , sizeof *(ELEMENT->ChiS_fIc));    // free
	ELEMENT->ChiS_fIs_sp = calloc(NP , sizeof *(ELEMENT->ChiS_fIs_sp)); // free
	ELEMENT->ChiS_fIc_sp = calloc(NP , sizeof *(ELEMENT->ChiS_fIc_sp)); // free

	ELEMENT->GradChiS_fIs = calloc(NP , sizeof *(ELEMENT->GradChiS_fIs)); // free
	ELEMENT->GradChiS_fIc = calloc(NP , sizeof *(ELEMENT->GradChiS_fIc)); // free

	ELEMENT->I_vGs_fS  = calloc(NP , sizeof *(ELEMENT->I_vGs_fS));  // free
	ELEMENT->I_vGs_fIs = calloc(NP , sizeof *(ELEMENT->I_vGs_fIs)); // free
	ELEMENT->I_vGs_fIc = calloc(NP , sizeof *(ELEMENT->I_vGs_fIc)); // free
	ELEMENT->I_vGc_fS  = calloc(NP , sizeof *(ELEMENT->I_vGc_fS));  // free
	ELEMENT->I_vGc_fIs = calloc(NP , sizeof *(ELEMENT->I_vGc_fIs)); // free
	ELEMENT->I_vGc_fIc = calloc(NP , sizeof *(ELEMENT->I_vGc_fIc)); // free
	ELEMENT->I_vCs_fS  = calloc(NP , sizeof *(ELEMENT->I_vCs_fS));  // free
	ELEMENT->I_vCs_fIs = calloc(NP , sizeof *(ELEMENT->I_vCs_fIs)); // free
	ELEMENT->I_vCs_fIc = calloc(NP , sizeof *(ELEMENT->I_vCs_fIc)); // free
	ELEMENT->I_vCc_fS  = calloc(NP , sizeof *(ELEMENT->I_vCc_fS));  // free
	ELEMENT->I_vCc_fIs = calloc(NP , sizeof *(ELEMENT->I_vCc_fIs)); // free
	ELEMENT->I_vCc_fIc = calloc(NP , sizeof *(ELEMENT->I_vCc_fIc)); // free

	ELEMENT->D_vGs_fIs = calloc(NP , sizeof *(ELEMENT->D_vGs_fIs)); // free
	ELEMENT->D_vGs_fIc = calloc(NP , sizeof *(ELEMENT->D_vGs_fIc)); // free
	ELEMENT->D_vGc_fIs = calloc(NP , sizeof *(ELEMENT->D_vGc_fIs)); // free
	ELEMENT->D_vGc_fIc = calloc(NP , sizeof *(ELEMENT->D_vGc_fIc)); // free

	ELEMENT->Is_Weak_VV    = calloc(NP , sizeof *(ELEMENT->Is_Weak_VV));    // free
	ELEMENT->Ic_Weak_VV    = calloc(NP , sizeof *(ELEMENT->Ic_Weak_VV));    // free
	ELEMENT->Ds_Weak_VV    = calloc(NP , sizeof *(ELEMENT->Ds_Weak_VV));    // free
	ELEMENT->Dc_Weak_VV    = calloc(NP , sizeof *(ELEMENT->Dc_Weak_VV));    // free
	ELEMENT->Ds_Weak_VV_sp = calloc(NP , sizeof *(ELEMENT->Ds_Weak_VV_sp)); // free
	ELEMENT->Dc_Weak_VV_sp = calloc(NP , sizeof *(ELEMENT->Dc_Weak_VV_sp)); // free
	ELEMENT->Is_Weak_FF    = calloc(NP , sizeof *(ELEMENT->Is_Weak_FF));    // free
	ELEMENT->Ic_Weak_FF    = calloc(NP , sizeof *(ELEMENT->Ic_Weak_FF));    // free
	ELEMENT->Is_Weak_FF_sp = calloc(NP , sizeof *(ELEMENT->Is_Weak_FF_sp)); // free
	ELEMENT->Ic_Weak_FF_sp = calloc(NP , sizeof *(ELEMENT->Ic_Weak_FF_sp)); // free

	ELEMENT->L2hat_vS_vS = calloc(NP , sizeof *(ELEMENT->L2hat_vS_vS)); // free
	ELEMENT->GfS_fIs    = calloc(NP , sizeof *(ELEMENT->GfS_fIs));    // free
	ELEMENT->GfS_fIc    = calloc(NP , sizeof *(ELEMENT->GfS_fIc));    // free

	ELEMENT->nOrd_fS   = calloc(NP , sizeof *(ELEMENT->nOrd_fS));  // free
	ELEMENT->nOrd_fIs  = calloc(NP , sizeof *(ELEMENT->nOrd_fIs)); // free
	ELEMENT->nOrd_fIc  = calloc(NP , sizeof *(ELEMENT->nOrd_fIc)); // free


	ELEMENT->I_vGs_vP[1]  = calloc(NP , sizeof **(ELEMENT->I_vGs_vP));
	ELEMENT->I_vGs_vGs[1] = calloc(NP , sizeof **(ELEMENT->I_vGs_vGs));
	ELEMENT->I_vGs_vGc[1] = calloc(NP , sizeof **(ELEMENT->I_vGs_vGc));
	ELEMENT->I_vGs_vCs[1] = calloc(NP , sizeof **(ELEMENT->I_vGs_vCs));
	ELEMENT->I_vGs_vIs[1] = calloc(NP , sizeof **(ELEMENT->I_vGs_vIs));
	ELEMENT->I_vGs_vIc[1] = calloc(NP , sizeof **(ELEMENT->I_vGs_vIc));
	ELEMENT->I_vGs_vS[1]  = calloc(NP , sizeof **(ELEMENT->I_vGs_vS));

	ELEMENT->D_vGs_vCs[1] = calloc(NP , sizeof **(ELEMENT->D_vGs_vCs));
	ELEMENT->D_vGs_vIs[1] = calloc(NP , sizeof **(ELEMENT->D_vGs_vIs));

	ELEMENT->I_vGs_fS[1]  = calloc(NP , sizeof **(ELEMENT->I_vGs_fS));
	ELEMENT->I_vGs_fIs[1] = calloc(NP , sizeof **(ELEMENT->I_vGs_fIs));
	ELEMENT->I_vGs_fIc[1] = calloc(NP , sizeof **(ELEMENT->I_vGs_fIc));

	ELEMENT->D_vGs_fIs[1] = calloc(NP , sizeof **(ELEMENT->D_vGs_fIs));
	ELEMENT->D_vGs_fIc[1] = calloc(NP , sizeof **(ELEMENT->D_vGs_fIc));

	for (P = 0; P < NP; P++) {
		ELEMENT->NfnS[P]  = calloc(NESUBCMAX , sizeof **(ELEMENT->NfnS));
		ELEMENT->NfnIs[P] = calloc(NESUBCMAX , sizeof **(ELEMENT->NfnIs));
		ELEMENT->NfnIc[P] = calloc(NESUBCMAX , sizeof **(ELEMENT->NfnIc));

		ELEMENT->ChiS_vP[P]    = calloc(NP , sizeof **(ELEMENT->ChiS_vP));
		ELEMENT->ChiS_vS[P]    = calloc(NP , sizeof **(ELEMENT->ChiS_vS));
		ELEMENT->ChiS_vIs[P]   = calloc(NP , sizeof **(ELEMENT->ChiS_vIs));
		ELEMENT->ChiS_vIc[P]   = calloc(NP , sizeof **(ELEMENT->ChiS_vIc));
		ELEMENT->ChiInvS_vS[P] = calloc(NP , sizeof **(ELEMENT->ChiInvS_vS));

		ELEMENT->ICs[P] = calloc(NP , sizeof **(ELEMENT->ICs));
		ELEMENT->ICc[P] = calloc(NP , sizeof **(ELEMENT->ICc));

		ELEMENT->I_vGc_vP[P]  = calloc(NP , sizeof **(ELEMENT->I_vGc_vP));
		ELEMENT->I_vGc_vCc[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_vCc));
		ELEMENT->I_vGc_vIs[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_vIs));
		ELEMENT->I_vGc_vIc[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_vIc));
		ELEMENT->I_vGc_vS[P]  = calloc(NP , sizeof **(ELEMENT->I_vGc_vS));
		ELEMENT->I_vCs_vS[P]  = calloc(NP , sizeof **(ELEMENT->I_vCs_vS));
		ELEMENT->I_vCs_vIs[P] = calloc(NP , sizeof **(ELEMENT->I_vCs_vIs));
		ELEMENT->I_vCs_vIc[P] = calloc(NP , sizeof **(ELEMENT->I_vCs_vIc));
		ELEMENT->I_vCc_vS[P]  = calloc(NP , sizeof **(ELEMENT->I_vCc_vS));
		ELEMENT->I_vCc_vIs[P] = calloc(NP , sizeof **(ELEMENT->I_vCc_vIs));
		ELEMENT->I_vCc_vIc[P] = calloc(NP , sizeof **(ELEMENT->I_vCc_vIc));

		ELEMENT->Ihat_vS_vS[P] = calloc(NP , sizeof **(ELEMENT->Ihat_vS_vS));

		ELEMENT->D_vGc_vCc[P] = calloc(NP , sizeof **(ELEMENT->D_vGc_vCc));
		ELEMENT->D_vGc_vIc[P] = calloc(NP , sizeof **(ELEMENT->D_vGc_vIc));
		ELEMENT->D_vCs_vCs[P] = calloc(NP , sizeof **(ELEMENT->D_vCs_vCs));
		ELEMENT->D_vCc_vCc[P] = calloc(NP , sizeof **(ELEMENT->D_vCc_vCc));

		ELEMENT->ChiS_fS[P]     = calloc(NP , sizeof **(ELEMENT->ChiS_fS));
		ELEMENT->ChiS_fIs[P]    = calloc(NP , sizeof **(ELEMENT->ChiS_fIs));
		ELEMENT->ChiS_fIc[P]    = calloc(NP , sizeof **(ELEMENT->ChiS_fIc));
		ELEMENT->ChiS_fIs_sp[P] = calloc(NP , sizeof **(ELEMENT->ChiS_fIs_sp));
		ELEMENT->ChiS_fIc_sp[P] = calloc(NP , sizeof **(ELEMENT->ChiS_fIc_sp));

		ELEMENT->GradChiS_fIs[P] = calloc(NP , sizeof **(ELEMENT->GradChiS_fIs));
		ELEMENT->GradChiS_fIc[P] = calloc(NP , sizeof **(ELEMENT->GradChiS_fIc));

		ELEMENT->I_vGc_fS[P]  = calloc(NP , sizeof **(ELEMENT->I_vGc_fS));
		ELEMENT->I_vGc_fIs[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_fIs));
		ELEMENT->I_vGc_fIc[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_fIc));
		ELEMENT->I_vCs_fS[P]  = calloc(NP , sizeof **(ELEMENT->I_vCs_fS));
		ELEMENT->I_vCs_fIs[P] = calloc(NP , sizeof **(ELEMENT->I_vCs_fIs));
		ELEMENT->I_vCs_fIc[P] = calloc(NP , sizeof **(ELEMENT->I_vCs_fIc));
		ELEMENT->I_vCc_fS[P]  = calloc(NP , sizeof **(ELEMENT->I_vCc_fS));
		ELEMENT->I_vCc_fIs[P] = calloc(NP , sizeof **(ELEMENT->I_vCc_fIs));
		ELEMENT->I_vCc_fIc[P] = calloc(NP , sizeof **(ELEMENT->I_vCc_fIc));

		ELEMENT->D_vGc_fIs[P] = calloc(NP , sizeof **(ELEMENT->D_vGc_fIs));
		ELEMENT->D_vGc_fIc[P] = calloc(NP , sizeof **(ELEMENT->D_vGc_fIc));

		ELEMENT->Is_Weak_VV[P]    = calloc(NP , sizeof **(ELEMENT->Is_Weak_VV));
		ELEMENT->Ic_Weak_VV[P]    = calloc(NP , sizeof **(ELEMENT->Ic_Weak_VV));
		ELEMENT->Ds_Weak_VV[P]    = calloc(NP , sizeof **(ELEMENT->Ds_Weak_VV));
		ELEMENT->Dc_Weak_VV[P]    = calloc(NP , sizeof **(ELEMENT->Dc_Weak_VV));
		ELEMENT->Ds_Weak_VV_sp[P] = calloc(NP , sizeof **(ELEMENT->Ds_Weak_VV_sp));
		ELEMENT->Dc_Weak_VV_sp[P] = calloc(NP , sizeof **(ELEMENT->Dc_Weak_VV_sp));
		ELEMENT->Is_Weak_FF[P]    = calloc(NP , sizeof **(ELEMENT->Is_Weak_FF));
		ELEMENT->Ic_Weak_FF[P]    = calloc(NP , sizeof **(ELEMENT->Ic_Weak_FF));
		ELEMENT->Is_Weak_FF_sp[P] = calloc(NP , sizeof **(ELEMENT->Is_Weak_FF_sp));
		ELEMENT->Ic_Weak_FF_sp[P] = calloc(NP , sizeof **(ELEMENT->Ic_Weak_FF_sp));

		ELEMENT->L2hat_vS_vS[P] = calloc(NP , sizeof **(ELEMENT->L2hat_vS_vS));
		ELEMENT->GfS_fIs[P]    = calloc(NP , sizeof **(ELEMENT->GfS_fIs));
		ELEMENT->GfS_fIc[P]    = calloc(NP , sizeof **(ELEMENT->GfS_fIc));

		if      (P == 0)    PbMin = P,   PbMax = P+1;
		else if (P == PMax) PbMin = P-1, PbMax = PMax;
		else                PbMin = P-1, PbMax = P+1;
		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			ELEMENT->ChiS_vS[P][Pb]    = calloc(NVREFSFMAX , sizeof ***(ELEMENT->ChiS_vS));
			ELEMENT->ChiS_vIs[P][Pb]   = calloc(NVREFSFMAX , sizeof ***(ELEMENT->ChiS_vIs));
			ELEMENT->ChiS_vIc[P][Pb]   = calloc(NVREFMAX   , sizeof ***(ELEMENT->ChiS_vIc));

			ELEMENT->I_vGc_vS[P][Pb]  = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGc_vS));
			ELEMENT->I_vGc_vIs[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGc_vIs));
			ELEMENT->I_vGc_vIc[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGc_vIc));
			ELEMENT->I_vCs_vS[P][Pb]  = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCs_vS));
			ELEMENT->I_vCs_vIs[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCs_vIs));
			ELEMENT->I_vCs_vIc[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCs_vIc));
			ELEMENT->I_vCc_vS[P][Pb]  = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCc_vS));
			ELEMENT->I_vCc_vIs[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCc_vIs));
			ELEMENT->I_vCc_vIc[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCc_vIc));

			ELEMENT->Ihat_vS_vS[P][Pb] = calloc(NVREFMAX , sizeof ***(ELEMENT->Ihat_vS_vS));

			ELEMENT->ChiS_vP[P][Pb]    = calloc(1          , sizeof ***(ELEMENT->ChiS_vP));

			ELEMENT->I_vGc_vP[P][Pb]  = calloc(1          , sizeof ***(ELEMENT->I_vGc_vP));

			ELEMENT->Is_Weak_VV[P][Pb]    = calloc(NVREFSFMAX , sizeof ***(ELEMENT->Is_Weak_VV));
			ELEMENT->Ic_Weak_VV[P][Pb]    = calloc(NVREFSFMAX , sizeof ***(ELEMENT->Ic_Weak_VV));
			if (P == Pb) {
				ELEMENT->ChiInvS_vS[P][Pb] = calloc(1          , sizeof ***(ELEMENT->ChiInvS_vS));

				ELEMENT->ICs[P][Pb] = calloc(1 , sizeof ***(ELEMENT->ICs));
				ELEMENT->ICc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->ICc));

				ELEMENT->I_vGs_vP[1][Pb]  = calloc(1          , sizeof ***(ELEMENT->I_vGs_vP));
				ELEMENT->I_vGs_vGs[1][Pb] = calloc(NVREFMAX   , sizeof ***(ELEMENT->I_vGs_vGs));
				ELEMENT->I_vGs_vGc[1][Pb] = calloc(1          , sizeof ***(ELEMENT->I_vGs_vGc));
				ELEMENT->I_vGs_vCs[1][Pb] = calloc(1          , sizeof ***(ELEMENT->I_vGs_vCs));
				ELEMENT->I_vGs_vIs[1][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGs_vIs));
				ELEMENT->I_vGs_vIc[1][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGs_vIc));
				ELEMENT->I_vGs_vS[1][Pb]  = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGs_vS));
				ELEMENT->I_vGc_vCc[P][Pb] = calloc(1          , sizeof ***(ELEMENT->I_vGc_vCc));

				ELEMENT->D_vGs_vCs[1][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vGs_vCs));
				ELEMENT->D_vGs_vIs[1][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vGs_vIs));
				ELEMENT->D_vGc_vCc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vGc_vCc));
				ELEMENT->D_vGc_vIc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vGc_vIc));
				ELEMENT->D_vCs_vCs[P][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vCs_vCs));
				ELEMENT->D_vCc_vCc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vCc_vCc));

				ELEMENT->D_vGs_vCs[1][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vGs_vCs));
				ELEMENT->D_vGs_vIs[1][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vGs_vIs));
				ELEMENT->D_vGc_vCc[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vGc_vCc));
				ELEMENT->D_vGc_vIc[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vGc_vIc));
				ELEMENT->D_vCs_vCs[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vCs_vCs));
				ELEMENT->D_vCc_vCc[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vCc_vCc));

				ELEMENT->Ds_Weak_VV[P][Pb]    = calloc(1 , sizeof ***(ELEMENT->Ds_Weak_VV));
				ELEMENT->Dc_Weak_VV[P][Pb]    = calloc(1 , sizeof ***(ELEMENT->Dc_Weak_VV));
				ELEMENT->Ds_Weak_VV_sp[P][Pb] = calloc(1 , sizeof ***(ELEMENT->Ds_Weak_VV_sp));
				ELEMENT->Dc_Weak_VV_sp[P][Pb] = calloc(1 , sizeof ***(ELEMENT->Dc_Weak_VV_sp));

				ELEMENT->Ds_Weak_VV[P][Pb][0]    = calloc(d , sizeof ****(ELEMENT->Ds_Weak_VV));
				ELEMENT->Dc_Weak_VV[P][Pb][0]    = calloc(d , sizeof ****(ELEMENT->Dc_Weak_VV));
				ELEMENT->Ds_Weak_VV_sp[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->Ds_Weak_VV_sp));
				ELEMENT->Dc_Weak_VV_sp[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->Dc_Weak_VV_sp));

				ELEMENT->I_vGs_fS[1][Pb]  = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGs_fS));
				ELEMENT->I_vGs_fIs[1][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGs_fIs));
				ELEMENT->I_vGs_fIc[1][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGs_fIc));

				ELEMENT->D_vGs_fIs[1][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->D_vGs_fIs));
				ELEMENT->D_vGs_fIc[1][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->D_vGs_fIc));

				for (Vf = 0; Vf < NFREFMAX*NFMAX; Vf++) {
					ELEMENT->D_vGs_fIs[1][Pb][Vf] = calloc(d , sizeof ****(ELEMENT->D_vGs_fIs));
					ELEMENT->D_vGs_fIc[1][Pb][Vf] = calloc(d , sizeof ****(ELEMENT->D_vGs_fIc));
				}
			}

			ELEMENT->ChiS_fS[P][Pb]     = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fS));
			ELEMENT->ChiS_fIs[P][Pb]    = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fIs));
			ELEMENT->ChiS_fIc[P][Pb]    = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fIc));
			ELEMENT->ChiS_fIs_sp[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fIs_sp));
			ELEMENT->ChiS_fIc_sp[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fIc_sp));

			ELEMENT->GradChiS_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->GradChiS_fIs));
			ELEMENT->GradChiS_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->GradChiS_fIc));

			for (Vf = 0; Vf < NFREFMAX*NFMAX; Vf++) {
				ELEMENT->GradChiS_fIs[P][Pb][Vf] = calloc(d , sizeof ***(ELEMENT->GradChiS_fIs));
				ELEMENT->GradChiS_fIc[P][Pb][Vf] = calloc(d , sizeof ***(ELEMENT->GradChiS_fIc));
			}

			ELEMENT->I_vGc_fS[P][Pb]  = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGc_fS));
			ELEMENT->I_vGc_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGc_fIs));
			ELEMENT->I_vGc_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGc_fIc));
			ELEMENT->I_vCs_fS[P][Pb]  = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCs_fS));
			ELEMENT->I_vCs_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCs_fIs));
			ELEMENT->I_vCs_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCs_fIc));
			ELEMENT->I_vCc_fS[P][Pb]  = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCc_fS));
			ELEMENT->I_vCc_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCc_fIs));
			ELEMENT->I_vCc_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCc_fIc));

			ELEMENT->D_vGc_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->D_vGc_fIs));
			ELEMENT->D_vGc_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->D_vGc_fIc));

			for (Vf = 0; Vf < NFREFMAX*NFMAX; Vf++) {
				ELEMENT->D_vGc_fIs[P][Pb][Vf] = calloc(d , sizeof ****(ELEMENT->D_vGc_fIs));
				ELEMENT->D_vGc_fIc[P][Pb][Vf] = calloc(d , sizeof ****(ELEMENT->D_vGc_fIc));
			}

			ELEMENT->Is_Weak_FF[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->Is_Weak_FF));
			ELEMENT->Ic_Weak_FF[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->Ic_Weak_FF));
			ELEMENT->Is_Weak_FF_sp[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->Is_Weak_FF_sp));
			ELEMENT->Ic_Weak_FF_sp[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->Ic_Weak_FF_sp));

			ELEMENT->L2hat_vS_vS[P][Pb] = calloc(NVREFMAX       , sizeof ***(ELEMENT->L2hat_vS_vS));
			ELEMENT->GfS_fIs[P][Pb]    = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->GfS_fIs));
			ELEMENT->GfS_fIc[P][Pb]    = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->GfS_fIc));
		}

		ELEMENT->nOrd_fS[P]  = calloc(NFORDMAX, sizeof **(ELEMENT->nOrd_fS));
		ELEMENT->nOrd_fIs[P] = calloc(NFORDMAX, sizeof **(ELEMENT->nOrd_fIs));
		ELEMENT->nOrd_fIc[P] = calloc(NFORDMAX, sizeof **(ELEMENT->nOrd_fIc));
	}

	ELEMENT->next = NULL;
	ELEMENT->ELEMENTclass  = calloc(NESUBCMAX  , sizeof *(ELEMENT->ELEMENTclass)); // free
	ELEMENT->ELEMENT_FACET = calloc(NFMIXEDMAX , sizeof *(ELEMENT->ELEMENTclass)); // free

	return ELEMENT;
}


struct S_VOLUME *New_VOLUME(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	struct S_VOLUME *VOLUME;
	VOLUME = malloc(sizeof *VOLUME); // free

	// Structures
	VOLUME->indexl = UINT_MAX;
	VOLUME->indexg = UINT_MAX;
	VOLUME->P      = UINT_MAX;
	VOLUME->type   = UINT_MAX;
	VOLUME->Eclass = UINT_MAX;
	VOLUME->update = 0;
	VOLUME->curved = 0;
	VOLUME->level  = 0;
	VOLUME->NsubF  = calloc(NFMAX , sizeof *(VOLUME->NsubF)); // free

	VOLUME->neigh   = calloc(NFMAX*NFREFMAX , sizeof *(VOLUME->neigh));   // free
	VOLUME->neigh_f = calloc(NFMAX*NFREFMAX , sizeof *(VOLUME->neigh_f)); // free

	for (unsigned int i = 0, iMax = NFMAX*NFREFMAX; i < iMax; i++) {
		VOLUME->neigh[i]   = UINT_MAX;
		VOLUME->neigh_f[i] = UINT_MAX;
	}

	VOLUME->XYZ_vC = NULL; // free

	// Geometry
	VOLUME->NvnG  = UINT_MAX;
	VOLUME->XYZ_S = NULL; // free
	VOLUME->XYZ   = NULL; // free

	VOLUME->detJV_vI = NULL; // free
	VOLUME->C_vC     = NULL; // free
	VOLUME->C_vI     = NULL; // free
	VOLUME->C_vf     = calloc(NFMAX , sizeof *(VOLUME->C_vf)); // free

	// Initialization
	VOLUME->NvnS = UINT_MAX;
	VOLUME->What = NULL; // free
	VOLUME->RES  = NULL; // free

	// Solving
	VOLUME->IndA  = UINT_MAX;
	VOLUME->nnz_d = UINT_MAX;
	VOLUME->nnz_o = UINT_MAX;
	VOLUME->RHS   = NULL; // free
	VOLUME->LHS   = NULL; // free
	VOLUME->MInv  = NULL; // free

	// Linearization testing
	VOLUME->What_c = NULL; // free
	VOLUME->RHS_c  = NULL; // free

	VOLUME->uhat_c = NULL; // free
	VOLUME->qhat_c = calloc(d , sizeof *(VOLUME->qhat_c)); // free

	// hp adaptivity
//	VOLUME->minRES = 0.0;
//	VOLUME->maxRES = 0.0;

	VOLUME->refine_current = 0; // ToBeDeleted: Potentially not needed.
	VOLUME->Vadapt         = 0;
	VOLUME->adapt_type     = UINT_MAX;
	VOLUME->PNew           = UINT_MAX;
	VOLUME->hrefine_type   = UINT_MAX;

	// Poisson
	VOLUME->uhat      = NULL; // free
	VOLUME->qhat      = calloc(d , sizeof *(VOLUME->qhat)); // free
	VOLUME->qhat_uhat = calloc(d , sizeof *(VOLUME->qhat)); // free

	VOLUME->DxyzChiS = calloc(d , sizeof *(VOLUME->DxyzChiS)); // free

	// structs
	VOLUME->next    = NULL;
	VOLUME->grpnext = NULL;
	VOLUME->child0  = NULL; // free (in memory_free_children)
	VOLUME->parent  = NULL; // need not be freed
	VOLUME->FACET   = calloc(NFMAX*NSUBFMAX , sizeof *(VOLUME->FACET)); // free

	return VOLUME;
}

struct S_FACET *New_FACET(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	struct S_FACET *FACET;
	FACET = malloc(sizeof *FACET); // free

	// Structures
	FACET->indexg = UINT_MAX;
	FACET->P      = UINT_MAX;
	FACET->type   = UINT_MAX;
	FACET->BC     = UINT_MAX;
	FACET->level  = 0;
	FACET->update = 0;
	FACET->adapt_type = UINT_MAX;

	FACET->VIn   = NULL; // free (in memory_destructor_V)
	FACET->VOut  = NULL; // free (in memory_destructor_V)
	FACET->VfIn  = UINT_MAX;
	FACET->VfOut = UINT_MAX;

	FACET->IndOrdInOut = UINT_MAX;
	FACET->IndOrdOutIn = UINT_MAX;

	// Geometry
	FACET->curved  = UINT_MAX;
	FACET->typeInt = UINT_MAX;

	FACET->XYZ_fI   = NULL; // free
	FACET->XYZ_fS   = NULL; // free
	FACET->n_fI     = NULL; // free
	FACET->n_fS     = NULL; // free
	FACET->detJF_fI = NULL; // free
	FACET->detJF_fS = NULL; // free

	FACET->detJVIn_fI  = NULL; // free
	FACET->detJVOut_fI = NULL; // free

	// Solving
	FACET->RHSIn  = NULL; // free (in finalize_RHS)
	FACET->RHSOut = NULL; // free (in finalize_RHS)

	FACET->LHSInIn   = NULL; // free (in finalize_LHS)
	FACET->LHSOutIn  = NULL; // free (in finalize_LHS)
	FACET->LHSInOut  = NULL; // free (in finalize_LHS)
	FACET->LHSOutOut = NULL; // free (in finalize_LHS)

	// Poisson
	FACET->qhatIn  = calloc(d , sizeof *(FACET->qhatIn));  // free
	FACET->qhatOut = calloc(d , sizeof *(FACET->qhatOut)); // free
	FACET->qhat_uhatInIn   = calloc(d , sizeof *(FACET->qhat_uhatInIn));   // free
	FACET->qhat_uhatOutIn  = calloc(d , sizeof *(FACET->qhat_uhatOutIn));  // free
	FACET->qhat_uhatInOut  = calloc(d , sizeof *(FACET->qhat_uhatInOut));  // free
	FACET->qhat_uhatOutOut = calloc(d , sizeof *(FACET->qhat_uhatOutOut)); // free

	// Linearization testing
	FACET->RHSIn_c  = NULL; // free (in finalize_RHS_c)
	FACET->RHSOut_c = NULL; // free (in finalize_RHS_c)

	FACET->qhatIn_c  = calloc(d , sizeof *(FACET->qhatIn_c));  // free
	FACET->qhatOut_c = calloc(d , sizeof *(FACET->qhatOut_c)); // free

	FACET->next   = NULL;
	FACET->child0 = NULL; // free (in memory_free_children)
	FACET->parent = NULL; // need not be freed

	return FACET;
}
