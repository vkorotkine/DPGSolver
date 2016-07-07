// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>

#include "database.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Allocate memory for and initialize new structures.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *
 */

struct S_ELEMENT *New_ELEMENT(void)
{
	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             PMax = DB.PMax,
	             NP   = DB.NP;

	// Standard datatypes
	unsigned int P, Pb, PbMin, PbMax;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = malloc(sizeof *ELEMENT); // free

	// Mesh
	ELEMENT->present = 0;
	ELEMENT->type    = 0;
	ELEMENT->Eclass  = 0;
	ELEMENT->d       = 0;
	ELEMENT->Nve     = 0;
	ELEMENT->Nf      = 0;
	ELEMENT->Nvref   = 0;
	ELEMENT->NvrefSF = 0;

	ELEMENT->Nfve    = calloc(NFMAX         , sizeof *(ELEMENT->Nfve));    // free
	ELEMENT->VeCGmsh = calloc(NVEMAX        , sizeof *(ELEMENT->VeCGmsh)); // free
	ELEMENT->VeFcon  = calloc(NFMAX*NFVEMAX , sizeof *(ELEMENT->VeFcon));  // free

	// Operators

	// h-refinement related
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

	ELEMENT->Ghat_vS_vS = calloc(NP , sizeof *(ELEMENT->Ghat_vS_vS)); // free
	ELEMENT->GvShat_fS  = calloc(NP , sizeof *(ELEMENT->GvShat_fS));  // free
	ELEMENT->GfS_fIs    = calloc(NP , sizeof *(ELEMENT->GfS_fIs));    // free
	ELEMENT->GfS_fIc    = calloc(NP , sizeof *(ELEMENT->GfS_fIc));    // free

	ELEMENT->nOrd_fS   = calloc(NP , sizeof *(ELEMENT->nOrd_fS));  // free
	ELEMENT->nOrd_fIs  = calloc(NP , sizeof *(ELEMENT->nOrd_fIs)); // free
	ELEMENT->nOrd_fIc  = calloc(NP , sizeof *(ELEMENT->nOrd_fIc)); // free


	ELEMENT->I_vGs_vP[1]  = calloc(NP , sizeof **(ELEMENT->I_vGs_vP));
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

		ELEMENT->I_vGc_fS[P]  = calloc(NP , sizeof **(ELEMENT->I_vGc_fS));
		ELEMENT->I_vGc_fIs[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_fIs));
		ELEMENT->I_vGc_fIc[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_fIc));
		ELEMENT->I_vCs_fS[P]  = calloc(NP , sizeof **(ELEMENT->I_vCs_fS));
		ELEMENT->I_vCs_fIs[P] = calloc(NP , sizeof **(ELEMENT->I_vCs_fIs));
		ELEMENT->I_vCs_fIc[P] = calloc(NP , sizeof **(ELEMENT->I_vCs_fIc));
		ELEMENT->I_vCc_fS[P]  = calloc(NP , sizeof **(ELEMENT->I_vCc_fS));
		ELEMENT->I_vCc_fIs[P] = calloc(NP , sizeof **(ELEMENT->I_vCc_fIs));
		ELEMENT->I_vCc_fIc[P] = calloc(NP , sizeof **(ELEMENT->I_vCc_fIc));

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

		ELEMENT->Ghat_vS_vS[P] = calloc(NP , sizeof **(ELEMENT->Ghat_vS_vS));
		ELEMENT->GvShat_fS[P]  = calloc(NP , sizeof **(ELEMENT->GvShat_fS));
		ELEMENT->GfS_fIs[P]    = calloc(NP , sizeof **(ELEMENT->GfS_fIs));
		ELEMENT->GfS_fIc[P]    = calloc(NP , sizeof **(ELEMENT->GfS_fIc));

		if      (P == 0)    PbMin = P,   PbMax = P+1;
		else if (P == PMax) PbMin = P-1, PbMax = PMax;
		else                PbMin = P-1, PbMax = P+1;
		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			ELEMENT->ChiS_vS[P][Pb]    = calloc(NVREFSFMAX , sizeof ***(ELEMENT->ChiS_vS));
			ELEMENT->ChiS_vIs[P][Pb]   = calloc(NVREFSFMAX , sizeof ***(ELEMENT->ChiS_vIs));
			ELEMENT->ChiS_vIc[P][Pb]   = calloc(NVREFSFMAX , sizeof ***(ELEMENT->ChiS_vIc));

			ELEMENT->I_vGc_vS[P][Pb]  = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGc_vS));
			ELEMENT->I_vGc_vIs[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGc_vIs));
			ELEMENT->I_vGc_vIc[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGc_vIc));
			ELEMENT->I_vCs_vS[P][Pb]  = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCs_vS));
			ELEMENT->I_vCs_vIs[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCs_vIs));
			ELEMENT->I_vCs_vIc[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCs_vIc));
			ELEMENT->I_vCc_vS[P][Pb]  = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCc_vS));
			ELEMENT->I_vCc_vIs[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCc_vIs));
			ELEMENT->I_vCc_vIc[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCc_vIc));

			ELEMENT->Ihat_vS_vS[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->Ihat_vS_vS));

			ELEMENT->ChiS_vP[P][Pb]    = calloc(1          , sizeof ***(ELEMENT->ChiS_vP));

			ELEMENT->I_vGc_vP[P][Pb]  = calloc(1          , sizeof ***(ELEMENT->I_vGc_vP));

			ELEMENT->Is_Weak_VV[P][Pb]    = calloc(1 , sizeof ***(ELEMENT->Is_Weak_VV));
			ELEMENT->Ic_Weak_VV[P][Pb]    = calloc(1 , sizeof ***(ELEMENT->Ic_Weak_VV));
			if (P == Pb) {
				ELEMENT->ChiInvS_vS[P][Pb] = calloc(1          , sizeof ***(ELEMENT->ChiInvS_vS));

				ELEMENT->ICs[P][Pb] = calloc(1 , sizeof ***(ELEMENT->ICs));
				ELEMENT->ICc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->ICc));

				ELEMENT->I_vGs_vIs[1][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGs_vIs));
				ELEMENT->I_vGs_vIc[1][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGs_vIc));
				ELEMENT->I_vGs_vP[1][Pb]  = calloc(1          , sizeof ***(ELEMENT->I_vGs_vP));
				ELEMENT->I_vGs_vGc[1][Pb] = calloc(1          , sizeof ***(ELEMENT->I_vGs_vGc));
				ELEMENT->I_vGs_vCs[1][Pb] = calloc(1          , sizeof ***(ELEMENT->I_vGs_vCs));
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
			}

			ELEMENT->ChiS_fS[P][Pb]     = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fS));
			ELEMENT->ChiS_fIs[P][Pb]    = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fIs));
			ELEMENT->ChiS_fIc[P][Pb]    = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fIc));
			ELEMENT->ChiS_fIs_sp[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fIs_sp));
			ELEMENT->ChiS_fIc_sp[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fIc_sp));

			ELEMENT->I_vGc_fS[P][Pb]  = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGc_fS));
			ELEMENT->I_vGc_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGc_fIs));
			ELEMENT->I_vGc_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGc_fIc));
			ELEMENT->I_vCs_fS[P][Pb]  = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCs_fS));
			ELEMENT->I_vCs_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCs_fIs));
			ELEMENT->I_vCs_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCs_fIc));
			ELEMENT->I_vCc_fS[P][Pb]  = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCc_fS));
			ELEMENT->I_vCc_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCc_fIs));
			ELEMENT->I_vCc_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCc_fIc));

			ELEMENT->Is_Weak_FF[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->Is_Weak_FF));
			ELEMENT->Ic_Weak_FF[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->Ic_Weak_FF));
			ELEMENT->Is_Weak_FF_sp[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->Is_Weak_FF_sp));
			ELEMENT->Ic_Weak_FF_sp[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->Ic_Weak_FF_sp));

			ELEMENT->Ghat_vS_vS[P][Pb] = calloc(NVREFMAX       , sizeof ***(ELEMENT->Ghat_vS_vS));
			ELEMENT->GvShat_fS[P][Pb]  = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->GvShat_fS));
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
	//unsigned int NP = DB.NP;

	struct S_VOLUME *VOLUME;
	VOLUME = malloc(sizeof *VOLUME); // free

	// Structures
	VOLUME->indexl = 0;
	VOLUME->indexg = 0;
	VOLUME->P      = 0;
	VOLUME->type   = 0;
	VOLUME->Eclass = 0;
	VOLUME->update = 0;
	VOLUME->curved = 0;
	VOLUME->hlevel = 0;

	VOLUME->neigh = calloc(NFMAX*NFREFMAX , sizeof *(VOLUME->neigh)); // free

	VOLUME->XYZ_vC = NULL; // free

	// Geometry
	VOLUME->NvnG  = 0;
	VOLUME->XYZ_S = NULL; // free
	VOLUME->XYZ   = NULL; // free

	VOLUME->detJV_vI = NULL; // free
	VOLUME->C_vC     = NULL; // free
	VOLUME->C_vI     = NULL; // free
	VOLUME->C_vf     = calloc(NFMAX , sizeof *(VOLUME->C_vf)); // free

	// Initialization
	VOLUME->NvnS = 0;
	VOLUME->What = NULL; // free
	VOLUME->RES  = NULL; // free

	// Solving
	VOLUME->RHS       = NULL; // free
	VOLUME->MInv      = NULL; // free

	// hp adaptivity
//	VOLUME->minRES = 0.0;
//	VOLUME->maxRES = 0.0;

	VOLUME->refine_current = 0; // ToBeDeleted: Potentially not needed.
	VOLUME->Vadapt         = 0;
	VOLUME->adapt_type     = 0;
	VOLUME->PNew           = 0;

	// structs
	VOLUME->next    = NULL;
	VOLUME->grpnext = NULL;

	return VOLUME;
}

struct S_FACET *New_FACET(void)
{
	struct S_FACET *FACET;
	FACET = malloc(sizeof *FACET); // free

	// Structures
	FACET->indexg = 0;
	FACET->P      = 0;
	FACET->type   = 0;
	FACET->VfIn   = 0;
	FACET->VfOut  = 0;
	FACET->BC     = 0;

	FACET->VIn   = NULL; // free (in memory_destructor_V)
	FACET->VOut  = NULL; // free (in memory_destructor_V)
	FACET->VfIn  = 0;
	FACET->VfOut = 0;

	FACET->IndOrdInOut = 0;
	FACET->IndOrdOutIn = 0;

	// Geometry
	FACET->curved  = 0;
	FACET->typeInt = 0;

	FACET->XYZ_fI   = NULL; // free
	FACET->XYZ_fS   = NULL; // free
	FACET->n_fI     = NULL; // free
	FACET->n_fS     = NULL; // free
	FACET->detJF_fI = NULL; // free
	FACET->detJF_fS = NULL; // free

	// Solving
	FACET->RHSIn  = NULL; // tbd
	FACET->RHSOut = NULL; // tbd

	FACET->next = NULL;

	return FACET;
}
