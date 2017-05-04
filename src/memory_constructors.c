// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "memory_constructors.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#include "Parameters.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "adaptation.h"

/*
 *	Purpose:
 *		Allocate memory for and initialize new structures.
 *
 *	Comments:
 *		Split this function into memory_constructor for each struct type to minimize dependencies. (ToBeDeleted)
 *		Change all initializations from 0 to UINT_MAX. (ToBeDeleted)
 *		GfS_fIs/c may not be needed based on FACE_info testing. See comments at the start of explicit_FACE_Info.
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
	ELEMENT->NveP2   = UINT_MAX;
	ELEMENT->Nf      = UINT_MAX;
	ELEMENT->Nvref   = UINT_MAX;
	ELEMENT->NvrefSF = UINT_MAX;

	ELEMENT->Nfve    = calloc(NFMAX         , sizeof *(ELEMENT->Nfve));    // free
	ELEMENT->VeCGmsh = calloc(NVEMAX        , sizeof *(ELEMENT->VeCGmsh)); // free
	ELEMENT->VeEcon  = calloc(NEMAX*NEVEMAX , sizeof *(ELEMENT->VeEcon));  // free
	ELEMENT->VeFcon  = calloc(NFMAX*NFVEMAX , sizeof *(ELEMENT->VeFcon));  // free
	ELEMENT->NrefV   = calloc(NREFVVARMAX   , sizeof *(ELEMENT->NrefV));   // free

	// Operators

	// h-refinement related
	ELEMENT->NEhref  = UINT_MAX;
	ELEMENT->type_h  = calloc(NEHREFMAX      , sizeof *(ELEMENT->type_h));  // free
	ELEMENT->Nfref   = calloc(NFMAX          , sizeof *(ELEMENT->Nfref));   // free
	ELEMENT->NfMixed = calloc(NFMIXEDMAX     , sizeof *(ELEMENT->NfMixed)); // free
	ELEMENT->VeE     = calloc(NEREFMAX*NEMAX , sizeof *(ELEMENT->VeE));     // free
	ELEMENT->VeF     = calloc(NFREFMAX*NFMAX , sizeof *(ELEMENT->VeF));     // free
	ELEMENT->Nvve    = calloc(NVREFMAX       , sizeof *(ELEMENT->Nvve));    // free
	ELEMENT->VeV     = calloc(NVREFMAX       , sizeof *(ELEMENT->VeV));     // free

	// Normals
	ELEMENT->nr = calloc(NFMAX*DMAX , sizeof *(ELEMENT->nr)); // free

	// Plotting
	ELEMENT->connectivity  = calloc(NP , sizeof *(ELEMENT->connectivity));  // free
	ELEMENT->connect_types = calloc(NP , sizeof *(ELEMENT->connect_types)); // free
	ELEMENT->connectivityE = calloc(NP , sizeof *(ELEMENT->connectivityE)); // free
	ELEMENT->NvnP          = calloc(NP , sizeof *(ELEMENT->NvnP));          // free
	ELEMENT->connect_NE    = calloc(NP , sizeof *(ELEMENT->connect_NE));    // free

	// Operators
	ELEMENT->NvnGs = calloc(NP , sizeof *(ELEMENT->NvnGs)); // free
	ELEMENT->NvnG2 = calloc(NP , sizeof *(ELEMENT->NvnG2)); // free
	ELEMENT->NvnGc = calloc(NP , sizeof *(ELEMENT->NvnGc)); // free
	ELEMENT->NvnCs = calloc(NP , sizeof *(ELEMENT->NvnCs)); // free
	ELEMENT->NvnCc = calloc(NP , sizeof *(ELEMENT->NvnCc)); // free
	ELEMENT->NvnIs = calloc(NP , sizeof *(ELEMENT->NvnIs)); // free
	ELEMENT->NvnIc = calloc(NP , sizeof *(ELEMENT->NvnIc)); // free
	ELEMENT->NvnS  = calloc(NP , sizeof *(ELEMENT->NvnS));  // free
	ELEMENT->NfnG2 = calloc(NP , sizeof *(ELEMENT->NfnG2)); // free
	ELEMENT->NfnGc = calloc(NP , sizeof *(ELEMENT->NfnGc)); // free
	ELEMENT->NfnS  = calloc(NP , sizeof *(ELEMENT->NfnS));  // free
	ELEMENT->NfnIs = calloc(NP , sizeof *(ELEMENT->NfnIs)); // free
	ELEMENT->NfnIc = calloc(NP , sizeof *(ELEMENT->NfnIc)); // free
	ELEMENT->NenG2 = calloc(NP , sizeof *(ELEMENT->NenG2)); // free
	ELEMENT->NenGc = calloc(NP , sizeof *(ELEMENT->NenGc)); // free

	ELEMENT->w_fIs = calloc(NP , sizeof *(ELEMENT->w_fIs)); // free
	ELEMENT->w_fIc = calloc(NP , sizeof *(ELEMENT->w_fIc)); // free

	ELEMENT->w_vIs = calloc(NP , sizeof *(ELEMENT->w_vIs)); // free
	ELEMENT->w_vIc = calloc(NP , sizeof *(ELEMENT->w_vIc)); // free

	ELEMENT->ChiS_vP      = calloc(NP , sizeof *(ELEMENT->ChiS_vP));      // free
	ELEMENT->ChiS_vS      = calloc(NP , sizeof *(ELEMENT->ChiS_vS));      // free
	ELEMENT->ChiS_vIs     = calloc(NP , sizeof *(ELEMENT->ChiS_vIs));     // free
	ELEMENT->ChiS_vIc     = calloc(NP , sizeof *(ELEMENT->ChiS_vIc));     // free
	ELEMENT->ChiInvS_vS   = calloc(NP , sizeof *(ELEMENT->ChiInvS_vS));   // free
	ELEMENT->ChiInvGs_vGs = calloc(NP , sizeof *(ELEMENT->ChiInvGs_vGs)); // free

	ELEMENT->ChiBezInvS_vS = calloc(NP , sizeof *(ELEMENT->ChiBezInvS_vS)); // free

	ELEMENT->IG2 = calloc(NP , sizeof *(ELEMENT->IG2)); // free
	ELEMENT->IGc = calloc(NP , sizeof *(ELEMENT->IGc)); // free
	ELEMENT->ICs = calloc(NP , sizeof *(ELEMENT->ICs)); // free
	ELEMENT->ICc = calloc(NP , sizeof *(ELEMENT->ICc)); // free

	ELEMENT->TGs   = calloc(NP , sizeof *(ELEMENT->TGs));   // free
	ELEMENT->TS    = calloc(NP , sizeof *(ELEMENT->TS));    // free
	ELEMENT->TS_vB = calloc(NP , sizeof *(ELEMENT->TS_vB)); // free

	ELEMENT->TInvS_vB = calloc(NP , sizeof *(ELEMENT->TInvS_vB)); // free

	ELEMENT->VeMask = calloc(NP , sizeof *(ELEMENT->VeMask)); // free

	ELEMENT->GradChiS_vS  = calloc(NP , sizeof *(ELEMENT->GradChiS_vS));  // free
	ELEMENT->GradChiS_vIs = calloc(NP , sizeof *(ELEMENT->GradChiS_vIs)); // free
	ELEMENT->GradChiS_vIc = calloc(NP , sizeof *(ELEMENT->GradChiS_vIc)); // free

	ELEMENT->I_vGs_vP  = calloc(NP , sizeof *(ELEMENT->I_vGs_vP));  // free
	ELEMENT->I_vGs_vGs = calloc(NP , sizeof *(ELEMENT->I_vGs_vGs)); // free
	ELEMENT->I_vGs_vG2 = calloc(NP , sizeof *(ELEMENT->I_vGs_vG2)); // free
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
	ELEMENT->Fmask     = calloc(NP , sizeof *(ELEMENT->Fmask));     // free
	ELEMENT->I_vGc_fGc = calloc(NP , sizeof *(ELEMENT->I_vGc_fGc)); // free
	ELEMENT->I_vG2_fG2 = calloc(NP , sizeof *(ELEMENT->I_vG2_fG2)); // free
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

	ELEMENT->I_vGc_eGc = calloc(NP , sizeof *(ELEMENT->I_vGc_eGc)); // free
	ELEMENT->I_vG2_eG2 = calloc(NP , sizeof *(ELEMENT->I_vG2_eG2)); // free

	ELEMENT->I_fGs_vG2 = calloc(NP , sizeof *(ELEMENT->I_fGs_vG2)); // free
	ELEMENT->I_fG2_vG2 = calloc(NP , sizeof *(ELEMENT->I_fG2_vG2)); // free
	ELEMENT->I_fGs_vGc = calloc(NP , sizeof *(ELEMENT->I_fGs_vGc)); // free
	ELEMENT->I_fGc_vGc = calloc(NP , sizeof *(ELEMENT->I_fGc_vGc)); // free

	ELEMENT->I_eGs_vG2 = calloc(NP , sizeof *(ELEMENT->I_eGs_vG2)); // free
	ELEMENT->I_eG2_vG2 = calloc(NP , sizeof *(ELEMENT->I_eG2_vG2)); // free
	ELEMENT->I_eGs_vGc = calloc(NP , sizeof *(ELEMENT->I_eGs_vGc)); // free
	ELEMENT->I_eGc_vGc = calloc(NP , sizeof *(ELEMENT->I_eGc_vGc)); // free

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


	ELEMENT->TGs[1] = calloc(NP , sizeof **(ELEMENT->TGs));

	ELEMENT->VeMask[1] = calloc(NP , sizeof **(ELEMENT->VeMask));

	ELEMENT->ChiInvGs_vGs[1] = calloc(NP , sizeof **(ELEMENT->ChiInvGs_vGs));

	ELEMENT->I_vGs_vP[1]  = calloc(NP , sizeof **(ELEMENT->I_vGs_vP));
	ELEMENT->I_vGs_vGs[1] = calloc(NP , sizeof **(ELEMENT->I_vGs_vGs));
	ELEMENT->I_vGs_vG2[1] = calloc(NP , sizeof **(ELEMENT->I_vGs_vG2));
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

	ELEMENT->I_fGs_vG2[1] = calloc(NP , sizeof **(ELEMENT->I_fGs_vG2));
	ELEMENT->I_fGs_vGc[1] = calloc(NP , sizeof **(ELEMENT->I_fGs_vGc));

	ELEMENT->I_eGs_vG2[1] = calloc(NP , sizeof **(ELEMENT->I_eGs_vG2));
	ELEMENT->I_eGs_vGc[1] = calloc(NP , sizeof **(ELEMENT->I_eGs_vGc));

	for (P = 0; P < NP; P++) {
		ELEMENT->NfnG2[P] = calloc(NESUBCMAX , sizeof **(ELEMENT->NfnG2));
		ELEMENT->NfnGc[P] = calloc(NESUBCMAX , sizeof **(ELEMENT->NfnGc));
		ELEMENT->NfnS[P]  = calloc(NESUBCMAX , sizeof **(ELEMENT->NfnS));
		ELEMENT->NfnIs[P] = calloc(NESUBCMAX , sizeof **(ELEMENT->NfnIs));
		ELEMENT->NfnIc[P] = calloc(NESUBCMAX , sizeof **(ELEMENT->NfnIc));

		ELEMENT->w_fIs[P] = calloc(NESUBCMAX , sizeof **(ELEMENT->w_fIs));
		ELEMENT->w_fIc[P] = calloc(NESUBCMAX , sizeof **(ELEMENT->w_fIc));

		ELEMENT->ChiS_vP[P]    = calloc(NP , sizeof **(ELEMENT->ChiS_vP));
		ELEMENT->ChiS_vS[P]    = calloc(NP , sizeof **(ELEMENT->ChiS_vS));
		ELEMENT->ChiS_vIs[P]   = calloc(NP , sizeof **(ELEMENT->ChiS_vIs));
		ELEMENT->ChiS_vIc[P]   = calloc(NP , sizeof **(ELEMENT->ChiS_vIc));
		ELEMENT->ChiInvS_vS[P] = calloc(NP , sizeof **(ELEMENT->ChiInvS_vS));

		ELEMENT->ChiBezInvS_vS[P] = calloc(NP , sizeof **(ELEMENT->ChiBezInvS_vS));

		ELEMENT->IG2[P] = calloc(NP , sizeof **(ELEMENT->IG2));
		ELEMENT->IGc[P] = calloc(NP , sizeof **(ELEMENT->IGc));
		ELEMENT->ICs[P] = calloc(NP , sizeof **(ELEMENT->ICs));
		ELEMENT->ICc[P] = calloc(NP , sizeof **(ELEMENT->ICc));

		ELEMENT->TS[P]    = calloc(NP , sizeof **(ELEMENT->TS));
		ELEMENT->TS_vB[P] = calloc(NP , sizeof **(ELEMENT->TS_vB));

		ELEMENT->TInvS_vB[P] = calloc(NP , sizeof **(ELEMENT->TInvS_vB));

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

		ELEMENT->GradChiS_vS[P]  = calloc(NP , sizeof **(ELEMENT->GradChiS_vS));
		ELEMENT->GradChiS_vIs[P] = calloc(NP , sizeof **(ELEMENT->GradChiS_vIs));
		ELEMENT->GradChiS_vIc[P] = calloc(NP , sizeof **(ELEMENT->GradChiS_vIc));

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

		ELEMENT->Fmask[P]     = calloc(NP , sizeof **(ELEMENT->Fmask));
		ELEMENT->I_vGc_fGc[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_fGc));
		ELEMENT->I_vG2_fG2[P] = calloc(NP , sizeof **(ELEMENT->I_vG2_fG2));
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

		ELEMENT->I_vGc_eGc[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_eGc));
		ELEMENT->I_vG2_eG2[P] = calloc(NP , sizeof **(ELEMENT->I_vG2_eG2));

		ELEMENT->I_fG2_vG2[P] = calloc(NP , sizeof **(ELEMENT->I_fG2_vG2));
		ELEMENT->I_fGc_vGc[P] = calloc(NP , sizeof **(ELEMENT->I_fGc_vGc));

		ELEMENT->I_eG2_vG2[P] = calloc(NP , sizeof **(ELEMENT->I_eG2_vG2));
		ELEMENT->I_eGc_vGc[P] = calloc(NP , sizeof **(ELEMENT->I_eGc_vGc));

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

		get_Pb_range(P,&PbMin,&PbMax);
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
				ELEMENT->ChiInvS_vS[P][Pb]   = calloc(1 , sizeof ***(ELEMENT->ChiInvS_vS));
				ELEMENT->ChiInvGs_vGs[1][Pb] = calloc(1 , sizeof ***(ELEMENT->ChiInvGs_vGs));

				ELEMENT->ChiBezInvS_vS[P][Pb] = calloc(1 , sizeof ***(ELEMENT->ChiBezInvS_vS));

				ELEMENT->IG2[P][Pb] = calloc(1 , sizeof ***(ELEMENT->IG2));
				ELEMENT->IGc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->IGc));
				ELEMENT->ICs[P][Pb] = calloc(1 , sizeof ***(ELEMENT->ICs));
				ELEMENT->ICc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->ICc));

				ELEMENT->TGs[1][Pb]   = calloc(1 , sizeof ***(ELEMENT->TGs));
				ELEMENT->TS[P][Pb]    = calloc(1 , sizeof ***(ELEMENT->TS));
				ELEMENT->TS_vB[P][Pb] = calloc(1 , sizeof ***(ELEMENT->TS_vB));

				ELEMENT->TInvS_vB[P][Pb] = calloc(1 , sizeof ***(ELEMENT->TInvS_vB));

				ELEMENT->VeMask[1][Pb] = calloc(NVREFMAX , sizeof ***(ELEMENT->VeMask));

				ELEMENT->I_vGs_vP[1][Pb]  = calloc(1          , sizeof ***(ELEMENT->I_vGs_vP));
				ELEMENT->I_vGs_vGs[1][Pb] = calloc(NVREFMAX   , sizeof ***(ELEMENT->I_vGs_vGs));
				ELEMENT->I_vGs_vG2[1][Pb] = calloc(NVREFMAX   , sizeof ***(ELEMENT->I_vGs_vG2));
				ELEMENT->I_vGs_vGc[1][Pb] = calloc(1          , sizeof ***(ELEMENT->I_vGs_vGc));
				ELEMENT->I_vGs_vCs[1][Pb] = calloc(1          , sizeof ***(ELEMENT->I_vGs_vCs));
				ELEMENT->I_vGs_vIs[1][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGs_vIs));
				ELEMENT->I_vGs_vIc[1][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGs_vIc));
				ELEMENT->I_vGs_vS[1][Pb]  = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGs_vS));
				ELEMENT->I_vGc_vCc[P][Pb] = calloc(1          , sizeof ***(ELEMENT->I_vGc_vCc));

				ELEMENT->GradChiS_vS[P][Pb]  = calloc(1 , sizeof ***(ELEMENT->GradChiS_vS));
				ELEMENT->GradChiS_vIs[P][Pb] = calloc(1 , sizeof ***(ELEMENT->GradChiS_vIs));
				ELEMENT->GradChiS_vIc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->GradChiS_vIc));

				ELEMENT->GradChiS_vS[P][Pb][0]  = calloc(d , sizeof ****(ELEMENT->GradChiS_vS));
				ELEMENT->GradChiS_vIs[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->GradChiS_vIs));
				ELEMENT->GradChiS_vIc[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->GradChiS_vIc));

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

				ELEMENT->I_fGs_vG2[1][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_fGs_vG2));
				ELEMENT->I_fGs_vGc[1][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_fGs_vGc));

				ELEMENT->I_eGs_vG2[1][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_eGs_vG2));
				ELEMENT->I_eGs_vGc[1][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_eGs_vGc));
			}
			ELEMENT->I_fG2_vG2[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_fG2_vG2));
			ELEMENT->I_fGc_vGc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_fGc_vGc));

			ELEMENT->I_eG2_vG2[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_eG2_vG2));
			ELEMENT->I_eGc_vGc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_eGc_vGc));

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

			ELEMENT->Fmask[P][Pb]     = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->Fmask));
			ELEMENT->I_vGc_fGc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGc_fGc));
			ELEMENT->I_vG2_fG2[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vG2_fG2));
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

			ELEMENT->I_vGc_eGc[P][Pb] = calloc(NEREFMAX*NEMAX , sizeof ***(ELEMENT->I_vGc_eGc));
			ELEMENT->I_vG2_eG2[P][Pb] = calloc(NEREFMAX*NEMAX , sizeof ***(ELEMENT->I_vG2_eG2));

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
	ELEMENT->ELEMENT_FACE = calloc(NFMIXEDMAX , sizeof *(ELEMENT->ELEMENTclass)); // free

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

	VOLUME->XYZ_vV  = NULL; // free
	VOLUME->XYZ_vVc = NULL; // free

	// Geometry
	VOLUME->VeInd  = calloc(NVEMAX         , sizeof *(VOLUME->VeInd));  // free
	VOLUME->VeInfo = calloc(NVEMAX*NVEINFO , sizeof *(VOLUME->VeInfo)); // free
	VOLUME->BC     = calloc(2     , sizeof *(VOLUME->BC));    // free
	VOLUME->BC[0]  = calloc(NFMAX , sizeof *(VOLUME->BC[0])); // free
	VOLUME->BC[1]  = calloc(NEMAX , sizeof *(VOLUME->BC[1])); // free
	VOLUME->NvnG   = UINT_MAX;
	VOLUME->XYZ_S  = NULL; // free
	VOLUME->XYZ    = NULL; // free

	VOLUME->detJV_vI = NULL; // free
	VOLUME->C_vC     = NULL; // free
	VOLUME->C_vI     = NULL; // free
	VOLUME->C_vf     = calloc(NFMAX , sizeof *(VOLUME->C_vf)); // free

	// Initialization
	VOLUME->NvnS = UINT_MAX;
	VOLUME->What = NULL; // free
	VOLUME->RES  = NULL; // free
	VOLUME->QhatV = calloc(d , sizeof *(VOLUME->QhatV)); // free
	VOLUME->Qhat  = calloc(d , sizeof *(VOLUME->Qhat));  // free

	VOLUME->QhatV_What = calloc(d , sizeof *(VOLUME->QhatV_What)); // free
	VOLUME->Qhat_What  = calloc(d , sizeof *(VOLUME->Qhat_What));  // free

	// Solving
	VOLUME->IndA  = UINT_MAX;
	VOLUME->nnz_d = UINT_MAX;
	VOLUME->nnz_o = UINT_MAX;
	VOLUME->RHS   = NULL; // free
	VOLUME->LHS   = NULL; // free
	VOLUME->LHSQ  = calloc(d , sizeof *(VOLUME->LHSQ)); // free

	VOLUME->MInv      = NULL; // free
	VOLUME->MInv_diag = NULL; // free

	// Linearization testing
	VOLUME->What_c = NULL; // free
	VOLUME->RHS_c  = NULL; // free

	VOLUME->qhat_c = calloc(d , sizeof *(VOLUME->qhat_c)); // free

	VOLUME->QhatV_c = calloc(d , sizeof *(VOLUME->QhatV_c)); // free
	VOLUME->Qhat_c  = calloc(d , sizeof *(VOLUME->Qhat_c));  // free

	// hp adaptivity
//	VOLUME->minRES = 0.0;
//	VOLUME->maxRES = 0.0;

	VOLUME->refine_current = 0; // ToBeDeleted: Potentially not needed.
	VOLUME->Vadapt         = 0;
	VOLUME->adapt_type     = UINT_MAX;
	VOLUME->PNew           = UINT_MAX;
	VOLUME->hrefine_type   = UINT_MAX;

	VOLUME->XYZ_vVP2 = NULL; // free

	// Poisson
	VOLUME->qhat      = calloc(d , sizeof *(VOLUME->qhat)); // free
	VOLUME->qhat_uhat = calloc(d , sizeof *(VOLUME->qhat)); // free

	// structs
	VOLUME->next    = NULL;
	VOLUME->grpnext = NULL;
	VOLUME->child0  = NULL; // free (in memory_free_children)
	VOLUME->parent  = NULL; // need not be freed
	VOLUME->FACE   = calloc(NFMAX*NSUBFMAX , sizeof *(VOLUME->FACE)); // free

	return VOLUME;
}

struct S_FACE *New_FACE(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	struct S_FACE *FACE;
	FACE = malloc(sizeof *FACE); // free

	// Structures
	FACE->indexg = UINT_MAX;
	FACE->P      = UINT_MAX;
	FACE->type   = UINT_MAX;
	FACE->BC     = UINT_MAX;
	FACE->Boundary = UINT_MAX;
	FACE->level  = 0;
	FACE->update = 0;
	FACE->adapt_type = UINT_MAX;

	FACE->VIn   = NULL; // free (in memory_destructor_V)
	FACE->VOut  = NULL; // free (in memory_destructor_V)
	FACE->VfIn  = UINT_MAX;
	FACE->VfOut = UINT_MAX;

	FACE->IndOrdInOut = UINT_MAX;
	FACE->IndOrdOutIn = UINT_MAX;

	// Geometry
	FACE->curved  = 0;
	FACE->typeInt = UINT_MAX;

	FACE->XYZ_fI   = NULL; // free
	FACE->XYZ_fS   = NULL; // free
	FACE->n_fI     = NULL; // free
	FACE->n_fS     = NULL; // free
	FACE->detJF_fI = NULL; // free
	FACE->detJF_fS = NULL; // free

	FACE->detJVIn_fI  = NULL; // free
	FACE->detJVOut_fI = NULL; // free

	// Solving
	FACE->RHSIn  = NULL; // free (in finalize_RHS)
	FACE->RHSOut = NULL; // free (in finalize_RHS)

	FACE->LHSInIn   = NULL; // free (in finalize_LHS)
	FACE->LHSOutIn  = NULL; // free (in finalize_LHS)
	FACE->LHSInOut  = NULL; // free (in finalize_LHS)
	FACE->LHSOutOut = NULL; // free (in finalize_LHS)

	FACE->QhatL = calloc(d , sizeof *(FACE->QhatL)); // free
	FACE->QhatR = calloc(d , sizeof *(FACE->QhatR)); // free

	FACE->Qhat_WhatLL = calloc(d , sizeof *(FACE->Qhat_WhatLL)); // free
	FACE->Qhat_WhatRL = calloc(d , sizeof *(FACE->Qhat_WhatRL)); // free
	FACE->Qhat_WhatLR = calloc(d , sizeof *(FACE->Qhat_WhatLR)); // free
	FACE->Qhat_WhatRR = calloc(d , sizeof *(FACE->Qhat_WhatRR)); // free

	// Poisson
	FACE->qhatIn  = calloc(d , sizeof *(FACE->qhatIn));  // free
	FACE->qhatOut = calloc(d , sizeof *(FACE->qhatOut)); // free
	FACE->qhat_uhatInIn   = calloc(d , sizeof *(FACE->qhat_uhatInIn));   // free
	FACE->qhat_uhatOutIn  = calloc(d , sizeof *(FACE->qhat_uhatOutIn));  // free
	FACE->qhat_uhatInOut  = calloc(d , sizeof *(FACE->qhat_uhatInOut));  // free
	FACE->qhat_uhatOutOut = calloc(d , sizeof *(FACE->qhat_uhatOutOut)); // free

	// Linearization testing
	FACE->RHSIn_c  = NULL; // free (in finalize_RHS_c)
	FACE->RHSOut_c = NULL; // free (in finalize_RHS_c)

	FACE->qhatIn_c  = calloc(d , sizeof *(FACE->qhatIn_c));  // free
	FACE->qhatOut_c = calloc(d , sizeof *(FACE->qhatOut_c)); // free

	FACE->QhatL_c = calloc(d , sizeof *(FACE->QhatL_c)); // free
	FACE->QhatR_c = calloc(d , sizeof *(FACE->QhatR_c)); // free

	FACE->next   = NULL;
	FACE->child0 = NULL; // free (in memory_free_children)
	FACE->parent = NULL; // need not be freed

	return FACE;
}
