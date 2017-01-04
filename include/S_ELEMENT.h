// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__S_ELEMENT_h__INCLUDED
#define DPG__S_ELEMENT_h__INCLUDED

#include "S_OpCSR.h"

struct S_ELEMENT {
	// Mesh
	unsigned int present, type, d, Nve, NveP2, Ne, Nf, Nvref, NvrefSF, Eclass, NEhref,
	             Neve, *Nfve, *VeCGmsh, *VeEcon, *VeFcon, *NrefV, *type_h;

	// Operators
	unsigned int *connect_NE, *NvnP, *Nvve,
	             *NvnGs, *NvnG2, *NvnGc, *NvnCs, *NvnCc, *NvnJs, *NvnJc, *NvnIs, *NvnIc, *NvnS,
	             **NfnG2, **NfnGc, **NfnS, **NfnIs, **NfnIc,
	             *NenG2, *NenGc,
	             Neref, *Nfref, *NfMixed,
	             **connectivity, **connect_types,
	             ***nOrd_fS, ***nOrd_fIs, ***nOrd_fIc, ****Fmask, ****VeMask;
	double       **VeE, **VeF, **VeV, *nr,
	             **w_vIs, **w_vIc, ***w_fIs, ***w_fIc,
	             ****ChiS_vP, ****ChiS_vS, ****ChiS_vIs, ****ChiS_vIc,
	             ****ChiInvS_vS, ****ChiInvGs_vGs,
	             ****IGc, ****ICs, ****ICc,
	             ****TGs,
	             ****I_vGs_vP, ****I_vGs_vGs, ****I_vGs_vG2, ****I_vGs_vGc, ****I_vGs_vCs, ****I_vGs_vS, ****I_vGs_vIs, ****I_vGs_vIc,
	             ****I_vGc_vP,                               ****I_vGc_vCc, ****I_vGc_vS, ****I_vGc_vIs, ****I_vGc_vIc,
	             ****I_vCs_vS, ****I_vCs_vIs, ****I_vCs_vIc,
	             ****I_vCc_vS, ****I_vCc_vIs, ****I_vCc_vIc,
	             ****Ihat_vS_vS,
	             *****GradChiS_vS, *****GradChiS_vIs, *****GradChiS_vIc,
	             *****D_vGs_vCs, *****D_vGs_vIs,
	             *****D_vGc_vCc, *****D_vGc_vIc,
	             *****D_vCs_vCs,
	             *****D_vCc_vCc,
	             ****ChiS_fS, ****ChiS_fIs, ****ChiS_fIc,
	             *****GradChiS_fIs, *****GradChiS_fIc,
	             ****I_vGs_fS, ****I_vGs_fIs, ****I_vGs_fIc,
	             ****I_vGc_fGc, ****I_vG2_fG2, ****I_vGc_fS, ****I_vGc_fIs, ****I_vGc_fIc,
	             ****I_vGc_eGc, ****I_vG2_eG2,
	             ****I_vCs_fS, ****I_vCs_fIs, ****I_vCs_fIc,
	             ****I_vCc_fS, ****I_vCc_fIs, ****I_vCc_fIc,
	             *****D_vGs_fIs, *****D_vGs_fIc,
	             *****D_vGc_fIs, *****D_vGc_fIc,
	             ****I_fGs_vGc, ****I_fGc_vGc, ****I_fGs_vG2, ****I_fG2_vG2,
	             ****I_eGs_vGc, ****I_eGc_vGc, ****I_eGs_vG2, ****I_eG2_vG2,
	             ****Is_Weak_VV, ****Ic_Weak_VV,
	             ****Is_Weak_FF, ****Ic_Weak_FF,
	             *****Ds_Weak_VV, *****Dc_Weak_VV,
	             ****L2hat_vS_vS,
	             ****GfS_fIs, ****GfS_fIc;

	struct S_OpCSR ****ChiS_fIs_sp, ****ChiS_fIc_sp,
	               *****Ds_Weak_VV_sp, *****Dc_Weak_VV_sp,
	               ****Is_Weak_FF_sp, ****Ic_Weak_FF_sp;

	struct S_ELEMENT *next;
	struct S_ELEMENT **ELEMENTclass, **ELEMENT_FACE;
};

#endif // DPG__S_ELEMENT_h__INCLUDED
