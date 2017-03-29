// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__explicit_VOLUME_info_h__INCLUDED
#define DPG__explicit_VOLUME_info_h__INCLUDED

#include "S_VOLUME.h"

struct S_OPERATORS {
	unsigned int NvnI, NvnS, NvnS_SF, NvnI_SF;
	double       *ChiS_vI, *ChiS_vI_SF, **D_Weak, *I_Weak;

	struct S_OpCSR **D_Weak_sp;
};

struct S_VDATA {
	unsigned int P, Eclass;

	struct S_OPERATORS **OPS;
	struct S_VOLUME    *VOLUME;
};

extern void explicit_VOLUME_info (void);
extern void init_ops_VOLUME      (struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass);
extern void init_VDATA           (struct S_VDATA *VDATA, struct S_VOLUME *VOLUME);
extern void compute_W_vI         (struct S_VDATA *VDATA, double *W_vI);

#endif // DPG__explicit_VOLUME_info_h__INCLUDED
