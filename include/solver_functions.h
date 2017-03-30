// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_functions_h__INCLUDED
#define DPG__solver_functions_h__INCLUDED

#include "S_VOLUME.h"
#include "S_OpCSR.h"

struct S_OPERATORS_V {
	unsigned int NvnI, NvnS, NvnS_SF, NvnI_SF;
	double       *ChiS_vI, **D_Weak, *I_Weak, *ChiS_vI_SF, **D_Weak_SF, *I_Weak_SF;

	struct S_OpCSR **D_Weak_sp;
};

struct S_VDATA {
	unsigned int P, Eclass;

	struct S_OPERATORS_V **OPS;
	struct S_VOLUME      *VOLUME;
};

extern void init_ops_VOLUME    (struct S_OPERATORS_V *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass);
extern void init_VDATA         (struct S_VDATA *VDATA, struct S_VOLUME *VOLUME);
extern void compute_W_vI       (struct S_VDATA *VDATA, double *W_vI);
extern void convert_between_rp (const unsigned int Nn, const unsigned int Nrc, const double *C, double *Ap, double *Ar,
                                const char *conv_type);

extern void finalize_VOLUME_Inviscid_Weak (const unsigned int Nrc, const double *Ar_vI, double *RLHS,
                                           const char *term_type, struct S_VDATA *VDATA);

#endif // DPG__solver_functions_h__INCLUDED
