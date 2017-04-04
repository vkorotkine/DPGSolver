// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_functions_h__INCLUDED
#define DPG__solver_functions_h__INCLUDED

#include "S_VOLUME.h"
#include "S_OpCSR.h"


// VOLUME structs and functions

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


// FACE structs and functions

struct S_OPERATORS_F {
	// ToBeModified
	unsigned int NvnS, NfnI, NfnS, NvnS_SF, NfnI_SF, NvnI_SF, *nOrdInOut, *nOrdOutIn;
	double       **ChiS_fI, **ChiS_vI, **I_Weak_FF, **I_Weak_VV, **ChiS_fS, **GfS_fI;

	struct S_OpCSR **ChiS_fI_sp, **I_Weak_FF_sp;
};

struct S_FDATA {
	unsigned int P, Vf, f, SpOp, Eclass, IndFType;

	struct S_OPERATORS_F **OPS;
	struct S_VOLUME      *VOLUME;
	struct S_FACE        *FACE;
};

struct S_NumericalFlux {
	const double *WL_fIL, *WR_fIL;

	double *nFluxNum_fIL, *dnFluxNumdWL_fIL, *dnFluxNumdWR_fIL;
};

extern void init_ops_FACE      (struct S_OPERATORS_F *OPS, const struct S_VOLUME *VOLUME, const struct S_FACE *FACE,
                                const unsigned int IndClass);
extern void init_FDATA         (struct S_FDATA *FDATA, struct S_FACE *FACE, const char side);
extern void compute_W_fI       (const struct S_FDATA *FDATA, double *W_fI);
extern void compute_WR_fIL     (const struct S_FDATA *FDATA, const double *WL_fIL, double *WR_fIL);
extern void compute_numerical_flux (const struct S_FDATA *FDATA, struct S_NumericalFlux *NFluxData,
                                    const char imex_type);

#endif // DPG__solver_functions_h__INCLUDED
