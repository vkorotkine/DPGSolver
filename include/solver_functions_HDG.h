// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_functions_HDG_h__INCLUDED
#define DPG__solver_functions_HDG_h__INCLUDED

#include "matrix_structs.h"
#include "fluxes_structs.h"
#include "S_VOLUME.h"

extern void convert_to_multiarray_V (struct S_VOLUME *const VOLUME, char const mem_op);

struct S_OPERATORS_V {
	struct S_MATRIX const *ChiS_vI,
	                      *const *D_Weak,
	                      *I_vG_vI;
};

struct S_VDATA {
	struct S_MULTI_ARRAY *W_vI, *XYZ_vI;

	struct S_OPERATORS_V const *OPS;
	struct S_VOLUME      const *VOLUME;
};

extern void set_VDATA (struct S_VDATA *const VDATA, struct S_VOLUME const *const VOLUME, char const mem_op);
extern void set_FLUXDATA (struct S_FLUX_MA *const FLUXDATA);
extern void coef_to_values_vI_MA (struct S_VDATA *const VDATA, char const coef_type, char const mem_op);
extern void compute_flux_inviscid_MA (struct S_VDATA *const VDATA, struct S_FLUX_MA *const FLUXDATA,
                                      char const imex_type, char const mem_op);
extern void compute_flux_ref_MA (struct S_MULTI_ARRAY const *const C, struct S_MULTI_ARRAY const *const Ap,
                                 struct S_MULTI_ARRAY **Ar, const char mem_op);
extern void finalize_VOLUME_Inviscid_Weak_MA (struct S_MULTI_ARRAY const *const Ar_vI, struct S_MULTI_ARRAY *const RLHS,
                                              char const imex_type, struct S_VDATA const *const VDATA);


#endif // DPG__solver_functions_HDG_h__INCLUDED
