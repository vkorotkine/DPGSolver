// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_functions_HDG_h__INCLUDED
#define DPG__solver_functions_HDG_h__INCLUDED

#include "matrix_structs.h"
#include "fluxes_structs.h"
#include "S_VOLUME.h"

extern void convert_to_mat_V (struct S_VOLUME *const VOLUME, char const mem_op);

struct S_OPERATORS_V {
	struct S_MATRIX const *ChiS_vI,
	                      *const *D_Weak,
	                      *I_vG_vI;
};

struct S_VDATA {
	struct S_MATRIX *W_vI, *XYZ_vI;

	struct S_OPERATORS_V const *OPS;
	struct S_VOLUME      const *VOLUME;
};

extern struct S_OPERATORS_V *init_mat_ops_VOLUME (struct S_VOLUME const *const VOLUME);
extern void coef_to_values_vI_M (struct S_VDATA *const VDATA, char const coef_type, char const mem_op);
extern void compute_flux_inviscid_M (struct S_VDATA *const VDATA, struct S_FLUX_M *const FLUXDATA, char const imex_type,
                                     char const mem_op);


#endif // DPG__solver_functions_HDG_h__INCLUDED
