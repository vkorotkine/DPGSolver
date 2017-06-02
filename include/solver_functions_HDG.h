// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_functions_HDG_h__INCLUDED
#define DPG__solver_functions_HDG_h__INCLUDED

#include "matrix_structs.h"
#include "S_VOLUME.h"

struct S_OPERATORS_V {
	struct S_MATRIX const *const *D_Weak;
};

extern struct S_OPERATORS_V *init_mat_ops_VOLUME (struct S_VOLUME const *const VOLUME);


#endif // DPG__solver_functions_HDG_h__INCLUDED
