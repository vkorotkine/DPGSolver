// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__compute_GradW_DG_h__INCLUDED
#define DPG__compute_GradW_DG_h__INCLUDED

#include "solver.h"

extern void compute_GradW_DG (const struct S_solver_info*const solver_info);
extern void free_GradW_DG    (const struct S_solver_info*const solver_info);

#endif // DPG__compute_GradW_DG_h__INCLUDED
