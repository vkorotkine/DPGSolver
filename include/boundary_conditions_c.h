// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__boundary_conditions_c_h__INCLUDED
#define DPG__boundary_conditions_c_h__INCLUDED

#include <complex.h>

#include "boundary_conditions.h"

extern void set_BC_from_BType            (struct S_BC *const BCdata, char const *const BType);
extern void correct_XYZ_for_exact_normal (struct S_BC *const BCdata, char const *const BType);

extern void compute_boundary_values_c    (struct S_BC *const BCdata);

#endif // DPG__boundary_conditions_c_h__INCLUDED
