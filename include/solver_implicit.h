// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_implicit_h__INCLUDED
#define DPG__solver_implicit_h__INCLUDED

#include <stdbool.h>
#include "petscksp.h"
#include "petscvec.h"


extern void solver_implicit (bool const PrintEnabled);
extern void setup_KSP       (Mat A, KSP ksp);

extern void solver_implicit_linear_system (Mat *A, Vec *b, Vec *x, KSP *ksp, unsigned int const iteration,
                                           bool const PrintEnabled);
extern void solver_implicit_update_What    (Vec x);

#endif // DPG__solver_implicit_h__INCLUDED
