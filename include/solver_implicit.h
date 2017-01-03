// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_implicit_h__INCLUDED
#define DPG__solver_implicit_h__INCLUDED

#include "petscksp.h"


extern void solver_implicit (void);
extern void setup_KSP       (Mat A, KSP ksp);

#endif // DPG__solver_implicit_h__INCLUDED
