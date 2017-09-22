// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_symmetric_functions_h__INCLUDED
#define DPG__solver_symmetric_functions_h__INCLUDED

#include "solver.h"
#include "S_VOLUME.h"

extern void correct_collocated_for_symmetry_local (struct S_LHS_info*const LHS_info);
extern void correct_collocated_for_symmetry   (void);
extern void correct_collocated_for_symmetry_c (struct S_VOLUME *const VOLUME_perturbed, bool const correct_all,
                                               bool const correct_V, bool const correct_F);

#endif // DPG__solver_symmetric_functions_h__INCLUDED
