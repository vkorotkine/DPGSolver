// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__finalize_RHS_c_h__INCLUDED
#define DPG__finalize_RHS_c_h__INCLUDED

#include <stdbool.h>
#include "S_VOLUME.h"

extern void finalize_RHS_c (struct S_VOLUME *const VOLUME_perturbed, const bool compute_all);

#endif // DPG__finalize_RHS_c_h__INCLUDED
