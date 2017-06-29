// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__explicit_GradW_c_h__INCLUDED
#define DPG__explicit_GradW_c_h__INCLUDED

#include <stdbool.h>
#include "S_VOLUME.h"

extern void explicit_GradW_c (struct S_VOLUME *const VOLUME_perturbed, bool const compute_all);

#endif // DPG__explicit_GradW_c_h__INCLUDED
