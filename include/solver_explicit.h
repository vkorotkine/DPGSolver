// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_explicit_h__INCLUDED
#define DPG__solver_explicit_h__INCLUDED

#include "S_VOLUME.h"

extern void solver_explicit (void);
extern void enforce_positivity_highorder (struct S_VOLUME *VOLUME);

#endif // DPG__solver_explicit_h__INCLUDED
