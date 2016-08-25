// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__finalize_RHS_c_h__INCLUDED
#define DPG__finalize_RHS_c_h__INCLUDED

#include "petscvec.h"


extern void finalize_RHS_c (void);
extern void assemble_RHS_c (Vec *b);

#endif // DPG__finalize_RHS_c_h__INCLUDED
