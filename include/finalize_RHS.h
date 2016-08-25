// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__finalize_RHS_h__INCLUDED
#define DPG__finalize_RHS_h__INCLUDED

#include "petscvec.h"


extern double finalize_RHS   (void);
extern void   finalize_Vec   (Vec *a, const unsigned int finalize_type);
extern void   assemble_RHS   (Vec *b);

#endif // DPG__finalize_RHS_h__INCLUDED
