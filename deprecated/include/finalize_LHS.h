// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__finalize_LHS_h__INCLUDED
#define DPG__finalize_LHS_h__INCLUDED

#include "petscvec.h"
#include "petscmat.h"

extern double finalize_LHS   (Mat *A, Vec *b, Vec *x, const unsigned int assemble_type);
extern void   initialize_KSP (Mat *A, Vec *b, Vec *x);
extern void   finalize_Vec   (Vec *a, const unsigned int finalize_type);
extern void   finalize_Mat   (Mat *A, const unsigned int finalize_type);
extern void   finalize_ksp   (Mat *A, Vec *b, Vec *x, const unsigned int finalize_type);
extern void   assemble_RHS   (Vec *b, Vec *x);

#endif // DPG__finalize_LHS_h__INCLUDED
