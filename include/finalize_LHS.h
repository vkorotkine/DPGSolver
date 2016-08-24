// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__finalize_RHS_h__INCLUDED
#define DPG__finalize_RHS_h__INCLUDED

#include "petscmat.h"

extern void finalize_LHS   (Mat *A, Vec *b, const unsigned int assemble_type);
extern void initialize_KSP (Mat *A, Vec *b);
extern void compute_dof    (void);
extern void finalize_Mat   (Mat *A, const unsigned int finalize_type);

#endif // DPG__finalize_RHS_h__INCLUDED
