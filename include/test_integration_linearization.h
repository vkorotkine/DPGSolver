// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__test_integration_linearization_h__INCLUDED
#define DPG__test_integration_linearization_h__INCLUDED

#include "petscmat.h"


extern void test_integration_linearization (int nargc, char **argv);
extern void compute_A_cs                   (Mat *A, Vec *b, Vec *x, const unsigned int assemble_type);

#endif // DPG__test_integration_linearization_h__INCLUDED
