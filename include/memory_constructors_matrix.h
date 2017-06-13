// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__memory_constructors_matrix_h__INCLUDED
#define DPG__memory_constructors_matrix_h__INCLUDED

#include "matrix_structs.h"

extern struct S_MATRIX *constructor1_mat  (size_t const N0);
extern struct S_MATRIX **constructor2_mat (size_t const N0, size_t const N1);

#endif // DPG__memory_constructors_matrix_h__INCLUDED
