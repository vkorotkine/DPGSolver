// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__memory_constructors_matrix_h__INCLUDED
#define DPG__memory_constructors_matrix_h__INCLUDED

#include "matrix_structs.h"
#include "setup_operators_support.h"

extern struct S_MATRIX **    constructor2_mat (size_t const N0);
extern struct S_MATRIX ***   constructor3_mat (size_t const N0, size_t const N1);
extern struct S_MATRIX ****  constructor4_mat (size_t const N0, size_t const N1, size_t const N2);
extern struct S_MATRIX ***** constructor5_mat (size_t const N0, size_t const N1, size_t const N2, size_t const N3);

extern void constructor_move4_mat (char const layout, char const format,
                                   unsigned int const *const NRows1, unsigned int const *const NCols1,
                                   unsigned int const *const *const NRows2, unsigned int const *const *const NCols2,
								   double ****A_d, struct S_MATRIX ****A, struct S_OP_RANGE *const op_range);

extern void constructor_move5_mat (char const layout, char const format, unsigned int const *const NRows,
                                   unsigned int const *const NCols, double *****A_d, struct S_MATRIX *****A,
                                   struct S_OP_RANGE *const op_range);

#endif // DPG__memory_constructors_matrix_h__INCLUDED
