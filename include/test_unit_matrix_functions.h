// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__test_unit_matrix_functions_h__INCLUDED
#define DPG__test_unit_matrix_functions_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>

#include "mkl.h"

#include "Parameters.h"
#include "Test.h"

#include "test_support.h"
#include "array_norm.h"
#include "matrix_functions.h"
#include "array_free.h"


extern void test_unit_matrix_diag     (void);
extern void test_unit_matrix_identity (void);
extern void test_unit_matrix_inverse  (void);
extern void test_unit_matrix_mm       (void);
extern void test_unit_convert_to_CSR  (void);

#endif // DPG__test_unit_matrix_functions_h__INCLUDED
