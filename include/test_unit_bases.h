// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__test_unit_bases_h__INCLUDED
#define DPG__test_unit_bases_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mkl.h"

#include "Parameters.h"
#include "Test.h"

#include "test_support.h"
#include "test_code_bases.h"
#include "bases.h"
#include "cubature.h"
#include "math_functions.h"
#include "matrix_functions.h"
#include "array_norm.h"
#include "array_print.h"


extern void test_unit_basis_TP  (void);
extern void test_unit_basis_SI  (void);
extern void test_unit_basis_PYR (void);

#endif // DPG__test_unit_bases_h__INCLUDED
