// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__initialize_test_case_h__INCLUDED
#define DPG__initialize_test_case_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
 
#include "Database.h"
#include "Parameters.h"
#include "Macros.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "variable_functions.h"
#include "exact_solutions.h"
#include "adaptation.h"
#include "compute_errors.h"


extern void initialize_test_case (const unsigned int adapt_update_MAX);

#endif // DPG__initialize_test_case_h__INCLUDED
