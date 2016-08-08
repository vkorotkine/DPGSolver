// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__compute_errors_h__INCLUDED
#define DPG__compute_errors_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "Database.h"
#include "Parameters.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "variable_functions.h"
#include "exact_solutions.h"


extern void compute_errors        (const struct S_VOLUME *VOLUME, double *L2Error2, double *Vol, unsigned int *DOF, const unsigned int solved);
extern void compute_errors_global (void);

#endif // DPG__compute_errors_h__INCLUDED
