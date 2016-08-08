// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__exact_solutions_h__INCLUDED
#define DPG__exact_solutions_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Database.h"
#include "Parameters.h"


extern void compute_exact_solution (const unsigned int Nn, double *XYZ, double *UEx, double *sEx, const unsigned int solved);

#endif // DPG__exact_solutions_h__INCLUDED
