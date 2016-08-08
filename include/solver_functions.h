// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__solver_functions_h__INCLUDED
#define DPG__solver_functions_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"


extern void get_facet_ordering (const unsigned int d, const unsigned int IndOrd, const unsigned int FType,
                                const unsigned int Ns, const unsigned int Nn, const unsigned int *symms,
                                const double *rst, unsigned int *nOrd);

#endif // DPG__solver_functions_h__INCLUDED
