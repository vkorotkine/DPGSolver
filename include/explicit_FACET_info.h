// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__explicit_FACET_info_h__INCLUDED
#define DPG__explicit_FACET_info_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Database.h"
#include "Parameters.h"

#include "element_functions.h"
#include "sum_factorization.h"
#include "matrix_functions.h"
#include "boundary_conditions.h"
#include "fluxes_inviscid.h"
#include "array_swap.h"


extern void explicit_FACET_info (void);

#endif // DPG__explicit_FACET_info_h__INCLUDED
