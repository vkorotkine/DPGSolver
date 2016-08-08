// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__setup_operators_h__INCLUDED
#define DPG__setup_operators_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
 
#include "mkl.h"

#include "Database.h"
#include "Parameters.h"
#include "Macros.h"

#include "bases.h"
#include "cubature.h"

#include "element_functions.h"
#include "sum_factorization.h"
#include "array_free.h"
#include "plotting_element_info.h"
#include "solver_functions.h"
#include "adaptation.h"
#include "memory_destructors.h"

extern void setup_operators (void);

#endif // DPG__setup_operators_h__INCLUDED
