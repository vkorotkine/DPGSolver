// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__setup_normals_h__INCLUDED
#define DPG__setup_normals_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
 
#include "Database.h"
#include "Parameters.h"

#include "element_functions.h"
#include "matrix_functions.h"


extern void setup_normals (struct S_FACET *FACET);

#endif // DPG__setup_normals_h__INCLUDED
