// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__output_to_paraview_h__INCLUDED
#define DPG__output_to_paraview_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
 
#include "Database.h"
#include "Parameters.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "variable_functions.h"


extern void output_to_paraview (const char *OutputType);

#endif // DPG__output_to_paraview_h__INCLUDED
