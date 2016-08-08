// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__setup_geometry_h__INCLUDED
#define DPG__setup_geometry_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
 
#include "mkl.h"

#include "Database.h"
#include "Parameters.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "output_to_paraview.h"
#include "setup_ToBeCurved.h"
#include "setup_geom_factors.h"
#include "setup_normals.h"
#include "vertices_to_exact_geom.h"


extern void setup_geometry  (void);
extern void setup_FACET_XYZ (struct S_FACET *FACET);

#endif // DPG__setup_geometry_h__INCLUDED
