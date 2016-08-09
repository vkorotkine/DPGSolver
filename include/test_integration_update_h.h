// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__test_integration_update_h_h__INCLUDED
#define DPG__test_integration_update_h_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "mkl.h"
#include <mpi.h>
#include <petscksp.h>

#include "Database.h"
#include "Parameters.h"
#include "Macros.h"
#include "Test.h"

#include "test_support.h"
#include "initialization.h"
#include "setup_parameters.h"
#include "setup_mesh.h"
#include "setup_operators.h"
#include "setup_structures.h"
#include "setup_geometry.h"
#include "initialize_test_case.h"
#include "adaptation.h"
#include "array_norm.h"
#include "array_print.h"
#include "memory_free.h"
#include "element_functions.h"
#include "matrix_functions.h"


extern void test_integration_update_h (int nargc, char **argv);

#endif // DPG__test_integration_update_h_h__INCLUDED
