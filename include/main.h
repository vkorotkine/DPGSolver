// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__main_h__INCLUDED
#define DPG__main_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <petscksp.h>
 
#include "database.h"
#include "parameters.h"
#include "test.h"

#include "initialization.h"
#include "setup_parameters.h"
#include "setup_mesh.h"
#include "setup_operators.h"
#include "setup_structures.h"
#include "setup_geometry.h"
#include "initialize_test_case.h"
#include "solver_explicit.h"
//#include "solver_implicit.h"
#include "compute_errors.h"
#include "memory_free.h"

#endif // DPG__main_h__INCLUDED
