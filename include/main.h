// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__main_h__INCLUDED
#define DPG__main_h__INCLUDED

#include "S_DB.h"
#include "Test.h"


#ifndef TEST


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h> // ToBeModified: Likely not use system headers for mpi/petsc
#include <petscksp.h>
 
#include "Parameters.h"

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


#else // Used if -DTEST is passed as a compilation flag


#include <stdio.h>
#include <time.h>

#include "test_unit_array_find_index.h"
#include "test_unit_array_norm.h"
#include "test_unit_array_sort.h"
#include "test_unit_array_swap.h"
#include "test_unit_math_functions.h"
#include "test_unit_matrix_functions.h"
#include "test_unit_bases.h"
#include "test_unit_grad_bases.h"
#include "test_unit_cubature.h"
#include "test_unit_find_periodic_connections.h"
#include "test_unit_sum_factorization.h"
#include "test_unit_plotting.h"
#include "test_unit_fluxes_inviscid.h"
#include "test_unit_jacobian_fluxes_inviscid.h"
#include "test_unit_get_facet_ordering.h"

#include "test_integration_L2_projections.h"
#include "test_integration_update_h.h"

#include "test_speed_array_swap.h"
#include "test_speed_mm_CTN.h"


#endif // End TEST

#endif // DPG__main_h__INCLUDED
