/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/** \file
 */

#include "solution_advection_default.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "macros.h"
#include "definitions_tol.h"

#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "file_processing.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solution.h"
#include "solution_advection.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

///\{ \name Constants used for the definition of the solution.
#define SOURCE_M 2.15 ///< 'M'ultiplication constant.
#define SOURCE_A 0.23 ///< 'A'ddition constant.
///\}

/** \brief Version of \ref constructor_sol_fptr used for the default linear advection test cases.
 *  \return See brief. */
static struct Multiarray_d* constructor_sol_advection_default
	(const struct const_Multiarray_d* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief Constructor for the source evaluated at the input xyz coordinates.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_source_advection_default
	(const struct const_Multiarray_d* xyz, ///< The xyz coordinates.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void set_sol_advection_default (const struct Simulation* sim, struct Solution_Container sol_cont)
{
	const struct const_Multiarray_d* xyz = constructor_xyz_sol(sim,&sol_cont); // destructed
	struct Multiarray_d* sol = constructor_sol_advection_default(xyz,sim); // destructed
	destructor_const_Multiarray_d(xyz);

	update_Solution_Container_sol(&sol_cont,sol);
	destructor_Multiarray_d(sol);
}

const struct const_Multiarray_d* constructor_const_sol_advection_default
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	struct Multiarray_d* sol = constructor_sol_advection_default(xyz,sim); // returned
	return (const struct const_Multiarray_d*) sol;
}

void compute_source_rhs_advection_default
	(const struct Simulation* sim, const struct Solver_Volume* s_vol, struct Multiarray_d* rhs)
{
	const struct const_Multiarray_d* xyz_vc    = constructor_xyz_vc_interp(s_vol,sim);             // destructed
	const struct const_Multiarray_d* source_vc = constructor_source_advection_default(xyz_vc,sim); // destructed
	destructor_const_Multiarray_d(xyz_vc);

	const struct const_Vector_d j_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)source_vc,&j_det_vc,false);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';
	const struct Operator* tw0_vt_vc = get_operator__tw0_vt_vc(s_vol);
	mm_NNC_Operator_Multiarray_d(1.0,1.0,tw0_vt_vc,source_vc,rhs,op_format,2,NULL,NULL);
	destructor_const_Multiarray_d(source_vc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief 1d version of \ref constructor_sol_advection_default.
 *  \return See brief. */
static struct Multiarray_d* constructor_sol_advection_default_1d
	(const struct const_Multiarray_d* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief 2d version of \ref constructor_sol_advection_default.
 *  \return See brief. */
static struct Multiarray_d* constructor_sol_advection_default_2d
	(const struct const_Multiarray_d* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief 1d version of \ref constructor_source_advection_default.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_source_advection_default_1d
	(const struct const_Multiarray_d* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief 2d version of \ref constructor_source_advection_default.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_source_advection_default_2d
	(const struct const_Multiarray_d* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

static struct Multiarray_d* constructor_sol_advection_default
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	static bool parameters_set = false;
	static mutable_constructor_sol_fptr constructor_sol = NULL;

	if (!parameters_set) {
		parameters_set = true;
		switch (DIM) {
			case 1: constructor_sol = constructor_sol_advection_default_1d; break;
			case 2: constructor_sol = constructor_sol_advection_default_2d; break;
			default: EXIT_UNSUPPORTED; break;
		}
	}
	return constructor_sol(xyz,sim);
}

static const struct const_Multiarray_d* constructor_source_advection_default
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	static bool parameters_set = false;
	static constructor_sol_fptr constructor_source = NULL;

	if (!parameters_set) {
		parameters_set = true;
		switch (DIM) {
			case 1: constructor_source = constructor_source_advection_default_1d; break;
			case 2: constructor_source = constructor_source_advection_default_2d; break;
			default: EXIT_UNSUPPORTED; break;
		}
	}
	return constructor_source(xyz,sim);
}

// Level 1 ********************************************************************************************************** //

static struct Multiarray_d* constructor_sol_advection_default_1d
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	assert(DIM == 1);
	assert(xyz->extents[1] == DIM);

	// Compute the solution
	const ptrdiff_t n_vs = xyz->extents[0];
	const int n_var = sim->test_case->n_var;

	struct Multiarray_d* sol = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const double* x = get_col_const_Multiarray_d(0,xyz);

	double* u = get_col_Multiarray_d(0,sol);
	for (int i = 0; i < n_vs; ++i) {
		u[i] = sin(SOURCE_M*x[i]+SOURCE_A);
	}

	return sol;
}

static struct Multiarray_d* constructor_sol_advection_default_2d
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	assert(DIM == 2);
	assert(xyz->extents[1] == DIM);

	// Compute the solution
	const ptrdiff_t n_vs = xyz->extents[0];
	const int n_var = sim->test_case->n_var;

	struct Multiarray_d* sol = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const double* x = get_col_const_Multiarray_d(0,xyz),
	            * y = get_col_const_Multiarray_d(1,xyz);

	double* u = get_col_Multiarray_d(0,sol);
	for (int i = 0; i < n_vs; ++i) {
		u[i] = sin(SOURCE_M*x[i]+SOURCE_A)
		     * sin(SOURCE_M*y[i]+SOURCE_A);
	}

	return sol;
}

static const struct const_Multiarray_d* constructor_source_advection_default_1d
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	assert(DIM == 1);
	assert(xyz->extents[1] == DIM);

	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(sim->input_path,&sol_data);
	}

	const ptrdiff_t n_vs = xyz->extents[0];
	const int n_var = sim->test_case->n_var;

	struct Multiarray_d* source = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const double* x = get_col_const_Multiarray_d(0,xyz);

	double* s = get_col_Multiarray_d(0,source);
	const double* b_adv = sol_data.b_adv;
	for (int i = 0; i < n_vs; ++i) {
		s[i] = b_adv[0]*SOURCE_M*cos(SOURCE_M*x[i]+SOURCE_A);
	}

	return (struct const_Multiarray_d*)source;
}

static const struct const_Multiarray_d* constructor_source_advection_default_2d
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	assert(DIM == 2);
	assert(xyz->extents[1] == DIM);

	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(sim->input_path,&sol_data);
	}

	const ptrdiff_t n_vs = xyz->extents[0];
	const int n_var = sim->test_case->n_var;

	struct Multiarray_d* source = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const double* x = get_col_const_Multiarray_d(0,xyz),
	            * y = get_col_const_Multiarray_d(1,xyz);

	double* s = get_col_Multiarray_d(0,source);
	const double* b_adv = sol_data.b_adv;
	for (int i = 0; i < n_vs; ++i) {
		s[i] = b_adv[0]*SOURCE_M*cos(SOURCE_M*x[i]+SOURCE_A)*sin(SOURCE_M*y[i]+SOURCE_A)
		     + b_adv[1]*SOURCE_M*sin(SOURCE_M*x[i]+SOURCE_A)*cos(SOURCE_M*y[i]+SOURCE_A);
	}

	return (struct const_Multiarray_d*)source;
}
