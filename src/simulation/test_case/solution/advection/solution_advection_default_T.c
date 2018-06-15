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

#include <assert.h>
#include <math.h>
#include <string.h>

#include "macros.h"
#include "definitions_tol.h"


#include "def_templates_solution.h"
#include "def_templates_solution_advection.h"

#include "def_templates_volume_solver.h"

#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_math_functions.h"
#include "def_templates_operators.h"
#include "def_templates_solve.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

///\{ \name Constants used for the definition of the solution.
#define SOURCE_M 2.15 ///< 'M'ultiplication constant.
#define SOURCE_A 0.23 ///< 'A'ddition constant.
///\}

/** \brief Version of \ref constructor_sol_fptr_T used for the default linear advection test cases.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_advection_default_T
	(const struct const_Multiarray_T* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief Constructor for the source evaluated at the input xyz coordinates.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_source_advection_default_T
	(const struct const_Multiarray_T* xyz, ///< The xyz coordinates.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void set_sol_advection_default_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	const struct const_Multiarray_T* xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed
	struct Multiarray_T* sol = constructor_sol_advection_default_T(xyz,sim); // destructed
	destructor_const_Multiarray_T(xyz);

	update_Solution_Container_sol_T(&sol_cont,sol,sim);
	destructor_Multiarray_T(sol);
}

const struct const_Multiarray_T* constructor_const_sol_advection_default_T
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim)
{
	struct Multiarray_T* sol = constructor_sol_advection_default_T(xyz,sim); // returned
	return (const struct const_Multiarray_T*) sol;
}

void compute_source_rhs_advection_default_T
	(const struct Simulation* sim, const struct Solver_Volume_T* s_vol, struct Multiarray_T* rhs)
{
	const struct const_Multiarray_T* xyz_vc    = constructor_xyz_vc_interp_T(s_vol,sim);             // destructed
	const struct const_Multiarray_T* source_vc = constructor_source_advection_default_T(xyz_vc,sim); // destructed
	destructor_const_Multiarray_T(xyz_vc);

	const struct const_Vector_T j_det_vc = interpret_const_Multiarray_as_Vector_T(s_vol->jacobian_det_vc);
	scale_Multiarray_by_Vector_T('L',1.0,(struct Multiarray_T*)source_vc,&j_det_vc,false);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';
	const struct Operator* tw0_vt_vc = get_operator__tw0_vt_vc_T(s_vol);
	mm_NNC_Operator_Multiarray_T(1.0,1.0,tw0_vt_vc,source_vc,rhs,op_format,2,NULL,NULL);
	destructor_const_Multiarray_T(source_vc);
}

void add_to_flux_imbalance_source_advection_default_T
	(const struct Simulation* sim, const struct Solver_Volume_T* s_vol, struct Multiarray_T* rhs)
{
	UNUSED(rhs);
	const struct const_Multiarray_T* xyz_vc    = constructor_xyz_vc_interp_T(s_vol,sim);             // destructed
	const struct const_Multiarray_T* source_vc = constructor_source_advection_default_T(xyz_vc,sim); // destructed
	destructor_const_Multiarray_T(xyz_vc);

	const struct const_Vector_T j_det_vc = interpret_const_Multiarray_as_Vector_T(s_vol->jacobian_det_vc);
	scale_Multiarray_by_Vector_T('L',1.0,(struct Multiarray_T*)source_vc,&j_det_vc,false);
	add_to_flux_imbalance_source_T(source_vc,s_vol,sim);
	destructor_const_Multiarray_T(source_vc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief 1d version of \ref constructor_sol_advection_default_T.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_advection_default_1d_T
	(const struct const_Multiarray_T* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief 2d version of \ref constructor_sol_advection_default_T.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_advection_default_2d_T
	(const struct const_Multiarray_T* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief 1d version of \ref constructor_source_advection_default_T.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_source_advection_default_1d_T
	(const struct const_Multiarray_T* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief 2d version of \ref constructor_source_advection_default_T.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_source_advection_default_2d_T
	(const struct const_Multiarray_T* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

static struct Multiarray_T* constructor_sol_advection_default_T
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim)
{
	static bool parameters_set = false;
	static mutable_constructor_sol_fptr_T constructor_sol = NULL;

	if (!parameters_set) {
		parameters_set = true;
		switch (DIM) {
			case 1: constructor_sol = constructor_sol_advection_default_1d_T; break;
			case 2: constructor_sol = constructor_sol_advection_default_2d_T; break;
			default: EXIT_UNSUPPORTED; break;
		}
	}
	return constructor_sol(xyz,sim);
}

static const struct const_Multiarray_T* constructor_source_advection_default_T
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim)
{
	static bool parameters_set = false;
	static constructor_sol_fptr_T constructor_source = NULL;

	if (!parameters_set) {
		parameters_set = true;
		switch (DIM) {
			case 1: constructor_source = constructor_source_advection_default_1d_T; break;
			case 2: constructor_source = constructor_source_advection_default_2d_T; break;
			default: EXIT_UNSUPPORTED; break;
		}
	}
	return constructor_source(xyz,sim);
}

// Level 1 ********************************************************************************************************** //

static struct Multiarray_T* constructor_sol_advection_default_1d_T
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim)
{
	assert(DIM == 1);
	assert(xyz->extents[1] == DIM);

	// Compute the solution
	const ptrdiff_t n_vs = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const Type* x = get_col_const_Multiarray_T(0,xyz);

	Type* u = get_col_Multiarray_T(0,sol);
	for (int i = 0; i < n_vs; ++i) {
		u[i] = sin(SOURCE_M*real_T(x[i])+SOURCE_A);
		/* u[i] = 1.0; UNUSED(x); */
	}

	return sol;
}

static struct Multiarray_T* constructor_sol_advection_default_2d_T
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim)
{
	assert(DIM == 2);
	assert(xyz->extents[1] == DIM);

	// Compute the solution
	const ptrdiff_t n_vs = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const Type* x = get_col_const_Multiarray_T(0,xyz),
	          * y = get_col_const_Multiarray_T(1,xyz);

	Type* u = get_col_Multiarray_T(0,sol);
	for (int i = 0; i < n_vs; ++i) {
		u[i] = sin(SOURCE_M*real_T(x[i])+SOURCE_A) *
		       sin(SOURCE_M*real_T(y[i])+SOURCE_A);
	}

	return sol;
}

static const struct const_Multiarray_T* constructor_source_advection_default_1d_T
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim)
{
	assert(DIM == 1);
	assert(xyz->extents[1] == DIM);

	static bool need_input = true;
	static struct Sol_Data__Advection_T sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection_T(&sol_data);
	}

	const ptrdiff_t n_vs = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* source = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const Type* x = get_col_const_Multiarray_T(0,xyz);

	Type* s = get_col_Multiarray_T(0,source);
	for (int i = 0; i < n_vs; ++i) {
		const Type xyz_n[DIM] = ARRAY_DIM(x[i],0,0);
		const Real*const b_adv = sol_data.compute_b_adv(xyz_n);
		s[i] = b_adv[0]*SOURCE_M*cos(SOURCE_M*real_T(x[i])+SOURCE_A);
		/* s[i] = 0.0; UNUSED(b_adv); UNUSED(x); */
	}

	return (struct const_Multiarray_T*)source;
}

static const struct const_Multiarray_T* constructor_source_advection_default_2d_T
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim)
{
	assert(DIM == 2);
	assert(xyz->extents[1] == DIM);

	static bool need_input = true;
	static struct Sol_Data__Advection_T sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection_T(&sol_data);
	}

	const ptrdiff_t n_vs = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* source = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const Type* x = get_col_const_Multiarray_T(0,xyz),
	          * y = get_col_const_Multiarray_T(1,xyz);

	Type* s = get_col_Multiarray_T(0,source);
	for (int i = 0; i < n_vs; ++i) {
		const Type xyz_n[DIM] = ARRAY_DIM(x[i],y[i],0);
		const Real*const b_adv = sol_data.compute_b_adv(xyz_n);
		s[i] = b_adv[0]*SOURCE_M*cos(SOURCE_M*real_T(x[i])+SOURCE_A)*sin(SOURCE_M*real_T(y[i])+SOURCE_A) +
		       b_adv[1]*SOURCE_M*sin(SOURCE_M*real_T(x[i])+SOURCE_A)*cos(SOURCE_M*real_T(y[i])+SOURCE_A);
	}

	return (struct const_Multiarray_T*)source;
}

#include "undef_templates_solution.h"
#include "undef_templates_solution_advection.h"

#include "undef_templates_volume_solver.h"

#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_math_functions.h"
#include "undef_templates_operators.h"
#include "undef_templates_solve.h"
#include "undef_templates_test_case.h"
