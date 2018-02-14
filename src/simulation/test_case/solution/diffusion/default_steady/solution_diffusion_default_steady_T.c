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
#include "definitions_math.h"
#include "definitions_tol.h"


#include "def_templates_solution.h"
#include "def_templates_solution_diffusion.h"

#include "def_templates_volume_solver.h"

#include "def_templates_multiarray.h"

#include "def_templates_operators.h"
#include "def_templates_solve.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Version of \ref constructor_sol_fptr_T used for the default steady diffusion test cases.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_diffusion_default_steady
	(const struct const_Multiarray_R* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief Version of \ref constructor_sol_fptr_T used for the default steady diffusion test cases.
 *  \return See brief. */
static struct Multiarray_T* constructor_grad_diffusion_default_steady
	(const struct const_Multiarray_R* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief Constructor for the source evaluated at the input xyz coordinates.
 *  \return See brief.
 *
 *  \note As the diffusion equation is being solved as if it were a conservation equation with \f$ f = -\nabla u \f$,
 *  the source term is given by \f$ s = -\nabla \cdot \nabla u \f$ (i.e. with a "-ve" sign).
 */
static const struct const_Multiarray_T* constructor_source_diffusion_default_steady
	(const struct const_Multiarray_R* xyz, ///< The xyz coordinates.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void set_sol_diffusion_default_steady_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	assert(sol_cont.node_kind != 'r');
	const struct const_Multiarray_R* xyz = constructor_xyz_sol_T(sim,&sol_cont);  // destructed
	struct Multiarray_T* sol = constructor_sol_diffusion_default_steady(xyz,sim); // destructed
	destructor_const_Multiarray_R(xyz);

	update_Solution_Container_sol_T(&sol_cont,sol,sim);
	destructor_Multiarray_T(sol);
}

void set_grad_diffusion_default_steady_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	assert(sol_cont.node_kind != 's');
	const struct const_Multiarray_R* xyz = constructor_xyz_sol_T(sim,&sol_cont);  // destructed
	struct Multiarray_T* grad = constructor_grad_diffusion_default_steady(xyz,sim); // destructed
	destructor_const_Multiarray_R(xyz);

	update_Solution_Container_grad_T(&sol_cont,grad,sim);
	destructor_Multiarray_T(grad);
}

const struct const_Multiarray_T* constructor_const_sol_diffusion_default_steady_T
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	struct Multiarray_T* sol = constructor_sol_diffusion_default_steady(xyz,sim); // returned
	return (const struct const_Multiarray_T*) sol;
}

const struct const_Multiarray_T* constructor_const_grad_diffusion_default_steady_T
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	struct Multiarray_T* grad = constructor_grad_diffusion_default_steady(xyz,sim); // returned
	return (const struct const_Multiarray_T*) grad;
}

void compute_source_rhs_diffusion_default_steady_T
	(const struct Simulation* sim, const struct Solver_Volume_T* s_vol, struct Multiarray_T* rhs)
{
	const struct const_Multiarray_R* xyz_vc    = constructor_xyz_vc_interp_T(s_vol,sim);                    // dest.
	const struct const_Multiarray_T* source_vc = constructor_source_diffusion_default_steady(xyz_vc,sim); // dest.
	destructor_const_Multiarray_R(xyz_vc);

	const struct const_Vector_d j_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
	scale_Multiarray_T_by_Vector_R('L',1.0,(struct Multiarray_T*)source_vc,&j_det_vc,false);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';
	const struct Operator* tw0_vt_vc = get_operator__tw0_vt_vc_T(s_vol);
	mm_NNC_Operator_Multiarray_T(1.0,1.0,tw0_vt_vc,source_vc,rhs,op_format,2,NULL,NULL);
	destructor_const_Multiarray_T(source_vc);
}

void add_to_flux_imbalance_source_diffusion_default_steady_T
	(const struct Simulation* sim, const struct Solver_Volume_T* s_vol, struct Multiarray_T* rhs)
{
	UNUSED(rhs);
	const struct const_Multiarray_R* xyz_vc    = constructor_xyz_vc_interp_T(s_vol,sim);                    // dest.
	const struct const_Multiarray_T* source_vc = constructor_source_diffusion_default_steady(xyz_vc,sim); // dest.
	destructor_const_Multiarray_R(xyz_vc);

	const struct const_Vector_d j_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
	scale_Multiarray_T_by_Vector_R('L',1.0,(struct Multiarray_T*)source_vc,&j_det_vc,false);
	add_to_flux_imbalance_source_T(source_vc,s_vol,sim);
	destructor_const_Multiarray_T(source_vc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief 1d version of \ref constructor_sol_diffusion_default_steady.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_diffusion_default_steady_1d
	(const struct const_Multiarray_R* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief 1d version of \ref constructor_grad_diffusion_default_steady.
 *  \return See brief. */
static struct Multiarray_T* constructor_grad_diffusion_default_steady_1d
	(const struct const_Multiarray_R* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief 1d version of \ref constructor_source_diffusion_default_steady.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_source_diffusion_default_steady_1d
	(const struct const_Multiarray_R* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief 2d version of \ref constructor_sol_diffusion_default_steady.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_diffusion_default_steady_2d
	(const struct const_Multiarray_R* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief 2d version of \ref constructor_grad_diffusion_default_steady.
 *  \return See brief. */
static struct Multiarray_T* constructor_grad_diffusion_default_steady_2d
	(const struct const_Multiarray_R* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

/** \brief 2d version of \ref constructor_source_diffusion_default_steady.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_source_diffusion_default_steady_2d
	(const struct const_Multiarray_R* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

static struct Multiarray_T* constructor_sol_diffusion_default_steady
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	static bool parameters_set = false;
	static mutable_constructor_sol_fptr_T constructor_sol = NULL;

	if (!parameters_set) {
		parameters_set = true;
		switch (DIM) {
			case 1: constructor_sol = constructor_sol_diffusion_default_steady_1d; break;
			case 2: constructor_sol = constructor_sol_diffusion_default_steady_2d; break;
			default: EXIT_UNSUPPORTED; break;
		}
	}
	return constructor_sol(xyz,sim);
}

static struct Multiarray_T* constructor_grad_diffusion_default_steady
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	static bool parameters_set = false;
	static mutable_constructor_sol_fptr_T constructor_grad = NULL;

	if (!parameters_set) {
		parameters_set = true;
		switch (DIM) {
			case 1: constructor_grad = constructor_grad_diffusion_default_steady_1d; break;
			case 2: constructor_grad = constructor_grad_diffusion_default_steady_2d; break;
			default: EXIT_UNSUPPORTED; break;
		}
	}
	return constructor_grad(xyz,sim);
}

static const struct const_Multiarray_T* constructor_source_diffusion_default_steady
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	static bool parameters_set = false;
	static constructor_sol_fptr_T constructor_source = NULL;

	if (!parameters_set) {
		parameters_set = true;
		switch (DIM) {
			case 1: constructor_source = constructor_source_diffusion_default_steady_1d; break;
			case 2: constructor_source = constructor_source_diffusion_default_steady_2d; break;
			default: EXIT_UNSUPPORTED; break;
		}
	}
	return constructor_source(xyz,sim);
}

// Level 1 ********************************************************************************************************** //

/// \brief Container for solution data relating to 'd'efault 'd'iffusion test case.
struct Sol_Data__dd {
	Real scale; ///< The scale used in the definition of the solution.
};

/** \brief Return the statically allocated \ref Sol_Data__dd container.
 *  \return See brief. */
static struct Sol_Data__dd get_sol_data
	( );

static struct Multiarray_T* constructor_sol_diffusion_default_steady_1d
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	assert(DIM == 1);
	assert(xyz->extents[1] == DIM);

	const struct Sol_Data__dd sol_data = get_sol_data();

	const ptrdiff_t n_vs = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const Real* x = get_col_const_Multiarray_R(0,xyz);

	Type* u = get_col_Multiarray_T(0,sol);
	const Real scale = sol_data.scale;
	for (int i = 0; i < n_vs; ++i)
		u[i] = cos(scale*PI*x[i]);

	return sol;
}

static struct Multiarray_T* constructor_grad_diffusion_default_steady_1d
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	assert(DIM == 1);
	assert(xyz->extents[1] == DIM);

	const struct Sol_Data__dd sol_data = get_sol_data(sim);

	const ptrdiff_t n_vs = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* grad = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_vs,n_var,DIM}); // returned

	const Real* x = get_col_const_Multiarray_R(0,xyz);

	Type* g[] = { get_col_Multiarray_T(0*n_vs*n_var,grad), };
	const Real scale = sol_data.scale;
	for (int i = 0; i < n_vs; ++i) {
		g[0][i] = -1.0*scale*PI*sin(scale*PI*x[i]);
	}

	return grad;
}

static const struct const_Multiarray_T* constructor_source_diffusion_default_steady_1d
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	assert(DIM == 1);
	assert(xyz->extents[1] == DIM);

	const struct Sol_Data__dd sol_data = get_sol_data(sim);

	const ptrdiff_t n_vs = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* source = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const Real* x = get_col_const_Multiarray_R(0,xyz);

	Type* s = get_col_Multiarray_T(0,source);
	const Real scale = sol_data.scale;
	for (int i = 0; i < n_vs; ++i)
		s[i] = DIM*pow(scale*PI,2.0)*cos(scale*PI*x[i]);

	return (struct const_Multiarray_T*)source;
}

static struct Multiarray_T* constructor_sol_diffusion_default_steady_2d
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	assert(DIM == 2);
	assert(xyz->extents[1] == DIM);

	const struct Sol_Data__dd sol_data = get_sol_data(sim);

	const ptrdiff_t n_vs = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const Real*const x = get_col_const_Multiarray_R(0,xyz),
	          *const y = get_col_const_Multiarray_R(1,xyz);

	Type* u = get_col_Multiarray_T(0,sol);
	const Real scale = sol_data.scale;
	for (int i = 0; i < n_vs; ++i)
		u[i] = cos(scale*PI*x[i])*cos(scale*PI*y[i]);

	return sol;
}

static struct Multiarray_T* constructor_grad_diffusion_default_steady_2d
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	assert(DIM == 2);
	assert(xyz->extents[1] == DIM);

	const struct Sol_Data__dd sol_data = get_sol_data(sim);

	const ptrdiff_t n_vs = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* grad = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_vs,n_var,DIM}); // returned

	const Real*const x = get_col_const_Multiarray_R(0,xyz),
	          *const y = get_col_const_Multiarray_R(1,xyz);

	Type*const g[] = { get_col_Multiarray_T(0*n_var,grad),
	                   get_col_Multiarray_T(1*n_var,grad), };
	const Real scale = sol_data.scale;
	for (int i = 0; i < n_vs; ++i) {
		g[0][i] = -1.0*scale*PI*sin(scale*PI*x[i])*cos(scale*PI*y[i]);
		g[1][i] = -1.0*scale*PI*cos(scale*PI*x[i])*sin(scale*PI*y[i]);
	}

	return grad;
}

static const struct const_Multiarray_T* constructor_source_diffusion_default_steady_2d
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	assert(DIM == 2);
	assert(xyz->extents[1] == DIM);

	const struct Sol_Data__dd sol_data = get_sol_data(sim);

	const ptrdiff_t n_vs = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* source = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const Real* x = get_col_const_Multiarray_R(0,xyz),
	          * y = get_col_const_Multiarray_R(1,xyz);

	Type* s = get_col_Multiarray_T(0,source);
	const Real scale = sol_data.scale;
	for (int i = 0; i < n_vs; ++i)
		s[i] = DIM*pow(scale*PI,2.0)*cos(scale*PI*x[i])*cos(scale*PI*y[i]);

	return (struct const_Multiarray_T*)source;
}

// Level 2 ********************************************************************************************************** //

/// \brief Read the required solution data into \ref Sol_Data__dd.
static void read_data_default_diffusion
	(struct Sol_Data__dd*const sol_data ///< \ref Sol_Data__dd.
	);

static struct Sol_Data__dd get_sol_data ( )
{
	static bool need_input = true;

	static struct Sol_Data__dd sol_data;
	if (need_input) {
		need_input = false;
		read_data_default_diffusion(&sol_data);
	}

	return sol_data;
}

// Level 3 ********************************************************************************************************** //

static void read_data_default_diffusion (struct Sol_Data__dd*const sol_data)
{
	const int count_to_find = 1;

	FILE* input_file = fopen_input('s',NULL,NULL); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("scale",&count_found,line,&sol_data->scale);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}
