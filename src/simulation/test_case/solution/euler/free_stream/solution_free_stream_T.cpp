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
#include "definitions_intrusive.h"
#include "definitions_math.h"
#include "definitions_physics.h"
#include "definitions_solution.h"
#include "definitions_tol.h"


#include "def_templates_solution.h"
#include "def_templates_solution_euler.h"

#include "def_templates_multiarray.h"

#include "def_templates_math_functions.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a \ref Multiarray_T\* container holding the solution values at the input coordinates.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_free_stream
	(const struct const_Multiarray_T* xyz, ///< xyz coordinates at which to evaluate the solution.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_T* constructor_const_sol_free_stream_T
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim)
{
	struct Multiarray_T* sol = constructor_sol_free_stream(xyz,sim); // returned
	return (const struct const_Multiarray_T*) sol;
}

void set_sol_free_stream_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	const struct const_Multiarray_T* xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed
	struct Multiarray_T* sol = constructor_sol_free_stream(xyz,sim); // destructed
	destructor_const_Multiarray_T(xyz);

	update_Solution_Container_sol_T(&sol_cont,sol,sim);
	destructor_Multiarray_T(sol);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

struct Sol_Data__fs;

/** \brief Version of \ref mutable_constructor_sol_fptr_T with additional \ref Sol_Data__fs input.
 *  \return See brief.
 *
 *  \param xyz      See brief.
 *  \param sim      See brief.
 *  \param sol_data \ref Sol_Data__fs.
 */
typedef struct Multiarray_T* (*mutable_constructor_sol_fs_fptr_T)
	(const struct const_Multiarray_T*const xyz,
	 const struct Simulation*const sim,
	 const struct Sol_Data__fs*const sol_data
	);

/// \brief Container for solution data relating to 'f'ree 's'tream.
struct Sol_Data__fs {
	// Read parameters
	Real rho,   ///< The free stream density.
	     p,     ///< The free stream pressure.
	     mach,  ///< The free stream 'mach' number.
	     theta; ///< The free stream flow angle in the xy-plane (in radians).

	Real pert_mag; ///< 'Pert'urbation 'mag'nitude.

	mutable_constructor_sol_fs_fptr_T constructor_sol; ///< Pointer to the function used to construct the solution.
};

/** \brief Return the statically allocated \ref Sol_Data__fs container.
 *  \return See brief. */
static struct Sol_Data__fs get_sol_data
	( );

static struct Multiarray_T* constructor_sol_free_stream
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim)
{
	assert(DIM >= 2);
	const struct Sol_Data__fs sol_data = get_sol_data();
	return sol_data.cpponstructor_sol(xyz,sim,&sol_data);
}

// Level 1 ********************************************************************************************************** //

/// \brief Read the required solution data into \ref Sol_Data__fs.
static void read_data_free_stream
	(struct Sol_Data__fs*const sol_data ///< \ref Sol_Data__fs.
	);

static struct Sol_Data__fs get_sol_data ( )
{
	static bool need_input = true;

	static struct Sol_Data__fs sol_data;
	if (need_input) {
		need_input = false;
		read_data_free_stream(&sol_data);
	}

	return sol_data;
}

// Level 2 ********************************************************************************************************** //

/** \brief Version of \ref constructor_sol_free_stream with a constant profile along the boundaries.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_free_stream_const
	(const struct const_Multiarray_T*const xyz, ///< See brief.
	 const struct Simulation*const sim,         ///< See brief.
	 const struct Sol_Data__fs*const sol_data   ///< Defined for \ref mutable_constructor_sol_fs_fptr_T.
	);

/** \brief Version of \ref constructor_sol_free_stream with a trigonometric profile along the x-coordinate boundaries.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_free_stream_trig_x
	(const struct const_Multiarray_T*const xyz, ///< See brief.
	 const struct Simulation*const sim,         ///< See brief.
	 const struct Sol_Data__fs*const sol_data   ///< Defined for \ref mutable_constructor_sol_fs_fptr_T.
	);

static void read_data_free_stream (struct Sol_Data__fs*const sol_data)
{
	const int count_to_find = 5;

	int boundary_pert = -1;

	FILE* input_file = fopen_input('s',NULL,NULL); // closed

	int count_found = 0,
	    count_dummy = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("density",  &count_found,line,&sol_data->rho);
		read_skip_string_count_d("pressure", &count_found,line,&sol_data->p);
		read_skip_string_count_d("mach",     &count_found,line,&sol_data->mach);
		read_skip_string_count_d("theta_deg",&count_found,line,&sol_data->theta);
		read_skip_string_count_d("magnitude_perturb",    &count_dummy,line,&sol_data->pert_mag);
		read_skip_convert_i(line,"boundary_perturb_type",&boundary_pert,&count_found);
	}
	fclose(input_file);

	sol_data->theta *= PI/180.0;

	switch (boundary_pert) {
	case BOUNDARY_PERTURB_TYPE_NONE:
		sol_data->constructor_sol = constructor_sol_free_stream_const;
		break;
	case BOUNDARY_PERTURB_TYPE_TRIG_X:
		sol_data->constructor_sol = constructor_sol_free_stream_trig_x;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",boundary_pert);
		break;
	}

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}

// Level 3 ********************************************************************************************************** //

static struct Multiarray_T* constructor_sol_free_stream_const
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim, const struct Sol_Data__fs*const sol_data)
{
	const ptrdiff_t n_n = xyz->extents[0];
	assert(DIM == xyz->extents[1]);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,n_var}); // returned

	const Real rho_fs = sol_data->rho,
	           p_fs   = sol_data->p,
	           c_fs   = sqrt(GAMMA*p_fs/rho_fs);

	const Real V_fs   = sol_data->mach*c_fs,
	           u_fs   = V_fs*cos(sol_data->theta),
	           v_fs   = V_fs*sin(sol_data->theta);

	Type* rho = get_col_Multiarray_T(0,sol),
	    * u   = get_col_Multiarray_T(1,sol),
	    * v   = get_col_Multiarray_T(2,sol),
	    * p   = get_col_Multiarray_T(n_var-1,sol);
	for (int i = 0; i < n_n; ++i) {
		rho[i] = rho_fs;
		u[i]   = u_fs;
		v[i]   = v_fs;
		p[i]   = p_fs;
	}

	if (DIM == 3) {
		Type* w = get_col_Multiarray_T(3,sol);
		for (int i = 0; i < n_n; ++i)
			w[i] = 0.0;
	}
	convert_variables_T(sol,'p','c');

	return sol;
}

static struct Multiarray_T* constructor_sol_free_stream_trig_x
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim, const struct Sol_Data__fs*const sol_data)
{
	const ptrdiff_t n_n = xyz->extents[0];
	assert(DIM == xyz->extents[1]);
	assert(DIM >= 2);
	assert(sol_data->pert_mag != 0.0);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	const Type*const x = get_col_const_Multiarray_T(0,xyz),
	          *const y = get_col_const_Multiarray_T(1,xyz);

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,n_var}); // returned

	const Real rho_fs = sol_data->rho,
	           p_fs   = sol_data->p,
	           c_fs   = sqrt(GAMMA*p_fs/rho_fs);
	const Real mag    = sol_data->pert_mag;

	const Real V_fs   = sol_data->mach*c_fs,
	           u_fs   = V_fs*cos(sol_data->theta),
	           v_fs   = V_fs*sin(sol_data->theta);
	assert(u_fs != 0.0); // For perturbation below.

	Type* rho = get_col_Multiarray_T(0,sol),
	    * u   = get_col_Multiarray_T(1,sol),
	    * v   = get_col_Multiarray_T(2,sol),
	    * p   = get_col_Multiarray_T(n_var-1,sol);
	for (int i = 0; i < n_n; ++i) {
		const Real scale = 1e-1*PI,
		           shift = -0.8712;
/// \todo Delete unused and clean up.
UNUSED(scale); UNUSED(shift); UNUSED(y); UNUSED(mag);
#if 0
		rho[i] = rho_fs + mag*sin(scale*y[i]+shift);
		u[i]   = u_fs;// + mag*cos(scale*y[i]+shift);
		v[i]   = v_fs;
//		p[i]   = p_fs + mag*sin(scale*y[i]+shift);
		p[i]   = pow(rho[i],GAMMA);
#else
UNUSED(v_fs);
UNUSED(x);
		rho[i] = rho_fs;//+ mag*sin(scale*y[i]+shift);
		p[i]   = p_fs;
		const Type c = sqrt_T(GAMMA*p[i]/rho[i]);
//		u[i]   = 2.0/GM1*c-5.4;
		u[i]   = -2.0/GM1*c+6.4;
#if TYPE_RC == TYPE_REAL
//printf("%f %f %f\n",c,u[i],u[i]/c);
#endif
		v[i]   = 0.0;
//		rho[i] = rho_fs;
//		p[i]   = p_fs;
#endif
	}

	if (DIM == 3) {
		EXIT_ADD_SUPPORT;
		Type* w = get_col_Multiarray_T(3,sol);
		for (int i = 0; i < n_n; ++i)
			w[i] = 0.0;
	}
	convert_variables_T(sol,'p','c');

	return sol;
}

#include "undef_templates_solution.h"
#include "undef_templates_solution_euler.h"

#include "undef_templates_multiarray.h"

#include "undef_templates_math_functions.h"
#include "undef_templates_test_case.h"
