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
	(const struct const_Multiarray_R* xyz, ///< xyz coordinates at which to evaluate the solution.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_T* constructor_const_sol_free_stream_T
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	struct Multiarray_T* sol = constructor_sol_free_stream(xyz,sim); // returned
	return (const struct const_Multiarray_T*) sol;
}

void set_sol_free_stream_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	const struct const_Multiarray_R* xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed
	struct Multiarray_T* sol = constructor_sol_free_stream(xyz,sim); // destructed
	destructor_const_Multiarray_R(xyz);

	update_Solution_Container_sol_T(&sol_cont,sol,sim);
	destructor_Multiarray_T(sol);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for solution data relating to 'f'ree 's'tream.
struct Sol_Data__fs {
	// Read parameters
	Real mach,  ///< The free stream 'mach' number.
	     theta; ///< The free stream flow angle in the xy-plane (in radians).
};

/** \brief Return the statically allocated \ref Sol_Data__fs container.
 *  \return See brief. */
static struct Sol_Data__fs get_sol_data
	( );

static struct Multiarray_T* constructor_sol_free_stream
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	assert(DIM >= 2);
	const struct Sol_Data__fs sol_data = get_sol_data();

	// Compute the solution
	const ptrdiff_t n_n = xyz->extents[0];
	assert(DIM == xyz->extents[1]);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,n_var}); // returned

	const Real rho_fs = 1.0,
	           p_fs   = pow(rho_fs,GAMMA),
	           c_fs   = sqrt(GAMMA*p_fs/rho_fs);

	const Real V_fs   = sol_data.mach*c_fs,
	           u_fs   = V_fs*cos(sol_data.theta),
	           v_fs   = V_fs*sin(sol_data.theta);

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

static void read_data_free_stream (struct Sol_Data__fs*const sol_data)
{
	const int count_to_find = 2;

	FILE* input_file = fopen_input('s',NULL,NULL); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("mach", &count_found,line,&sol_data->mach);
		read_skip_string_count_d("theta_deg",&count_found,line,&sol_data->theta);
	}
	fclose(input_file);

	sol_data->theta *= PI/180.0;

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}
