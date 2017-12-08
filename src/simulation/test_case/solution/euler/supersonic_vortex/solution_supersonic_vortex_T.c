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
#include "definitions_test_case.h"
#include "definitions_tol.h"


#include "def_templates_solution.h"
#include "def_templates_solution_euler.h"

#include "def_templates_multiarray.h"

#include "def_templates_math_functions.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a \ref Multiarray_T\* container holding the solution values at the input coordinates.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_supersonic_vortex
	(const struct const_Multiarray_R* xyz, ///< xyz coordinates at which to evaluate the solution.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_T* constructor_const_sol_supersonic_vortex_T
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	struct Multiarray_T* sol = constructor_sol_supersonic_vortex(xyz,sim); // returned
	return (const struct const_Multiarray_T*) sol;
}

void set_sol_supersonic_vortex_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	const struct const_Multiarray_R* xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed
	struct Multiarray_T* sol = constructor_sol_supersonic_vortex(xyz,sim); // destructed
	destructor_const_Multiarray_R(xyz);

	update_Solution_Container_sol_T(&sol_cont,sol);
	destructor_Multiarray_T(sol);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for solution data relating to 's'upersonic 'v'ortex.
struct Sol_Data__sv {
	// Read parameters
	double r_i,   ///< The 'r'adius of the 'i'nternal cylinder.
	       m_i,   ///< The 'm'ach number at the 'i'nternal cylinder face.
	       rho_i, ///< The density (\f$ rho \f$) at the 'i'nternal cylinder face.
	       V_i;   ///< The magnitude of the 'V'elocity at the 'i'nternal cylinder face.
};

/** \brief Return the statically allocated \ref Sol_Data__sv container.
 *  \return See brief. */
static struct Sol_Data__sv get_sol_data
	(const struct Simulation* sim ///< \ref Simulation.
	);

static struct Multiarray_T* constructor_sol_supersonic_vortex
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	assert(DIM >= 2);
	const struct Sol_Data__sv sol_data = get_sol_data(sim);

	// Compute the solution
	const ptrdiff_t n_n = xyz->extents[0];
	assert(DIM == xyz->extents[1]);

	const Real* x = get_col_const_Multiarray_R(0,xyz),
	          * y = get_col_const_Multiarray_R(1,xyz);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,n_var}); // returned

	Type* rho = get_col_Multiarray_T(0,sol),
	    * u   = get_col_Multiarray_T(1,sol),
	    * v   = get_col_Multiarray_T(2,sol),
	    * p   = get_col_Multiarray_T(n_var-1,sol);
	for (int i = 0; i < n_n; ++i) {
		const Real r_i   = sol_data.r_i,
		           m_i   = sol_data.m_i,
		           rho_i = sol_data.rho_i,
		           V_i   = sol_data.V_i;

		const Real r  = sqrt(x[i]*x[i]+y[i]*y[i]),
		           t  = atan2(y[i],x[i]),
		           Vt = -V_i/r;

		rho[i] = rho_i*pow(1.0+0.5*GM1*m_i*m_i*(1.0-pow(r_i/r,2.0)),1.0/GM1);
		u[i]   = -sin(t)*Vt;
		v[i]   =  cos(t)*Vt;
		p[i]   = pow_T(rho[i],GAMMA)/GAMMA;
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

/// \brief Read the required solution data into \ref Sol_Data__sv.
static void read_data_supersonic_vortex
	(const char*const input_path,       ///< Defined in \ref fopen_input.
	 struct Sol_Data__sv*const sol_data ///< \ref Sol_Data__sv.
	);

static struct Sol_Data__sv get_sol_data (const struct Simulation* sim)
{
	static bool need_input = true;

	static struct Sol_Data__sv sol_data;
	if (need_input) {
		need_input = false;
		read_data_supersonic_vortex(sim->input_path,&sol_data);
	}

	return sol_data;
}

// Level 2 ********************************************************************************************************** //

static void read_data_supersonic_vortex (const char*const input_path, struct Sol_Data__sv*const sol_data)
{
	const int count_to_find = 4;

	FILE* input_file = fopen_input(input_path,'s'); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("r_i",  &count_found,line,&sol_data->r_i);
		read_skip_string_count_d("m_i",  &count_found,line,&sol_data->m_i);
		read_skip_string_count_d("rho_i",&count_found,line,&sol_data->rho_i);
		read_skip_string_count_d("V_i",  &count_found,line,&sol_data->V_i);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}
