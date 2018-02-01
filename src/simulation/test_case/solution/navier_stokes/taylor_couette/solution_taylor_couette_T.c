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
static struct Multiarray_T* constructor_sol_taylor_couette
	(const struct const_Multiarray_R* xyz, ///< xyz coordinates at which to evaluate the solution.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_T* constructor_const_sol_taylor_couette_T
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	struct Multiarray_T* sol = constructor_sol_taylor_couette(xyz,sim); // returned
	return (const struct const_Multiarray_T*) sol;
}

void set_sol_taylor_couette_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	const struct const_Multiarray_R* xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed
	struct Multiarray_T* sol = constructor_sol_taylor_couette(xyz,sim); // destructed
	destructor_const_Multiarray_R(xyz);

	update_Solution_Container_sol_T(&sol_cont,sol,sim);
	destructor_Multiarray_T(sol);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for solution data relating to 't'aylor-'c'ouette.
struct Sol_Data__tc {
	Real r_i,   ///< The 'r'adius of the 'i'nner cylinder.
	     r_o,   ///< The 'r'adius of the 'o'uter cylinder.
	     omega, ///< The angular velocity of the rotating cylinder.
	     t_i,   ///< The temperature of the 'i'nner cylinder.
	     p_i,   ///< The pressure of the 'i'nner cylinder.
	     mu,    ///< The constant viscosity value.
	     kappa, ///< The coefficient of thermal conductivity.
	     r_s;   ///< The specific gas constant.
};

/** \brief Return the statically allocated \ref Sol_Data__tc container.
 *  \return See brief. */
static struct Sol_Data__tc get_sol_data
	(const struct Simulation* sim ///< \ref Simulation.
	);

static struct Multiarray_T* constructor_sol_taylor_couette
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	assert(DIM == 2);
	const struct Sol_Data__tc sol_data = get_sol_data(sim);

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

	const Real r_i   = sol_data.r_i,
	           r_o   = sol_data.r_o,
	           omega = sol_data.omega,
	           t_i   = sol_data.t_i,
	           p_i   = sol_data.p_i,
	           mu    = sol_data.mu,
	           kappa = sol_data.kappa,
	           r_s   = sol_data.r_s,

	           c     = omega/(1.0/(r_i*r_i)-1.0/(r_o*r_o));

EXIT_ADD_SUPPORT; // Different solutions depending on the boundary condition.
	for (int i = 0; i < n_n; ++i) {
		const Real r  = sqrt(x[i]*x[i]+y[i]*y[i]),
		           th = atan2(y[i],x[i]),
		           Vt = c*(1.0/r-r/(r_o*r_o));
		           t  = t_i - 2.0*c*c/(r_o*r_o)*mu/kappa*log(r/r_i) - c*c*mu/kappa*(1.0/(r*r)-1.0/(r_i*r_i));

		rho[i] = rho_i*pow(1.0+0.5*GM1*m_i*m_i*(1.0-pow(r_i/r,2.0)),1.0/GM1);
		u[i]   = -sin(th)*Vt;
		v[i]   =  cos(th)*Vt;
EXIT_ADD_SUPPORT; // Proper reference below.
		// Illingworth(1950), p.8 notes that the pressure is nearly uniform => set p = p_i and compute rho using the
		// the ideal gas law.
		p[i]   = p_i;
		rho[i] = p_i/(r_s*t);
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

/// \brief Read the required solution data into \ref Sol_Data__tc.
static void read_data_taylor_couette
	(const char*const input_path,       ///< Defined in \ref fopen_input.
	 struct Sol_Data__tc*const sol_data ///< \ref Sol_Data__tc.
	);

static struct Sol_Data__tc get_sol_data (const struct Simulation* sim)
{
	static bool need_input = true;

	static struct Sol_Data__tc sol_data;
	if (need_input) {
		need_input = false;
		read_data_taylor_couette(sim->input_path,&sol_data);
	}

	return sol_data;
}

// Level 2 ********************************************************************************************************** //

static void read_data_taylor_couette (const char*const input_path, struct Sol_Data__tc*const sol_data)
{
	const int count_to_find = 8;

	FILE* input_file = fopen_input(input_path,'s',NULL); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("r_i",  &count_found,line,&sol_data->r_i);
		read_skip_string_count_d("r_o",  &count_found,line,&sol_data->r_o);
		read_skip_string_count_d("omega",&count_found,line,&sol_data->omega);
		read_skip_string_count_d("t_i",  &count_found,line,&sol_data->t_i);
		read_skip_string_count_d("p_i",  &count_found,line,&sol_data->p_i);
		read_skip_string_count_d("mu",   &count_found,line,&sol_data->mu);
		read_skip_string_count_d("kappa",&count_found,line,&sol_data->kappa);
		read_skip_string_count_d("r_s",  &count_found,line,&sol_data->r_s);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}
