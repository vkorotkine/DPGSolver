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

#include "solution_periodic_vortex.h"

#include <assert.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_math.h"
#include "definitions_tol.h"

#include "multiarray.h"

#include "solver_volume.h"

#include "simulation.h"
#include "solution.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for solution data relating to 'p'eriodic 'v'ortex.
struct Sol_Data__pv {
	const double x_c, ///< The x-coordinate of the 'c'enter of the vortex.
	             y_c, ///< The y-coordinate of the 'c'enter of the vortex.
	             r_v, ///< A parameter related to the radial decay of the vortex.
};

/// \brief Read the required solution data into \ref Sol_Data__pv.
static void read_data_periodic_vortex
	(const char*const input_path,      ///< Defined in \ref fopen_input.
	 struct Sol_Data_pv*const sol_data ///< \ref Sol_Data__pv.
	);

/// \brief Set the centre xy coordinates of the periodic vortex at the given time.
void set_xy_c
	(double* x_c,                        ///< \ref Sol_Data__pv::x_c.
	 double* y_c,                        ///< \ref Sol_Data__pv::y_c.
	 const struct Sol_Data_pv* sol_data, ///< \ref Sol_Data__pv.
	 const double time                   ///< \ref Test_Case::time.
	);

// Interface functions ********************************************************************************************** //

void compute_sol_coef_v_periodic_vortex (const struct Simulation* sim, struct Solver_Volume* volume)
{
	// Set solution data
	static bool need_input = true;

	static struct Sol_Data__pv sol_data;
	if (need_input) {
		need_input = false;
		read_data_periodic_vortex(sim->input_path,&sol_data);
	}

	// Set the coordinates of the vortex centre depending on the time.
	double x_c = 0.0,
	       y_c = 0.0;

	set_xy_c(&x_c,&y_c,sol_data,sim->test_case->time);

	// Compute the solution coefficients
	const struct const_Multiarray_d* xyz_vs = constructor_xyz_vs(sim,volume); // destructed

	const ptrdiff_t n_vs = xyz_vs->extents[0],
	                d    = xyz_vs->extents[1];

	assert(d >= 2);

	const double* x = get_col_const_Multiarray_d(0,xyz_vs),
	            * y = get_col_const_Multiarray_d(1,xyz_vs);

	const int n_var = sim->test_case->n_var;

	struct Multiarray_d* sol_vs = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_vs,n_var}); // destructed

	double* rho = get_col_Multiarray_d(0,sol_vs),
	      * u   = get_col_Multiarray_d(1,sol_vs),
	      * v   = get_col_Multiarray_d(2,sol_vs),
	      * p   = get_col_Multiarray_d(n_var-1,sol_vs);
	for (int i = 0; i < n_vs; ++i) {
		const double r2 = (pow(x[i]-x_c,2.0)+pow(y[i]-y_c,2.0))/(r_c*r_c);

		rho[i] = rho_inf;
		u[i]   = u_inf - con*(y[i]-y_c)/(r_v*r_v)*exp(-0.5*r2);
		v[i]   = v_inf + con*(x[i]-x_c)/(r_v*r_v)*exp(-0.5*r2);
		p[i]   = p_inf - rho_inf*(con*con)/(2*r_v*r_v)*exp(-r2);
	}

	if (d == 3) {
		double* w = get_col_Multiarray_d(3,sol_vs);
		for (int i = 0; i < n_vs; ++i)
			w[i] = 0.0;
	}
// convert to conservative (in-place)
// val to coef operator.

	destructor_const_Multiarray_d(xyz_vs);
	destructor_Multiarray_d(sol_vs);

	UNUSED(volume);
EXIT_UNSUPPORTED;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void read_data_periodic_vortex (const char*const input_path, struct Sol_Data__pv*const sol_data)
{
	int       count_found   = 0;
	const int count_to_find = 3;

	FILE* input_file = fopen_input(input_path,'s'); // closed

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		if (strstr(line,"x_c")) {
			++count_found;
			read_skip_const_d(line,&sol_data->x_c,2,false);
		} else if (strstr(line,"y_c")) {
			++count_found;
			read_skip_const_d(line,&sol_data->y_c,2,false);
		} else if (strstr(line,"r_v")) {
			++count_found;
			read_skip_const_d(line,&sol_data->r_v,2,false);
		}
	}
	fclose(input_file);

	assert(fmod(sol_data->theta,PI/4.0) < EPS);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}

void set_xy_c (double* x_c, double* y_c, const struct Sol_Data_pv* sol_data, const double time)
{

}
