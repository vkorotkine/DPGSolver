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

#include "multiarray.h"

#include "element_solution.h"
#include "volume_solver.h"

#include "file_processing.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solution.h"
#include "solution_euler.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for solution data relating to 'p'eriodic 'v'ortex.
struct Sol_Data__pv {
	// Read parameters
	double r_v,   ///< A parameter related to the radial decay of the vortex.
	       theta, ///< The angle at which the vortex propagates in the xy-plane.
	       V_inf, ///< The reference speed of the vortex in the theta-direction.
	       p_inf, ///< The reference pressure.
	       t_inf, ///< The reference temperature.
	       Rg;    ///< The gas constant.

	// Additional parameters
	double period_l, ///< The 'l'ength of the period.
	       rho_inf,  ///< The reference density.
	       u_inf,    ///< The speed of the vortex in the x-direction.
	       v_inf,    ///< The speed of the vortex in the y-direction.
	       con;      ///< A scaling constant.
};

/** \brief Return the statically allocated \ref Sol_Data__pv container.
 *  \return See brief. */
static struct Sol_Data__pv get_sol_data
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Set the centre xy coordinates of the periodic vortex at the given time.
void set_xy_c
	(double* x_c,                         ///< The x-coordinate of the vortex centre.
	 double* y_c,                         ///< The y-coordinate of the vortex centre.
	 const struct Sol_Data__pv* sol_data, ///< \ref Sol_Data__pv.
	 const double time                    ///< \ref Test_Case::time.
	);

/// \brief Read the required solution data into \ref Sol_Data__pv.
static void read_data_periodic_vortex
	(const char*const input_path,       ///< Defined in \ref fopen_input.
	 struct Sol_Data__pv*const sol_data ///< \ref Sol_Data__pv.
	);

/// \brief Set the remaining required solution data of \ref Sol_Data__pv based on the read values.
static void set_data_periodic_vortex
	(struct Sol_Data__pv*const sol_data ///< \ref Sol_Data__pv.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_periodic_vortex_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Sol_Data__pv get_sol_data (const struct Simulation* sim)
{
	static bool need_input = true;

	static struct Sol_Data__pv sol_data;
	if (need_input) {
		need_input = false;
		read_data_periodic_vortex(sim->input_path,&sol_data);
		set_data_periodic_vortex(&sol_data);
	}

	return sol_data;
}

void set_xy_c (double* x_c, double* y_c, const struct Sol_Data__pv* sol_data, const double time)
{
	const double theta    = sol_data->theta;
	const double period_l = sol_data->period_l;
	const double V_inf    = sol_data->V_inf;

	const double period_frac = fmod(time*V_inf+0.5*period_l,period_l)/period_l;

	// As the solution is not actually specified by a periodic function, the solution should only be evaluated when it
	// is close to the centre of the domain.
	assert(fabs(period_frac-0.5) <= 0.05);

	*x_c = (period_frac-0.5)*period_l*cos(theta);
	*y_c = (period_frac-0.5)*period_l*sin(theta);
}

static void read_data_periodic_vortex (const char*const input_path, struct Sol_Data__pv*const sol_data)
{
	const int count_to_find = 6;

	FILE* input_file = fopen_input(input_path,'s'); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		if (strstr(line,"r_v")) {
			++count_found;
			read_skip_d(line,&sol_data->r_v,1,false);
		} else if (strstr(line,"theta")) {
			++count_found;
			read_skip_d(line,&sol_data->theta,1,false);
		} else if (strstr(line,"V_inf")) {
			++count_found;
			read_skip_d(line,&sol_data->V_inf,1,false);
		} else if (strstr(line,"p_inf")) {
			++count_found;
			read_skip_d(line,&sol_data->p_inf,1,false);
		} else if (strstr(line,"t_inf")) {
			++count_found;
			read_skip_d(line,&sol_data->t_inf,1,false);
		} else if (strstr(line,"Rg")) {
			++count_found;
			read_skip_d(line,&sol_data->Rg,1,false);
		}
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}

static void set_data_periodic_vortex (struct Sol_Data__pv*const sol_data)
{
	const double theta = sol_data->theta;
	assert(fmod(theta,PI/4.0) < EPS);

	if (fmod(sol_data->theta,PI/2.0) < EPS) {
		sol_data->period_l = 2.0;
	} else {
		sol_data->period_l = 2.0*sqrt(2.0);
	}

	sol_data->rho_inf = (sol_data->p_inf)/(sol_data->Rg*sol_data->t_inf);

	const double V_inf = sol_data->V_inf;
	sol_data->u_inf = V_inf*cos(theta);
	sol_data->v_inf = V_inf*sin(theta);

	if (fabs(sol_data->u_inf) < 1e-1*EPS)
		sol_data->u_inf = 1e-1*EPS;

	if (fabs(sol_data->v_inf) < 1e-1*EPS)
		sol_data->v_inf = 1e-1*EPS;

	sol_data->con = 0.1*V_inf;
}
