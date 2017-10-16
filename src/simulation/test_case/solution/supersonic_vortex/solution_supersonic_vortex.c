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

#include "solution_supersonic_vortex.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_math.h"
#include "definitions_tol.h"

#include "multiarray.h"

#include "element_solution.h"
#include "solver_volume.h"

#include "file_processing.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solution.h"
#include "solution_euler.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for solution data relating to 'p'eriodic 'v'ortex.
struct Sol_Data__sv {
	// Read parameters

	// Additional parameters
};

/// \brief Read the required solution data into \ref Sol_Data__sv.
static void read_data_supersonic_vortex
	(const char*const input_path,       ///< Defined in \ref fopen_input.
	 struct Sol_Data__sv*const sol_data ///< \ref Sol_Data__sv.
	);

/// \brief Set the remaining required solution data of \ref Sol_Data__sv based on the read values.
static void set_data_supersonic_vortex
	(struct Sol_Data__sv*const sol_data ///< \ref Sol_Data__sv.
	);

/// \brief Set the centre xy coordinates of the supersonic vortex at the given time.
void set_xy_c
	(double* x_c,                         ///< The x-coordinate of the vortex centre.
	 double* y_c,                         ///< The y-coordinate of the vortex centre.
	 const struct Sol_Data__sv* sol_data, ///< \ref Sol_Data__sv.
	 const double time                    ///< \ref Test_Case::time.
	);

// Interface functions ********************************************************************************************** //

void compute_sol_coef_v_supersonic_vortex (const struct Simulation* sim, struct Solver_Volume* volume)
{
	assert(sim->d >= 2);

	const char type_out = 'c'; // 'c'onservative

	// Set solution data
	static bool need_input = true;

	static struct Sol_Data__sv sol_data;
	if (need_input) {
		need_input = false;
		read_data_supersonic_vortex(sim->input_path,&sol_data);
		set_data_supersonic_vortex(&sol_data);
	}

	// Compute the solution coefficients
	const struct const_Multiarray_d* xyz_vs = constructor_xyz_vs(sim,volume); // destructed

	const ptrdiff_t n_vs = xyz_vs->extents[0],
	                d    = xyz_vs->extents[1];

	const double* x = get_col_const_Multiarray_d(0,xyz_vs),
	            * y = get_col_const_Multiarray_d(1,xyz_vs);

	const int n_var = sim->test_case->n_var;

	struct Multiarray_d* sol_vs = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_vs,n_var}); // destructed

	double* rho = get_col_Multiarray_d(0,sol_vs),
	      * u   = get_col_Multiarray_d(1,sol_vs),
	      * v   = get_col_Multiarray_d(2,sol_vs),
	      * p   = get_col_Multiarray_d(n_var-1,sol_vs);
	for (int i = 0; i < n_vs; ++i) {
		UNUSED(rho); UNUSED(u); UNUSED(v); UNUSED(p); UNUSED(x); UNUSED(y);
EXIT_ADD_SUPPORT;
	}

	if (d == 3) {
		double* w = get_col_Multiarray_d(3,sol_vs);
		for (int i = 0; i < n_vs; ++i)
			w[i] = 0.0;
	}
	destructor_const_Multiarray_d(xyz_vs);

	convert_variables(sol_vs,'p',type_out);

// external function here
EXIT_ADD_SUPPORT;
	// Convert to coefficients
	const char op_format = 'd';

	struct Volume* base_volume = (struct Volume*) volume;
	struct const_Solution_Element* element = (struct const_Solution_Element*) base_volume->element;

	const int p_ref = volume->p_ref;

	const struct Operator* vc0_vs_vs = get_Multiarray_Operator(element->vc0_vs_vs,(ptrdiff_t[]){0,0,p_ref,p_ref});

	struct Multiarray_d* sol_coef = volume->sol_coef;

	resize_Multiarray_d(sol_coef,sol_vs->order,sol_vs->extents);
	mm_NN1C_Operator_Multiarray_d(
		vc0_vs_vs,(const struct const_Multiarray_d*)sol_vs,sol_coef,op_format,sol_coef->order,NULL,NULL);

	destructor_Multiarray_d(sol_vs);
}

void compute_sol_coef_f_supersonic_vortex (const struct Simulation* sim, struct Solver_Face* face)
{
	switch (sim->method) {
	case METHOD_DG:
		compute_sol_coef_f_do_nothing(sim,face);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void read_data_supersonic_vortex (const char*const input_path, struct Sol_Data__sv*const sol_data)
{
	const int count_to_find = 0;

	FILE* input_file = fopen_input(input_path,'s'); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
UNUSED(sol_data);
EXIT_ADD_SUPPORT;
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}

static void set_data_supersonic_vortex (struct Sol_Data__sv*const sol_data)
{
UNUSED(sol_data);
EXIT_ADD_SUPPORT;
}
