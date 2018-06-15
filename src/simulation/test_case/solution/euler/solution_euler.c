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

#include "solution_euler.h"

#include "definitions_math.h"

#include "multiarray.h"

#include "boundary.h"
#include "compute_error_euler.h"
#include "const_cast.h"
#include "file_processing.h"
#include "flux_euler.h"
#include "geometry.h"
#include "geometry_parametric.h"
#include "numerical_flux_euler.h"
#include "simulation.h"
#include "solution.h"
#include "test_case.h"

#include "periodic_vortex/solution_periodic_vortex.h"
#include "supersonic_vortex/solution_supersonic_vortex.h"
#include "free_stream/solution_free_stream.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for solution data relating to 'c'oefficients of 'd'rag and 'l'ift.
struct Sol_Data__c_dl {
	// Read parameters
	double rho,      ///< The free stream density.
	       p,        ///< The free stream pressure.
	       mach,     ///< The free stream 'mach' number.
	       theta,    ///< The free stream flow angle in the xy-plane (in radians).
	       area_ref; ///< The reference area used to scale the coefficients.
};

/** \brief Return the statically allocated \ref Sol_Data__c_dl container.
 *  \return See brief. */
static struct Sol_Data__c_dl get_sol_data__c_dl
	( );

/** \brief Return the static variable holding the specific gas constant for the current test case.
 *  \return See brief. */
static double get_specific_gas_constant
	( );

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_euler_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "solution_euler_T.c"
#include "undef_templates_type.h"

void compute_entropy (struct Multiarray_d* s, const struct const_Multiarray_d* vars, const char var_type)
{
	assert(s->extents[0] == vars->extents[0]);
	assert(s->extents[1] == 1);
	assert(vars->layout == 'C');

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,var_type,'p');

	const double* rho = get_col_const_Multiarray_d(0,vars),
	            * p   = get_col_const_Multiarray_d(vars->extents[1]-1,vars);

	const ptrdiff_t ext_0 = s->extents[0];
	for (ptrdiff_t i = 0; i < ext_0; ++i)
//		s->data[i] = p[i]*pow(rho[i],-GAMMA);
		s->data[i] = log(p[i]*pow(rho[i],-GAMMA));

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,'p',var_type);
}

void compute_mach (struct Multiarray_d* mach, const struct const_Multiarray_d* vars, const char var_type)
{
	assert(mach->extents[0] == vars->extents[0]);
	assert(mach->extents[1] == 1);
	assert(vars->layout == 'C');

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,var_type,'p');

	const double* rho = get_col_const_Multiarray_d(0,vars),
		      *const uvw[DIM] = ARRAY_DIM( get_col_const_Multiarray_d(1,vars),
		                                   get_col_const_Multiarray_d(2,vars),
		                                   get_col_const_Multiarray_d(3,vars) ),
	            * p   = get_col_const_Multiarray_d(vars->extents[1]-1,vars);

	const ptrdiff_t ext_0 = mach->extents[0];
	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		double V2 = 0.0;
		for (int d = 0; d < DIM; ++d)
			V2 += uvw[d][i]*uvw[d][i];
		const double c2 = GAMMA*p[i]/rho[i];
		mach->data[i] = sqrt(V2/c2);
	}

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,'p',var_type);
}

void compute_temperature
	(struct Multiarray_d*const t, const struct const_Multiarray_d*const vars, const char var_type)
{
	assert(t->extents[0] == vars->extents[0]);
	assert(t->extents[1] == 1);
	assert(vars->layout == 'C');

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,var_type,'p');

	const double r_s = get_specific_gas_constant();
	const double* rho = get_col_const_Multiarray_d(0,vars),
	            * p   = get_col_const_Multiarray_d(vars->extents[1]-1,vars);

	const ptrdiff_t ext_0 = t->extents[0];
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		t->data[i] = p[i]/(r_s*rho[i]);

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,'p',var_type);
}

void compute_max_wavespeed (struct Multiarray_d* V_p_c, const struct const_Multiarray_d* vars, const char var_type)
{
	assert(V_p_c->extents[0] == vars->extents[0]);
	assert(V_p_c->extents[1] == 1);
	assert(vars->extents[1]-2 == DIM);
	assert(vars->layout == 'C');

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,var_type,'p');

	const double*const rho      = get_col_const_Multiarray_d(0,vars),
	            *const uvw[DIM] = ARRAY_DIM( get_col_const_Multiarray_d(1,vars),
		                                   get_col_const_Multiarray_d(2,vars),
		                                   get_col_const_Multiarray_d(3,vars) ),
	            *const p        = get_col_const_Multiarray_d(vars->extents[1]-1,vars);

	const ptrdiff_t ext_0 = V_p_c->extents[0];
	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		double V2 = 0.0;
		for (int d = 0; d < DIM; ++d)
			V2 += uvw[d][i]*uvw[d][i];
		const double c2 = GAMMA*p[i]/rho[i];
		V_p_c->data[i] = sqrt(V2)+sqrt(c2);
	}

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,'p',var_type);
}

void compute_cd_cl_values
	(struct Multiarray_d* c_dl, const struct const_Multiarray_d* vars, const char var_type,
	 const struct const_Multiarray_d*const normals)
{
	assert(c_dl->extents[0] == vars->extents[0]);
	assert(c_dl->extents[1] == 2);
	assert(vars->extents[1]-2 == DIM);
	assert(vars->layout == 'C');

	assert(DIM == 2); // Currently assuming that the free-stream flow direction is in the x-y plane. Add support.

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,var_type,'p');

	const struct Sol_Data__c_dl sol_data = get_sol_data__c_dl();

	const double rho_fs   = sol_data.rho,
	             p_fs     = sol_data.p,
	             c_fs     = sqrt(GAMMA*p_fs/rho_fs),
	             V_fs     = sol_data.mach*c_fs,
	             denom    = 0.5*rho_fs*V_fs*V_fs*sol_data.area_ref,
	             theta_fs = sol_data.theta;

	const double*const p = get_col_const_Multiarray_d(vars->extents[1]-1,vars);

	double*const cd = get_col_Multiarray_d(0,c_dl),
	      *const cl = get_col_Multiarray_d(1,c_dl);
	const ptrdiff_t ext_0 = c_dl->extents[0];
	for (int i = 0; i < ext_0; ++i) {
		const double*const n = get_row_const_Multiarray_d(i,normals);
		const double theta_n  = atan2(n[1],n[0]);
		const double theta_nf = theta_n-theta_fs;

		cd[i] = p[i]*cos(theta_nf)/denom;
		cl[i] = p[i]*sin(theta_nf)/denom;
	}

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,'p',var_type);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Read the required solution data into \ref Sol_Data__c_dl.
static void read_data_c_dl
	(struct Sol_Data__c_dl*const sol_data ///< \ref Sol_Data__c_dl.
	);

/// \brief Read the required data from the solution data file.
static void read_data_specific_gas
	(double*const r_s ///< Variable to be set to the value.
	);

static struct Sol_Data__c_dl get_sol_data__c_dl ( )
{
	static bool need_input = true;

	static struct Sol_Data__c_dl sol_data;
	if (need_input) {
		need_input = false;
		read_data_c_dl(&sol_data);
	}
	return sol_data;
}

static double get_specific_gas_constant ( )
{
	static double r_s = 0.0;

	static bool need_input = true;
	if (need_input) {
		need_input = false;
		read_data_specific_gas(&r_s);
	}
	return r_s;
}

// Level 1 ********************************************************************************************************** //

static void read_data_c_dl (struct Sol_Data__c_dl*const sol_data)
{
	const int count_to_find = 5;
	int count_found = 0;
	char line[STRLEN_MAX];
	FILE* input_file = NULL;

	input_file = fopen_input('s',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("density",  &count_found,line,&sol_data->rho);
		read_skip_string_count_d("pressure", &count_found,line,&sol_data->p);
		read_skip_string_count_d("mach",     &count_found,line,&sol_data->mach);
		read_skip_string_count_d("theta_deg",&count_found,line,&sol_data->theta);
	}
	fclose(input_file);

	sol_data->theta *= PI/180.0;

	input_file = fopen_input('g',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_c_style_d("area_ref",&count_found,line,&sol_data->area_ref);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}

static void read_data_specific_gas (double*const r_s)
{
	const int count_to_find = 1;

	FILE* input_file = fopen_input('s',NULL,NULL); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("r_s",&count_found,line,r_s);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}
