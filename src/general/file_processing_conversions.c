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

#include "file_processing_conversions.h"

#include <assert.h>
#include <string.h>

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_bc.h"
#include "definitions_dpg.h"
#include "definitions_geometry.h"
#include "definitions_physics.h"
#include "definitions_numerical_flux.h"
#include "definitions_test_case.h"

#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/** \brief Get the value of the macro defined constant associated with the string variable input of the definition type.
 *  \return See brief. */
static int get_define
	(const char*const def_str, ///< String associated with macro variable definition.
	 const char*const def_type ///< Defined for \ref read_skip_convert_const_i.
	);

// Interface functions ********************************************************************************************** //

void read_skip_convert_i
	(const char*const line, const char*const def_type, int*const var, int*const count_found)
{
	if (!strstr(line,def_type))
		return;

	if (count_found)
		++*count_found;

	char def_str[STRLEN_MAX];
	sscanf(line,"%*s %s",def_str);

	*var = get_define(def_str,def_type);
}

void read_skip_convert_const_i
	(const char*const line, const char*const def_type, const int*const var, int*const count_found)
{
	read_skip_convert_i(line,def_type,(int*)var,count_found);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static int get_define (const char*const def_str, const char*const def_type)
{
	int def_i = -1;
	if (strcmp(def_type,"solver_proc") == 0) {
		if      (strcmp(def_str,"explicit")           == 0) def_i = SOLVER_E;
		else if (strcmp(def_str,"implicit")           == 0) def_i = SOLVER_I;
		else if (strcmp(def_str,"explicit->implicit") == 0) def_i = SOLVER_EI;
		else
			EXIT_ERROR("Unsupported: %s\n",def_str);
	} else if (strcmp(def_type,"solver_type_e") == 0) {
		if      (strcmp(def_str,"forward_euler") == 0) def_i = SOLVER_E_EULER;
		else if (strcmp(def_str,"ssp_rk_33")     == 0) def_i = SOLVER_E_SSP_RK_33;
		else if (strcmp(def_str,"ls_rk_54")      == 0) def_i = SOLVER_E_LS_RK_54;
		else
			EXIT_ERROR("Unsupported: %s\n",def_str);
	} else if (strcmp(def_type,"solver_type_i") == 0) {
		if      (strcmp(def_str,"direct")    == 0) def_i = SOLVER_I_DIRECT;
		else if (strcmp(def_str,"iterative") == 0) def_i = SOLVER_I_ITERATIVE;
		else
			EXIT_ERROR("Unsupported: %s\n",def_str);
	} else if (strcmp(def_type,"lhs_terms") == 0) {
		if      (strcmp(def_str,"full_newton") == 0) def_i = LHS_FULL_NEWTON;
		else if (strcmp(def_str,"cfl_ramping") == 0) def_i = LHS_CFL_RAMPING;
		else
			EXIT_ERROR("Unsupported: %s\n",def_str);
	} else if (strcmp(def_type,"num_flux_1st") == 0) {
		if      (strcmp(def_str,"upwind")    == 0) def_i = NUM_FLUX_UPWIND;
		else if (strcmp(def_str,"Roe-Pike")  == 0) def_i = NUM_FLUX_ROE_PIKE;
		else                                       def_i = NUM_FLUX_INVALID;
	} else if (strcmp(def_type,"num_flux_2nd") == 0) {
		if      (strcmp(def_str,"BR2_stable") == 0) def_i = NUM_FLUX_BR2_STABLE;
		else if (strcmp(def_str,"CDG2") == 0)       def_i = NUM_FLUX_CDG2;
		else                                        def_i = NUM_FLUX_INVALID;
	} else if (strcmp(def_type,"test_norm") == 0) {
		if      (strcmp(def_str,"H0")        == 0) def_i = TEST_NORM_H0;
		else if (strcmp(def_str,"H1")        == 0) def_i = TEST_NORM_H1;
		else if (strcmp(def_str,"H1_upwind") == 0) def_i = TEST_NORM_H1_UPWIND;
		else                                       def_i = TEST_NORM_INVALID;
	} else if (strcmp(def_type,"geom_parametrization") == 0) {
		if      (strcmp(def_str,"radial_proj") == 0) def_i = GEOM_PRM_RADIAL_PROJ;
		else if (strcmp(def_str,"arc_length")  == 0) def_i = GEOM_PRM_ARC_LENGTH;
		else if (strcmp(def_str,"normal_proj") == 0) def_i = GEOM_PRM_NORMAL_PROJ;
		else
			EXIT_ERROR("Unsupported: %s\n",def_str);
	} else if (strcmp(def_type,"conservation") == 0) {
		if      (strcmp(def_str,"not_enforced")        == 0) def_i = CONSERVATION_NOT_ENFORCED;
		else if (strcmp(def_str,"lagrange_multiplier") == 0) def_i = CONSERVATION_LAGRANGE_MULT;
		else
			EXIT_ERROR("Unsupported: %s\n",def_str);
	} else if (strcmp(def_type,"viscosity_type") == 0) {
		if      (strcmp(def_str,"constant")   == 0) def_i = VISCOSITY_CONSTANT;
		else if (strcmp(def_str,"sutherland") == 0) def_i = VISCOSITY_SUTHERLAND;
		else
			EXIT_ERROR("Unsupported: %s\n",def_str);
	} else if (strcmp(def_type,"diabatic_flux_type") == 0) {
		if      (strcmp(def_str,"constant_zero") == 0) def_i = DIABATIC_FLUX_CONSTANT_ZERO;
		else if (strcmp(def_str,"constant")      == 0) def_i = DIABATIC_FLUX_CONSTANT;
		else
			EXIT_ERROR("Unsupported: %s\n",def_str);
	} else if (strcmp(def_type,"no_slip_type") == 0) {
		if      (strcmp(def_str,"rotating") == 0) def_i = NO_SLIP_ROTATING;
		else
			EXIT_ERROR("Unsupported: %s\n",def_str);
	} else if (strcmp(def_type,"no_slip_all_type") == 0) {
		if      (strcmp(def_str,"density_energy") == 0) def_i = NO_SLIP_ALL__RHO_E;
		else
			EXIT_ERROR("Unsupported: %s\n",def_str);
	} else {
		EXIT_ERROR("Unsupported: %s\n",def_type);
	}
	return def_i;
}
