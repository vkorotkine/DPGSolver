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

#include "compute_error_navier_stokes.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_error.h"

#include "face.h"
#include "volume.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "compute_error.h"
#include "compute_error_euler.h"
#include "const_cast.h"
#include "element.h"
#include "intrusive.h"
#include "simulation.h"
#include "solution_euler.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a statically allocated `char*` holding the specific header for all of the velocity and temperature
 *         variables.
 *  \return See brief. */
static const char* compute_header_spec_navier_stokes_uvwt
	();

/** \brief Return a statically allocated `char*` holding the specific header for the velocity variables.
 *  \return See brief. */
static const char* compute_header_spec_navier_stokes_uvw
	();

// Interface functions ********************************************************************************************** //

struct Error_CE* constructor_Error_CE_navier_stokes_uvwt (const struct Simulation* sim)
{
	const int n_out = DIM+1; // Can add support for gradients if desired.

	struct Error_CE_Helper* e_ce_h = constructor_Error_CE_Helper(sim,n_out);
	e_ce_h->header_spec = compute_header_spec_navier_stokes_uvwt();
	const_cast_i(&e_ce_h->error_type,ERROR_STANDARD);

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		e_ce_h->s_vol[0] = (struct Solver_Volume*) curr;
		struct Error_CE_Data* e_ce_d = constructor_Error_CE_Data(e_ce_h,sim); // destructed

		for (int i = 0; i < 2; ++i) {
			convert_variables(e_ce_d->sol[i],'c','p');
		}
		add_euler_variable_Error_CE_Data('t',e_ce_d,sim);

		for (int i = 0; i < 2; ++i) {
			remove_col_Multiarray_d(DIM+1,e_ce_d->sol[i]); // Remove pressure.
			remove_col_Multiarray_d(0,    e_ce_d->sol[i]); // Remove density.
		}

		increment_sol_L2(e_ce_h,e_ce_d);
		destructor_Error_CE_Data(e_ce_d);

		update_domain_order(e_ce_h);
	}

	struct Error_CE* error_ce = constructor_Error_CE(e_ce_h,sim); // returned
	destructor_Error_CE_Helper(e_ce_h);

	return error_ce;
}

struct Error_CE* constructor_Error_CE_navier_stokes_uvw_zero_surface (const struct Simulation*const sim)
{
	const int n_out = DIM;

	struct Error_CE_Helper* e_ce_h = constructor_Error_CE_Helper(sim,n_out);
	e_ce_h->header_spec = compute_header_spec_navier_stokes_uvw();
	const_cast_i(&e_ce_h->error_type,ERROR_STANDARD);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face*const face = (struct Face*) curr;
		if (!is_face_wall_boundary(face))
			continue;

		const struct Solver_Face*const s_face = (struct Solver_Face*) curr;
		e_ce_h->s_face = s_face;
		e_ce_h->s_vol[0] = (struct Solver_Volume*) face->neigh_info[0].volume;

		struct Boundary_Value_Input bv_i;
		constructor_Boundary_Value_Input_face_s_fcl_interp(&bv_i,s_face,sim); // destructed
		bv_i.g = NULL;

		struct Error_CE_Data e_ce_d;
		e_ce_d.sol[0] = (struct Multiarray_d*) constructor_copy_const_Multiarray_d(bv_i.s); // destructed
		destructor_Boundary_Value_Input(&bv_i);

		e_ce_d.sol[1] = constructor_copy_Multiarray_d(e_ce_d.sol[0]); // destructed
		set_to_value_Multiarray_d(e_ce_d.sol[1],0.0);

		for (int i = 0; i < 2; ++i) {
			remove_col_Multiarray_d(DIM+1,e_ce_d.sol[i]); // Remove pressure.
			remove_col_Multiarray_d(0,    e_ce_d.sol[i]); // Remove density.
		}

		increment_sol_face_L2(e_ce_h,&e_ce_d);
		update_domain_order(e_ce_h);

		for (int i = 0; i < 2; ++i)
			destructor_Multiarray_d(e_ce_d.sol[i]);
	}

	struct Error_CE* error_ce = constructor_Error_CE(e_ce_h,sim); // returned
	destructor_Error_CE_Helper(e_ce_h);

	return error_ce;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const char* compute_header_spec_navier_stokes_uvwt ( )
{
	static char header_spec[STRLEN_MAX];

	int index = sprintf(header_spec,"%-14s","$u$");
	if (DIM >= 2)
		index += sprintf(header_spec+index,"%-14s","$v$");
	if (DIM >= 3)
		index += sprintf(header_spec+index,"%-14s","$w$");
	sprintf(header_spec+index,"%-14s","$T$");

	return header_spec;
}

static const char* compute_header_spec_navier_stokes_uvw ( )
{
	static char header_spec[STRLEN_MAX];

	int index = sprintf(header_spec,"%-14s","$u$");
	if (DIM >= 2)
		index += sprintf(header_spec+index,"%-14s","$v$");
	if (DIM >= 3)
		index += sprintf(header_spec+index,"%-14s","$w$");
	return header_spec;
}
