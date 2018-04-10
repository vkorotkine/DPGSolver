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

#include "compute_error_euler.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_error.h"

#include "face_solver.h"
#include "volume.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "boundary.h"
#include "compute_error.h"
#include "const_cast.h"
#include "element.h"
#include "intrusive.h"
#include "simulation.h"
#include "solution_euler.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a statically allocated `char*` holding the specific header for all of the Euler variables.
 *  \return See brief. */
static const char* compute_header_spec_euler_all
	();

/** \brief Return a statically allocated `char*` holding the specific header for all of the Euler variables and
 *         their residuals.
 *  \return See brief. */
static const char* compute_header_spec_euler_all_p_rhs
	();

/** \brief Return a statically allocated `char*` holding the specific header for the entropy variable.
 *  \return See brief. */
static const char* compute_header_spec_euler_entropy
	();

/** \brief Version of \ref constructor_Error_CE_fptr checking the error of pressure drag and lift coefficients where
 *         specified functionals may be removed.
 *  \return See brief. */
struct Error_CE* constructor_Error_CE_functionals__cd_cl_general
	(const struct Simulation*const sim, ///< \ref Simulation.
	 const int remove_cd_cl             ///< Flag for which functional to remove (if any).
	);

// Interface functions ********************************************************************************************** //

struct Error_CE* constructor_Error_CE_euler_all (const struct Simulation* sim)
{
	const int n_out = DIM+2+1;

	struct Error_CE_Helper* e_ce_h = constructor_Error_CE_Helper(sim,n_out);
	e_ce_h->header_spec = compute_header_spec_euler_all();
	const_cast_i(&e_ce_h->error_type,ERROR_STANDARD);

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		e_ce_h->s_vol[0] = (struct Solver_Volume*) curr;
		struct Error_CE_Data* e_ce_d = constructor_Error_CE_Data(e_ce_h,sim); // destructed

		for (int i = 0; i < 2; ++i)
			convert_variables(e_ce_d->sol[i],'c','p');
		add_euler_variable_Error_CE_Data('s',e_ce_d,sim);

		increment_sol_L2(e_ce_h,e_ce_d);
		destructor_Error_CE_Data(e_ce_d);

		update_domain_order(e_ce_h);
	}

	struct Error_CE* error_ce = constructor_Error_CE(e_ce_h,sim); // returned
	destructor_Error_CE_Helper(e_ce_h);

	return error_ce;
}

struct Error_CE* constructor_Error_CE_euler_all_p_rhs (const struct Simulation* sim)
{
	const int n_out = 2*(DIM+2)+1;

	struct Error_CE_Helper* e_ce_h = constructor_Error_CE_Helper(sim,n_out);
	e_ce_h->header_spec = compute_header_spec_euler_all_p_rhs();
	const_cast_i(&e_ce_h->error_type,ERROR_STANDARD);

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		e_ce_h->s_vol[0] = (struct Solver_Volume*) curr;
		struct Error_CE_Data* e_ce_d = constructor_Error_CE_Data(e_ce_h,sim); // destructed

		for (int i = 0; i < 2; ++i)
			convert_variables(e_ce_d->sol[i],'c','p');
		add_euler_variable_Error_CE_Data('s',e_ce_d,sim);
		add_rhs_Error_CE_Data(e_ce_d,sim);

		increment_sol_L2(e_ce_h,e_ce_d);
		destructor_Error_CE_Data(e_ce_d);

		update_domain_order(e_ce_h);
	}

	struct Error_CE* error_ce = constructor_Error_CE(e_ce_h,sim); // returned
	destructor_Error_CE_Helper(e_ce_h);

	return error_ce;
}

struct Error_CE* constructor_Error_CE_euler_entropy (const struct Simulation* sim)
{
	const int n_out = 1;

	struct Error_CE_Helper* e_ce_h = constructor_Error_CE_Helper(sim,n_out);
	e_ce_h->header_spec = compute_header_spec_euler_entropy();
	const_cast_i(&e_ce_h->error_type,ERROR_STANDARD);

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		e_ce_h->s_vol[0] = (struct Solver_Volume*) curr;
		struct Error_CE_Data* e_ce_d = constructor_Error_CE_Data(e_ce_h,sim); // destructed

		for (int i = 0; i < 2; ++i)
			convert_variables(e_ce_d->sol[i],'c','p');
		add_euler_variable_Error_CE_Data('s',e_ce_d,sim);
		for (int i = 0; i < 2; ++i) {
			for (int vr = 0; vr < DIM+2; ++vr) // Remove all except entropy.
				remove_col_Multiarray_d(0,e_ce_d->sol[i]);
		}

		increment_sol_L2(e_ce_h,e_ce_d);
		destructor_Error_CE_Data(e_ce_d);

		update_domain_order(e_ce_h);
	}

	struct Error_CE* error_ce = constructor_Error_CE(e_ce_h,sim); // returned
	destructor_Error_CE_Helper(e_ce_h);

	return error_ce;
}

void add_euler_variable_Error_CE_Data
	(const char var_type, struct Error_CE_Data*const e_ce_d, const struct Simulation*const sim)
{
	UNUSED(sim);
	const ptrdiff_t ext_0 = e_ce_d->sol[0]->extents[0];
	struct Multiarray_d*const var = constructor_move_Multiarray_d_d('C',2,(ptrdiff_t[]){ext_0,1},false,NULL); // dest.
	for (int i = 0; i < 2; ++i) {
		struct Multiarray_d* sol = e_ce_d->sol[i];

		const ptrdiff_t ext_1_new = sol->extents[1]+1;
		resize_Multiarray_d(sol,sol->order,(ptrdiff_t[]){ext_0,ext_1_new});

		sol->extents[1] = DIM+2;
		var->data = get_col_Multiarray_d(ext_1_new-1,sol);
		switch (var_type) {
			case 's': compute_entropy(var,(const struct const_Multiarray_d*)sol,'p'); break;
			case 't': compute_temperature(var,(const struct const_Multiarray_d*)sol,'p'); break;
			default:  EXIT_ERROR("Unsupported: %c.",var_type); break;
		}
		sol->extents[1] = ext_1_new;
	}
	destructor_Multiarray_d(var);
}

struct Error_CE* constructor_Error_CE_functionals__cd_cl (const struct Simulation*const sim)
{
	return constructor_Error_CE_functionals__cd_cl_general(sim,REMOVE_CD_CL__NONE);
}

struct Error_CE* constructor_Error_CE_functionals__cl (const struct Simulation*const sim)
{
	return constructor_Error_CE_functionals__cd_cl_general(sim,REMOVE_CD_CL__D);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Return a statically allocated `char*` holding the specific header for the drag/lift functionals.
 *  \return See brief. */
static const char* compute_header_spec_cd_cl
	();

/** \brief Return a statically allocated `char*` holding the specific header for the lift functional.
 *  \return See brief. */
static const char* compute_header_spec_cl
	();

static const char* compute_header_spec_euler_all ( )
{
	static char header_spec[STRLEN_MAX];

	int index = sprintf(header_spec,"%-14s%-14s","$\\rho$","$u$");
	if (DIM >= 2)
		index += sprintf(header_spec+index,"%-14s","$v$");
	if (DIM >= 3)
		index += sprintf(header_spec+index,"%-14s","$w$");
	sprintf(header_spec+index,"%-14s%-14s","$p$","$s$");

	return header_spec;
}

static const char* compute_header_spec_euler_all_p_rhs ( )
{
	static char header_spec[STRLEN_MAX];

	int index = sprintf(header_spec,"%-14s%-14s","$\\rho$","$u$");
	if (DIM >= 2)
		index += sprintf(header_spec+index,"%-14s","$v$");
	if (DIM >= 3)
		index += sprintf(header_spec+index,"%-14s","$w$");
	index += sprintf(header_spec+index,"%-14s%-14s","$p$","$s$");

	index += sprintf(header_spec+index,"%-14s%-14s","$\\rho_{res}$","$u_{res}$");
	if (DIM >= 2)
		index += sprintf(header_spec+index,"%-14s","$v_{res}$");
	if (DIM >= 3)
		index += sprintf(header_spec+index,"%-14s","$w_{res}$");
	sprintf(header_spec+index,"%-14s","$E_{res}$");

	return header_spec;
}

static const char* compute_header_spec_euler_entropy ( )
{
	static char header_spec[STRLEN_MAX];
	sprintf(header_spec,"%-14s","$s$");
	return header_spec;
}

struct Error_CE* constructor_Error_CE_functionals__cd_cl_general
	(const struct Simulation*const sim, const int remove_cd_cl)
{
	const int n_out = (remove_cd_cl ? 1 : 2);

	struct Error_CE_Helper* e_ce_h = constructor_Error_CE_Helper(sim,n_out);
	switch (remove_cd_cl) {
		case REMOVE_CD_CL__NONE: e_ce_h->header_spec = compute_header_spec_cd_cl(); break;
		case REMOVE_CD_CL__D:    e_ce_h->header_spec = compute_header_spec_cl();    break;
		default: EXIT_ERROR("Unsupported: %d\n",remove_cd_cl); break;
	}
	const_cast_i(&e_ce_h->error_type,ERROR_FUNCTIONAL);

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

		const struct const_Multiarray_d*const xyz_fc = s_face->xyz_fc;
		const ptrdiff_t n_n = bv_i.xyz->extents[0];

		struct Error_CE_Data e_ce_d;
		e_ce_d.sol[0] = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_n,2}); // destructed
		compute_cd_cl_values(e_ce_d.sol[0],bv_i.s,'c',bv_i.normals);
		destructor_Boundary_Value_Input(&bv_i);

		e_ce_d.sol[1] = (struct Multiarray_d*)constructor_const_functionals_cd_cl_zero(xyz_fc,sim); // destructed

		for (int i = 0; i < 2; ++i) {
			switch (remove_cd_cl) {
			case REMOVE_CD_CL__NONE:
				break; // do nothing.
			case REMOVE_CD_CL__D:
				remove_col_Multiarray_d(0,e_ce_d.sol[i]); // Remove drag.
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",remove_cd_cl);
				break;
			}
		}

		increment_sol_integrated_face(e_ce_h,&e_ce_d);
		update_domain_order(e_ce_h);

		for (int i = 0; i < 2; ++i)
			destructor_Multiarray_d(e_ce_d.sol[i]);
	}

	struct Error_CE* error_ce = constructor_Error_CE(e_ce_h,sim); // returned
	destructor_Error_CE_Helper(e_ce_h);

	return error_ce;
}

// Level 1 ********************************************************************************************************** //

static const char* compute_header_spec_cd_cl ( )
{
	static char header_spec[STRLEN_MAX];
	sprintf(header_spec,"%-14s%-14s","$C_D$","$C_L$");
	return header_spec;
}

static const char* compute_header_spec_cl ( )
{
	static char header_spec[STRLEN_MAX];
	sprintf(header_spec,"%-14s","$C_L$");
	return header_spec;
}
