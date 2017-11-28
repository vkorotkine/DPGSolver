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
 *  \todo clean-up.
 */

#include "boundary_euler.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_test_case.h"

#include "multiarray.h"

#include "boundary.h"
#include "simulation.h"
#include "solution_euler.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

#define NEQ  NEQ_EULER  ///< Number of equations.
#define NVAR NVAR_EULER ///< Number of variables.

/** \brief Compute the normal velocity.
 *  \return See brief. */
static double compute_Vn
	(const double*const n,  ///< Array of unit normal vector components.
	 const double*const uvw ///< Array of velocity components.
	);

/// \brief Set the velocity components from the input values.
static void set_uvw
	(const double*const uvw_i, ///< Input velocity components.
	 double*const*const uvw    ///< Components to be set.
	);

/** \brief Compute the tangential velocity components.
 *  \return See brief. */
static void compute_Vt
	(const double*const n,   ///< Array of unit normal vector components.
	 const double Vn,        ///< Normal velocity.
	 const double*const uvw, ///< Total velocity components.
	 double*const uvw_t      ///< Set to tangential velocity components.
	);

/// \brief Compute the velocity components from the normal and tangential components.
static void compute_uvw
	(const double*const n,     ///< Array of unit normal vector components.
	 const double Vn,          ///< Normal velocity.
	 const double*const uvw_t, ///< Tangential velocity components.
	 double*const*const uvw    ///< Set to total velocity components.
	);

/// \brief Compute the velocity components as the velocity with the same tangential but opposite normal from the input.
static void compute_opposite_normal_uvw
	(const double*const n,     ///< Array of unit normal vector components.
	 const double Vn,          ///< Normal velocity.
	 const double*const uvw_i, ///< Input velocity components.
	 double*const*const uvw    ///< Set to total velocity components.
	);

// Interface functions ********************************************************************************************** //

void constructor_Boundary_Value_euler_riemann
	(struct Boundary_Value* bv, const struct Boundary_Value_Input* bv_i, const struct Solver_Face* face,
	 const struct Simulation* sim)
{
	UNUSED(face);

	const bool* c_m = bv_i->compute_member;
	assert(c_m[0] == true);

	const struct const_Multiarray_d* xyz = bv_i->xyz;
	const struct const_Multiarray_d* sol_l = bv_i->s;
	const struct const_Multiarray_d* sol_r = sim->test_case->constructor_sol(xyz,sim); // destructed

	const struct const_Multiarray_d* normals = bv_i->normals;
	assert(normals->layout == 'R');

	convert_variables((struct Multiarray_d*)sol_l,'c','p');
	convert_variables((struct Multiarray_d*)sol_r,'c','p');

	const ptrdiff_t n_n = xyz->extents[0];

	struct Multiarray_d* sol = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_n,NVAR}); // keep

	const double*const rho_l = get_col_const_Multiarray_d(0,sol_l),
	            *const u_l   = get_col_const_Multiarray_d(1,sol_l),
		      *const v_l   = (DIM > 1 ? get_col_const_Multiarray_d(2,sol_l) : NULL),
		      *const w_l   = (DIM > 2 ? get_col_const_Multiarray_d(3,sol_l) : NULL),
		      *const p_l   = get_col_const_Multiarray_d(NVAR-1,sol_l),

	            *const rho_r = get_col_const_Multiarray_d(0,sol_r),
	            *const u_r   = get_col_const_Multiarray_d(1,sol_r),
		      *const v_r   = (DIM > 1 ? get_col_const_Multiarray_d(2,sol_r) : NULL),
		      *const w_r   = (DIM > 2 ? get_col_const_Multiarray_d(3,sol_r) : NULL),
		      *const p_r   = get_col_const_Multiarray_d(NVAR-1,sol_r);

	double*const rho = get_col_Multiarray_d(0,sol),
	      *const u   = get_col_Multiarray_d(1,sol),
		*const v   = (DIM > 1 ? get_col_Multiarray_d(2,sol) : NULL),
		*const w   = (DIM > 2 ? get_col_Multiarray_d(3,sol) : NULL),
		*const p   = get_col_Multiarray_d(NVAR-1,sol);

	for (int n = 0; n < n_n; n++) {
		const double uvw_l[] = { u_l[n], (DIM > 1 ? v_l[n] : 0.0), (DIM > 2 ? w_l[n] : 0.0), },
		             uvw_r[] = { u_r[n], (DIM > 1 ? v_r[n] : 0.0), (DIM > 2 ? w_r[n] : 0.0), };
		const double c_l = sqrt(GAMMA*p_l[n]/rho_l[n]),
		             c_r = sqrt(GAMMA*p_r[n]/rho_r[n]);

		const double* data_n = get_row_const_Multiarray_d(n,normals);
		const double Vn_l = compute_Vn(data_n,uvw_l),
		             Vn_r = compute_Vn(data_n,uvw_r);

		// Riemann invariants
		const double R_l = Vn_l + 2.0/GM1*c_l, // Outgoing
		             R_r = Vn_r - 2.0/GM1*c_r; // Incoming

		const double Vn = 0.5*(R_l+R_r),
		             c  = 0.25*GM1*(R_l-R_r);

		double* uvw[] = { &u[n], (DIM > 1 ? &v[n] : NULL), (DIM > 2 ? &w[n] : NULL), };
		if (fabs(Vn) >= c) { // Supersonic
			if (Vn < 0.0) { // Inlet
				rho[n] = rho_r[n];
				set_uvw(uvw_r,uvw);
				p[n]   = p_r[n];
			} else { // Outlet
				rho[n] = rho_l[n];
				set_uvw(uvw_l,uvw);
				p[n]   = p_l[n];
			}
		} else { // Subsonic
			double uvw_t[] = { 0.0, 0.0, 0.0, };
			if (Vn < 0.0) { // Inlet
				const double s_r = sqrt(p_r[n]/pow(rho_r[n],GAMMA));
				compute_Vt(data_n,Vn_r,uvw_r,uvw_t);
				rho[n] = pow(1.0/GAMMA*c*c/(s_r*s_r),1.0/GM1);
			} else { // Outlet
				const double s_l = sqrt(p_l[n]/pow(rho_l[n],GAMMA));
				compute_Vt(data_n,Vn_l,uvw_l,uvw_t);
				rho[n] = pow(1.0/GAMMA*c*c/(s_l*s_l),1.0/GM1);
			}
			compute_uvw(data_n,Vn,uvw_t,uvw);

			p[n] = 1.0/GAMMA*c*c*rho[n];
		}
	}
	convert_variables(sol,'p','c');
	bv->s = (struct const_Multiarray_d*)sol;

	convert_variables((struct Multiarray_d*)sol_l,'p','c');
	destructor_const_Multiarray_d(sol_r);

	if (c_m[1] == true) {
		EXIT_ADD_SUPPORT;
	}
	assert(c_m[2] == false);
}

void constructor_Boundary_Value_euler_slipwall
	(struct Boundary_Value* bv, const struct Boundary_Value_Input* bv_i, const struct Solver_Face* face,
	 const struct Simulation* sim)
{
	UNUSED(face);
	UNUSED(sim);
	const bool* c_m = bv_i->compute_member;

	assert(c_m[0] == true);

	const struct const_Multiarray_d* sol_l = bv_i->s;

	const double*const rhou_l = get_col_const_Multiarray_d(1,sol_l),
		      *const rhov_l = (DIM > 1 ? get_col_const_Multiarray_d(2,sol_l) : NULL),
		      *const rhow_l = (DIM > 2 ? get_col_const_Multiarray_d(3,sol_l) : NULL);

	const ptrdiff_t n_n = sol_l->extents[0];
	struct Multiarray_d* sol = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_n,NVAR}); // keep

	double*const rhou = get_col_Multiarray_d(1,sol),
		*const rhov = (DIM > 1 ? get_col_Multiarray_d(2,sol) : NULL),
		*const rhow = (DIM > 2 ? get_col_Multiarray_d(3,sol) : NULL);

	const struct const_Multiarray_d* normals = bv_i->normals;
	assert(normals->layout == 'R');

	for (int n = 0; n < n_n; n++) {
		const double* data_n = get_row_const_Multiarray_d(n,normals);
		const double rhouvw_l[] = { rhou_l[n], (DIM > 1 ? rhov_l[n] : 0.0), (DIM > 2 ? rhow_l[n] : 0.0), };
		const double rhoVn_l = compute_Vn(data_n,rhouvw_l);

		double* rhouvw[] = { &rhou[n], (DIM > 1 ? &rhov[n] : NULL), (DIM > 2 ? &rhow[n] : NULL), };
		compute_opposite_normal_uvw(data_n,rhoVn_l,rhouvw_l,rhouvw);
	}
	bv->s = (struct const_Multiarray_d*)sol;

	if (c_m[1] == true) {
		struct Multiarray_d* ds_ds = constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){n_n,NVAR,NVAR}); // keep

		double *ds_ds_ptr[NEQ*NEQ];
		for (int eq  = 0; eq  < NEQ;  eq++)  {
		for (int var = 0; var < NVAR; var++) {
			ds_ds_ptr[eq*NVAR+var] = &ds_ds->data[(eq*NVAR+var)*n_n];
		}}
		const double* n_ptr = get_row_const_Multiarray_d(0,normals);

		if (DIM == 3) {
			for (int n = 0; n < n_n; n++) {
				const double n1 = *n_ptr++,
				             n2 = *n_ptr++,
				             n3 = *n_ptr++;

				int Indds_ds = 0;

				// *** var 1 ***
				*ds_ds_ptr[Indds_ds++] = 1.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;

				// *** var 2 ***
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 1.0 - 2.0*n1*n1;
				*ds_ds_ptr[Indds_ds++] =     - 2.0*n2*n1;
				*ds_ds_ptr[Indds_ds++] =     - 2.0*n3*n1;
				*ds_ds_ptr[Indds_ds++] = 0.0;

				// *** var 3 ***
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] =     - 2.0*n1*n2;
				*ds_ds_ptr[Indds_ds++] = 1.0 - 2.0*n2*n2;
				*ds_ds_ptr[Indds_ds++] =     - 2.0*n3*n2;
				*ds_ds_ptr[Indds_ds++] = 0.0;

				// *** var 4 ***
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] =     - 2.0*n1*n3;
				*ds_ds_ptr[Indds_ds++] =     - 2.0*n2*n3;
				*ds_ds_ptr[Indds_ds++] = 1.0 - 2.0*n3*n3;
				*ds_ds_ptr[Indds_ds++] = 0.0;

				// *** var 5 ***
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 1.0;

				for (int i = 0, iMax = NEQ*NVAR; i < iMax; i++)
					ds_ds_ptr[i]++;
			}
		} else if (DIM == 2) {
			for (int n = 0; n < n_n; n++) {
				const double n1 = *n_ptr++,
				             n2 = *n_ptr++;

				int Indds_ds = 0;

				// *** var 1 ***
				*ds_ds_ptr[Indds_ds++] = 1.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;

				// *** var 2 ***
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 1.0 - 2.0*n1*n1;
				*ds_ds_ptr[Indds_ds++] =     - 2.0*n2*n1;
				*ds_ds_ptr[Indds_ds++] = 0.0;

				// *** var 3 ***
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] =     - 2.0*n1*n2;
				*ds_ds_ptr[Indds_ds++] = 1.0 - 2.0*n2*n2;
				*ds_ds_ptr[Indds_ds++] = 0.0;

				// *** var 4 ***
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 1.0;

				for (int i = 0, iMax = NEQ*NVAR; i < iMax; i++)
					ds_ds_ptr[i]++;
			}
		} else if (DIM == 1) {
			for (int n = 0; n < n_n; n++) {
				const double n1 = *n_ptr++;

				int Indds_ds = 0;

				// *** var 1 ***
				*ds_ds_ptr[Indds_ds++] = 1.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;

				// *** var 2 ***
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 1.0 - 2.0*n1*n1;
				*ds_ds_ptr[Indds_ds++] = 0.0;

				// *** var 3 ***
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 0.0;
				*ds_ds_ptr[Indds_ds++] = 1.0;

				for (int i = 0, iMax = NEQ*NVAR; i < iMax; i++)
					ds_ds_ptr[i]++;
			}
		}
	}
	assert(c_m[2] == false);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static double compute_Vn (const double*const n, const double*const uvw)
{
	double Vn = 0.0;
	for (int d = 0; d < DIM; ++d)
		Vn += n[d]*uvw[d];
	return Vn;
}

static void set_uvw (const double*const uvw_i, double*const*const uvw)
{
	for (int d = 0; d < DIM; ++d)
		*uvw[d] = uvw_i[d];
}

static void compute_Vt
	(const double*const n, const double Vn, const double*const uvw, double*const uvw_t)
{
	for (int d = 0; d < DIM; ++d)
		uvw_t[d] = uvw[d] - Vn*n[d];
}

static void compute_uvw
	(const double*const n, const double Vn, const double*const uvw_t, double*const*const uvw)
{
	for (int d = 0; d < DIM; ++d)
		*uvw[d] = Vn*n[d] + uvw_t[d];
}

static void compute_opposite_normal_uvw
	(const double*const n, const double Vn, const double*const uvw_i, double*const*const uvw)
{
	for (int d = 0; d < DIM; ++d)
		*uvw[d] = uvw_i[d] - 2.0*Vn*n[d];
}
