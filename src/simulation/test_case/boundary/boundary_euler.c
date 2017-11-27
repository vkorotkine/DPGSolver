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

/** \brief Pointer to functions computing the normal velocity.
 *
 *  \param i      The index of the current dof.
 *  \param data_n The data for the unit normal vector for the current dof.
 *  \param u      The pointer to the u-velocity data.
 *  \param v      The pointer to the v-velocity data (NULL for d < 2).
 *  \param w      The pointer to the w-velocity data (NULL for d < 3).
 */
typedef double (*compute_Vn_fptr)
	(const int i,
	 const double*const data_n,
	 const double*const u,
	 const double*const v,
	 const double*const w
	);

/** \brief Pointer to functions computing the tangential velocity components.
 *
 *  \param i      The index of the current dof.
 *  \param data_n The data for the unit normal vector for the current dof.
 *  \param u      The pointer to the u-velocity data.
 *  \param v      The pointer to the v-velocity data (NULL for d < 2).
 *  \param w      The pointer to the w-velocity data (NULL for d < 3).
 *  \param Vn     The normal velocity.
 *  \param ut     The pointer to the tangential u-velocity data.
 *  \param vt     The pointer to the tangential v-velocity data.
 *  \param wt     The pointer to the tangential w-velocity data.
 */
typedef void (*compute_Vt_fptr)
	(const int i,
	 const double*const data_n,
	 const double*const u,
	 const double*const v,
	 const double*const w,
	 const double Vn,
	 double*const ut,
	 double*const vt,
	 double*const wt
	);

/** \brief Pointer to functions computing uvw values from normal and tangential components.
 *
 *  \param i      The index of the current dof.
 *  \param data_n The data for the unit normal vector for the current dof.
 *  \param Vn     The normal velocity.
 *  \param ut     The tangential u-velocity data.
 *  \param vt     The tangential v-velocity data.
 *  \param wt     The tangential w-velocity data.
 *  \param u      The pointer to the u-velocity data.
 *  \param v      The pointer to the v-velocity data (NULL for d < 2).
 *  \param w      The pointer to the w-velocity data (NULL for d < 3).
 */
typedef void (*compute_uvw_fptr)
	(const int i,
	 const double*const data_n,
	 const double Vn,
	 const double ut,
	 const double vt,
	 const double wt,
	 double*const u,
	 double*const v,
	 double*const w
	);

/** \brief Pointer to functions setting uvw values from input values.
 *
 *  \param i      The index of the current dof.
 *  \param u_i    The pointer to input the u-velocity data.
 *  \param v_i    The pointer to input the v-velocity data (NULL for d < 2).
 *  \param w_i    The pointer to input the w-velocity data (NULL for d < 3).
 *  \param u      The pointer to the u-velocity data.
 *  \param v      The pointer to the v-velocity data (NULL for d < 2).
 *  \param w      The pointer to the w-velocity data (NULL for d < 3).
 */
typedef void (*set_uvw_fptr)
	(const int i,
	 const double*const u_i,
	 const double*const v_i,
	 const double*const w_i,
	 double*const u,
	 double*const v,
	 double*const w
	);

/// \brief 1D version of \ref compute_Vn_fptr.
static double compute_Vn_1d
	(const int i,               ///< See brief.
	 const double*const data_n, ///< See brief.
	 const double*const u,      ///< See brief.
	 const double*const v,      ///< See brief.
	 const double*const w       ///< See brief.
	);

/// \brief 2D version of \ref compute_Vn_fptr.
static double compute_Vn_2d
	(const int i,               ///< See brief.
	 const double*const data_n, ///< See brief.
	 const double*const u,      ///< See brief.
	 const double*const v,      ///< See brief.
	 const double*const w       ///< See brief.
	);

/// \brief 3D version of \ref compute_Vn_fptr.
static double compute_Vn_3d
	(const int i,               ///< See brief.
	 const double*const data_n, ///< See brief.
	 const double*const u,      ///< See brief.
	 const double*const v,      ///< See brief.
	 const double*const w       ///< See brief.
	);

/// \brief 1D version of \ref compute_Vt_fptr.
static void compute_Vt_1d
	(const int i,               ///< See brief.
	 const double*const data_n, ///< See brief.
	 const double*const u,      ///< See brief.
	 const double*const v,      ///< See brief.
	 const double*const w,      ///< See brief.
	 const double Vn,           ///< See brief.
	 double*const ut,           ///< See brief.
	 double*const vt,           ///< See brief.
	 double*const wt            ///< See brief.
	);

/// \brief 2D version of \ref compute_Vt_fptr.
static void compute_Vt_2d
	(const int i,               ///< See brief.
	 const double*const data_n, ///< See brief.
	 const double*const u,      ///< See brief.
	 const double*const v,      ///< See brief.
	 const double*const w,      ///< See brief.
	 const double Vn,           ///< See brief.
	 double*const ut,           ///< See brief.
	 double*const vt,           ///< See brief.
	 double*const wt            ///< See brief.
	);

/// \brief 3D version of \ref compute_Vt_fptr.
static void compute_Vt_3d
	(const int i,               ///< See brief.
	 const double*const data_n, ///< See brief.
	 const double*const u,      ///< See brief.
	 const double*const v,      ///< See brief.
	 const double*const w,      ///< See brief.
	 const double Vn,           ///< See brief.
	 double*const ut,           ///< See brief.
	 double*const vt,           ///< See brief.
	 double*const wt            ///< See brief.
	);

/// \brief 1D version of \ref compute_uvw_fptr.
static void compute_uvw_1d
	(const int i,               ///< See brief.
	 const double*const data_n, ///< See brief.
	 const double Vn,           ///< See brief.
	 const double ut,           ///< See brief.
	 const double vt,           ///< See brief.
	 const double wt,           ///< See brief.
	 double*const u,            ///< See brief.
	 double*const v,            ///< See brief.
	 double*const w             ///< See brief.
	);

/// \brief 2D version of \ref compute_uvw_fptr.
static void compute_uvw_2d
	(const int i,               ///< See brief.
	 const double*const data_n, ///< See brief.
	 const double Vn,           ///< See brief.
	 const double ut,           ///< See brief.
	 const double vt,           ///< See brief.
	 const double wt,           ///< See brief.
	 double*const u,            ///< See brief.
	 double*const v,            ///< See brief.
	 double*const w             ///< See brief.
	);

/// \brief 3D version of \ref compute_uvw_fptr.
static void compute_uvw_3d
	(const int i,               ///< See brief.
	 const double*const data_n, ///< See brief.
	 const double Vn,           ///< See brief.
	 const double ut,           ///< See brief.
	 const double vt,           ///< See brief.
	 const double wt,           ///< See brief.
	 double*const u,            ///< See brief.
	 double*const v,            ///< See brief.
	 double*const w             ///< See brief.
	);

/// \brief 1D version of \ref set_uvw_fptr.
static void set_uvw_1d
	(const int i,            ///< See brief.
	 const double*const u_i, ///< See brief.
	 const double*const v_i, ///< See brief.
	 const double*const w_i, ///< See brief.
	 double*const u,         ///< See brief.
	 double*const v,         ///< See brief.
	 double*const w          ///< See brief.
	);

/// \brief 2D version of \ref set_uvw_fptr.
static void set_uvw_2d
	(const int i,            ///< See brief.
	 const double*const u_i, ///< See brief.
	 const double*const v_i, ///< See brief.
	 const double*const w_i, ///< See brief.
	 double*const u,         ///< See brief.
	 double*const v,         ///< See brief.
	 double*const w          ///< See brief.
	);

/// \brief 3D version of \ref set_uvw_fptr.
static void set_uvw_3d
	(const int i,            ///< See brief.
	 const double*const u_i, ///< See brief.
	 const double*const v_i, ///< See brief.
	 const double*const w_i, ///< See brief.
	 double*const u,         ///< See brief.
	 double*const v,         ///< See brief.
	 double*const w          ///< See brief.
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

	// set function pointers.
	static bool fptrs_set = false;
	static compute_Vn_fptr compute_Vn   = NULL;
	static compute_Vt_fptr compute_Vt   = NULL;
	static compute_uvw_fptr compute_uvw = NULL;
	static set_uvw_fptr set_uvw         = NULL;
	if (!fptrs_set) {
		fptrs_set = true;
		switch (sim->d) {
		case 1:
			compute_Vn  = compute_Vn_1d;
			compute_Vt  = compute_Vt_1d;
			compute_uvw = compute_uvw_1d;
			set_uvw     = set_uvw_1d;
			break;
		case 2:
			compute_Vn  = compute_Vn_2d;
			compute_Vt  = compute_Vt_2d;
			compute_uvw = compute_uvw_2d;
			set_uvw     = set_uvw_2d;
			break;
		case 3:
			compute_Vn  = compute_Vn_3d;
			compute_Vt  = compute_Vt_3d;
			compute_uvw = compute_uvw_3d;
			set_uvw     = set_uvw_3d;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	}

	convert_variables((struct Multiarray_d*)sol_l,'c','p');
	convert_variables((struct Multiarray_d*)sol_r,'c','p');

	const int d     = sim->d,
	          n_var = sim->test_case->n_var;
	const ptrdiff_t n_n = xyz->extents[0];

	struct Multiarray_d* sol = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_n,NVAR3D}); // keep

	const double*const rho_l = get_col_const_Multiarray_d(0,sol_l),
	            *const u_l   = get_col_const_Multiarray_d(1,sol_l),
		      *const v_l   = (d > 1 ? get_col_const_Multiarray_d(2,sol_l) : NULL),
		      *const w_l   = (d > 2 ? get_col_const_Multiarray_d(3,sol_l) : NULL),
		      *const p_l   = get_col_const_Multiarray_d(n_var-1,sol_l),

	            *const rho_r = get_col_const_Multiarray_d(0,sol_r),
	            *const u_r   = get_col_const_Multiarray_d(1,sol_r),
		      *const v_r   = (d > 1 ? get_col_const_Multiarray_d(2,sol_r) : NULL),
		      *const w_r   = (d > 2 ? get_col_const_Multiarray_d(3,sol_r) : NULL),
		      *const p_r   = get_col_const_Multiarray_d(n_var-1,sol_r);

	double*const rho = get_col_Multiarray_d(0,sol),
	      *const u   = get_col_Multiarray_d(1,sol),
		*const v   = (d > 1 ? get_col_Multiarray_d(2,sol) : NULL),
		*const w   = (d > 2 ? get_col_Multiarray_d(3,sol) : NULL),
		*const p   = get_col_Multiarray_d(n_var-1,sol);

	for (int n = 0; n < n_n; n++) {
		const double c_l = sqrt(GAMMA*p_l[n]/rho_l[n]),
		             c_r = sqrt(GAMMA*p_r[n]/rho_r[n]);

		const double* data_n = get_row_const_Multiarray_d(n,normals);
		const double Vn_l = compute_Vn(n,data_n,u_l,v_l,w_l),
		             Vn_r = compute_Vn(n,data_n,u_r,v_r,w_r);

		// Riemann invariants
		const double R_l = Vn_l + 2.0/GM1*c_l, // Outgoing
		             R_r = Vn_r - 2.0/GM1*c_r; // Incoming

		const double Vn = 0.5*(R_l+R_r),
		             c  = 0.25*GM1*(R_l-R_r);

		if (fabs(Vn) >= c) { // Supersonic
			if (Vn < 0.0) { // Inlet
				rho[n] = rho_r[n];
				set_uvw(n,u_r,v_r,w_r,u,v,w);
				p[n]   = p_r[n];
			} else { // Outlet
				rho[n] = rho_l[n];
				set_uvw(n,u_l,v_l,w_l,u,v,w);
				p[n]   = p_l[n];
			}
		} else { // Subsonic
			double ut = 0.0, vt = 0.0, wt = 0.0;
			if (Vn < 0.0) { // Inlet
				const double s_r = sqrt(p_r[n]/pow(rho_r[n],GAMMA));
				compute_Vt(n,data_n,u_r,v_r,w_r,Vn_r,&ut,&vt,&wt);
				rho[n] = pow(1.0/GAMMA*c*c/(s_r*s_r),1.0/GM1);
			} else { // Outlet
				const double s_l = sqrt(p_l[n]/pow(rho_l[n],GAMMA));
				compute_Vt(n,data_n,u_l,v_l,w_l,Vn_l,&ut,&vt,&wt);
				rho[n] = pow(1.0/GAMMA*c*c/(s_l*s_l),1.0/GM1);
			}
			compute_uvw(n,data_n,Vn,ut,vt,wt,u,v,w);

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

	const double*const rho_l  = get_col_const_Multiarray_d(0,sol_l),
	            *const rhou_l = get_col_const_Multiarray_d(1,sol_l),
		      *const rhov_l = (d > 1 ? get_col_const_Multiarray_d(2,sol_l) : NULL),
		      *const rhow_l = (d > 2 ? get_col_const_Multiarray_d(3,sol_l) : NULL),
		      *const E_l    = get_col_const_Multiarray_d(n_var-1,sol_l),

	struct Multiarray_d* sol = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_n,NVAR3D}); // keep

	double*const rho  = get_col_Multiarray_d(0,sol),
	      *const rhou = get_col_Multiarray_d(1,sol),
		*const rhov = (d > 1 ? get_col_Multiarray_d(2,sol) : NULL),
		*const rhow = (d > 2 ? get_col_Multiarray_d(3,sol) : NULL),
		*const E    = get_col_Multiarray_d(n_var-1,sol);

	const int d = sim->d;
	const ptrdiff_t n_n = xyz->extents[0];
	for (n = 0; n < n_n; n++) {
		const double* data_n = get_row_const_Multiarray_d(n,normals);
		const double rhoVn_l = compute_Vn(n,data_n,rhou_l,rhov_l,rhow_l);

		compute_opposite_normal_uvw(n,data_n,rhoVn_l,rhou_l,rhov_l,rhow_l,rhou,rhov,rhow);
//		rhou[n] = rhou_l[n]-2.0*rhoVn_l*data_n[0];
//		rhov[n] = rhov_l[n]-2.0*rhoVn_l*data_n[1];
//		rhow[n] = rhow_l[n]-2.0*rhoVn_l*data_n[2];
	}
EXIT_ADD_SUPPORT;
UNUSED(bv);

	if (c_m[1] == true) {
		EXIT_ADD_SUPPORT;
	}

	assert(c_m[2] == false);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static double compute_Vn_1d
	(const int i, const double*const data_n, const double*const u, const double*const v, const double*const w)
{
	assert(v == NULL);
	assert(w == NULL);
	return data_n[0]*u[i];
}

static double compute_Vn_2d
	(const int i, const double*const data_n, const double*const u, const double*const v, const double*const w)
{
	assert(w == NULL);
	return data_n[0]*u[i]+data_n[1]*v[i];
}

static double compute_Vn_3d
	(const int i, const double*const data_n, const double*const u, const double*const v, const double*const w)
{
	return data_n[0]*u[i]+data_n[1]*v[i]+data_n[2]*w[i];
}

static void compute_Vt_1d
	(const int i, const double*const data_n, const double*const u, const double*const v, const double*const w,
	 const double Vn, double*const ut, double*const vt, double*const wt)
{
	assert(v == NULL);
	assert(w == NULL);

	*ut = u[i] - Vn*data_n[0];
	UNUSED(vt);
	UNUSED(wt);
}

static void compute_Vt_2d
	(const int i, const double*const data_n, const double*const u, const double*const v, const double*const w,
	 const double Vn, double*const ut, double*const vt, double*const wt)
{
	assert(w == NULL);

	*ut = u[i] - Vn*data_n[0];
	*vt = v[i] - Vn*data_n[1];
	UNUSED(wt);
}

static void compute_Vt_3d
	(const int i, const double*const data_n, const double*const u, const double*const v, const double*const w,
	 const double Vn, double*const ut, double*const vt, double*const wt)
{
	*ut = u[i] - Vn*data_n[0];
	*vt = v[i] - Vn*data_n[1];
	*wt = w[i] - Vn*data_n[2];
}

static void compute_uvw_1d
	(const int i, const double*const data_n, const double Vn, const double ut, const double vt, const double wt,
	 double*const u, double*const v, double*const w)
{
	assert(v == NULL);
	assert(w == NULL);
	u[i] = Vn*data_n[0] + ut;
	UNUSED(vt);
	UNUSED(wt);
}

static void compute_uvw_2d
	(const int i, const double*const data_n, const double Vn, const double ut, const double vt, const double wt,
	 double*const u, double*const v, double*const w)
{
	assert(w == NULL);
	u[i] = Vn*data_n[0] + ut;
	v[i] = Vn*data_n[1] + vt;
	UNUSED(wt);
}

static void compute_uvw_3d
	(const int i, const double*const data_n, const double Vn, const double ut, const double vt, const double wt,
	 double*const u, double*const v, double*const w)
{
	u[i] = Vn*data_n[0] + ut;
	v[i] = Vn*data_n[1] + vt;
	w[i] = Vn*data_n[2] + wt;
}

static void set_uvw_1d
	(const int i, const double*const u_i, const double*const v_i, const double*const w_i,
	 double*const u, double*const v, double*const w)
{
	assert(v == NULL);
	assert(w == NULL);
	assert(v_i == NULL);
	assert(w_i == NULL);
	u[i] = u_i[i];
}

static void set_uvw_2d
	(const int i, const double*const u_i, const double*const v_i, const double*const w_i,
	 double*const u, double*const v, double*const w)
{
	assert(w == NULL);
	assert(w_i == NULL);
	u[i] = u_i[i];
	v[i] = v_i[i];
}

static void set_uvw_3d
	(const int i, const double*const u_i, const double*const v_i, const double*const w_i,
	 double*const u, double*const v, double*const w)
{
	u[i] = u_i[i];
	v[i] = v_i[i];
	w[i] = w_i[i];
}

static void compute_opposite_normal_uvw_2d
	(const int i, const double*const data_n, const double Vn, const double*const u_i, const double*const v_i,
	 const double*const w_i, double*const u, double*const v, double*const w)
{
	assert(v_i == NULL);
	assert(w_i == NULL);
	assert(v == NULL);
	assert(w == NULL);
	u[n] = u_i[n] - 2.0*Vn*data_n[0];
}

static void compute_opposite_normal_uvw_2d
	(const int i, const double*const data_n, const double Vn, const double*const u_i, const double*const v_i,
	 const double*const w_i, double*const u, double*const v, double*const w)
{
	assert(w_i == NULL);
	assert(w == NULL);
	u[n] = u_i[n] - 2.0*Vn*data_n[0];
	v[n] = v_i[n] - 2.0*Vn*data_n[1];
}

static void compute_opposite_normal_uvw_3d
	(const int i, const double*const data_n, const double Vn, const double*const u_i, const double*const v_i,
	 const double*const w_i, double*const u, double*const v, double*const w)
{
	u[n] = u_i[n] - 2.0*Vn*data_n[0];
	v[n] = v_i[n] - 2.0*Vn*data_n[1];
	w[n] = w_i[n] - 2.0*Vn*data_n[2];
}
