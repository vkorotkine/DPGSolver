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
 *  \brief Provides the templated Euler boundary condition functions.
 *  \todo clean-up.
 */

#include <assert.h>
#include <math.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_physics.h"


#include "def_templates_boundary.h"
#include "boundary_pde_T.c"

#include "def_templates_multiarray.h"

#include "def_templates_face_solver.h"

#include "def_templates_math_functions.h"
#include "def_templates_solution.h"

// Static function declarations ************************************************************************************* //

#define NEQ  NEQ_EULER  ///< Number of equations.
#define NVAR NVAR_EULER ///< Number of variables.

/** \brief Compute the normal velocity.
 *  \return See brief. */
static Type compute_Vn
	(const double*const n, ///< Array of unit normal vector components.
	 const Type*const uvw  ///< Array of velocity components.
	);

/// \brief Set the velocity components from the input values.
static void set_uvw
	(const Type*const uvw_i, ///< Input velocity components.
	 Type*const*const uvw    ///< Components to be set.
	);

/** \brief Compute the tangential velocity components.
 *  \return See brief. */
static void compute_Vt
	(const double*const n, ///< Array of unit normal vector components.
	 const Type Vn,        ///< Normal velocity.
	 const Type*const uvw, ///< Total velocity components.
	 Type*const uvw_t      ///< Set to tangential velocity components.
	);

/// \brief Compute the velocity components from the normal and tangential components.
static void compute_uvw
	(const double*const n,   ///< Array of unit normal vector components.
	 const Type Vn,          ///< Normal velocity.
	 const Type*const uvw_t, ///< Tangential velocity components.
	 Type*const*const uvw    ///< Set to total velocity components.
	);

/// \brief Compute the velocity components as the velocity with the same tangential but opposite normal from the input.
static void compute_opposite_normal_uvw
	(const double*const n,   ///< Array of unit normal vector components.
	 const Type Vn,          ///< Normal velocity.
	 const Type*const uvw_i, ///< Input velocity components.
	 Type*const*const uvw    ///< Set to total velocity components.
	);

// Interface functions ********************************************************************************************** //

void constructor_Boundary_Value_T_euler_riemann
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
	 const struct Simulation* sim)
{
	UNUSED(face);
	const bool* c_m = bv_i->compute_member;
	assert(c_m[0] == true);

	const struct const_Multiarray_d* xyz = bv_i->xyz;
	const struct const_Multiarray_T* sol_l = bv_i->s;
	const struct const_Multiarray_T* sol_r = constructor_sol_bv(xyz,sim); // destructed

	const struct const_Multiarray_d* normals = bv_i->normals;
	assert(normals->layout == 'R');

	convert_variables_T((struct Multiarray_T*)sol_l,'c','p');
	convert_variables_T((struct Multiarray_T*)sol_r,'c','p');

	const ptrdiff_t n_n = xyz->extents[0];

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,NVAR}); // keep

	const Type*const rho_l = get_col_const_Multiarray_T(0,sol_l),
	          *const u_l   = get_col_const_Multiarray_T(1,sol_l),
	          *const v_l   = (DIM > 1 ? get_col_const_Multiarray_T(2,sol_l) : NULL),
	          *const w_l   = (DIM > 2 ? get_col_const_Multiarray_T(3,sol_l) : NULL),
	          *const p_l   = get_col_const_Multiarray_T(NVAR-1,sol_l),

	          *const rho_r = get_col_const_Multiarray_T(0,sol_r),
	          *const u_r   = get_col_const_Multiarray_T(1,sol_r),
	          *const v_r   = (DIM > 1 ? get_col_const_Multiarray_T(2,sol_r) : NULL),
	          *const w_r   = (DIM > 2 ? get_col_const_Multiarray_T(3,sol_r) : NULL),
	          *const p_r   = get_col_const_Multiarray_T(NVAR-1,sol_r);

	Type*const rho = get_col_Multiarray_T(0,sol),
	    *const u   = get_col_Multiarray_T(1,sol),
	    *const v   = (DIM > 1 ? get_col_Multiarray_T(2,sol) : NULL),
	    *const w   = (DIM > 2 ? get_col_Multiarray_T(3,sol) : NULL),
	    *const p   = get_col_Multiarray_T(NVAR-1,sol);

	for (int n = 0; n < n_n; n++) {
		const Type uvw_l[] = { u_l[n], (DIM > 1 ? v_l[n] : 0.0), (DIM > 2 ? w_l[n] : 0.0), },
		           uvw_r[] = { u_r[n], (DIM > 1 ? v_r[n] : 0.0), (DIM > 2 ? w_r[n] : 0.0), };
		const Type c_l = sqrt_T(GAMMA*p_l[n]/rho_l[n]),
		           c_r = sqrt_T(GAMMA*p_r[n]/rho_r[n]);

		const double* data_n = get_row_const_Multiarray_d(n,normals);
		const Type Vn_l = compute_Vn(data_n,uvw_l),
		           Vn_r = compute_Vn(data_n,uvw_r);

		// Riemann invariants
		const Type R_l = Vn_l + 2.0/GM1*c_l, // Outgoing
		           R_r = Vn_r - 2.0/GM1*c_r; // Incoming

		const Type Vn = 0.5*(R_l+R_r),
		           c  = 0.25*GM1*(R_l-R_r);

		Type* uvw[] = { &u[n], (DIM > 1 ? &v[n] : NULL), (DIM > 2 ? &w[n] : NULL), };
		if (abs_T(Vn) >= abs_T(c)) { // Supersonic
			if (real_T(Vn) < 0.0) { // Inlet
				rho[n] = rho_r[n];
				set_uvw(uvw_r,uvw);
				p[n]   = p_r[n];
			} else { // Outlet
				rho[n] = rho_l[n];
				set_uvw(uvw_l,uvw);
				p[n]   = p_l[n];
			}
		} else { // Subsonic
			Type uvw_t[] = { 0.0, 0.0, 0.0, };
			if (real_T(Vn) < 0.0) { // Inlet
				const Type s_r = p_r[n]/pow_T(rho_r[n],GAMMA);
				rho[n] = pow_T(1.0/GAMMA*c*c/(s_r),1.0/GM1);
				compute_Vt(data_n,Vn_r,uvw_r,uvw_t);
			} else { // Outlet
				const Type s_l = p_l[n]/pow_T(rho_l[n],GAMMA);
				rho[n] = pow_T(1.0/GAMMA*c*c/(s_l),1.0/GM1);
				compute_Vt(data_n,Vn_l,uvw_l,uvw_t);
			}
			compute_uvw(data_n,Vn,uvw_t,uvw);

			p[n] = 1.0/GAMMA*c*c*rho[n];
		}
	}

	if (c_m[1] == true) {

	struct Multiarray_T* ds_ds = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_n,NVAR,NVAR}); // keep

	// Standard datatypes
	unsigned int i, iMax, n, eq, var, ind_ds_ds;
	Type       rhoL, rhoL_inv, uL, vL, wL, V2L, pL, rhoR, uR, vR, wR, pR,
	           cL, RL, VnL, cR, RR, VnR, c, Vn;
	double     n1, n2, n3;

	// silence
	n2 = n3 = 0.0;

	const double* n_ptr = normals->data;

	Type *ds_ds_ptr[NVAR*NVAR];
	for (int ind = 0, vr_l = 0; vr_l < NVAR; vr_l++) {
	for (int vr_r = 0; vr_r < NVAR; vr_r++) {
		ds_ds_ptr[ind] = &ds_ds->data[n_n*(ind)];
		++ind;
	}}

	for (n = 0; n < n_n; n++) {
		ind_ds_ds = 0;

		// Inner VOLUME
		rhoL     = rho_l[n];
		rhoL_inv = 1.0/rhoL;

		uL = u_l[n];
		vL = ( DIM > 1 ? v_l[n] : 0.0 );
		wL = ( DIM > 2 ? w_l[n] : 0.0 );
		pL = p_l[n];

		V2L = uL*uL+vL*vL+wL*wL;

		n1 = *n_ptr++;
		if      (DIM == 3) { n2 = *n_ptr++; n3 = *n_ptr++; }
		else if (DIM == 2) { n2 = *n_ptr++; n3 = 0.0;      }
		else if (DIM == 1) { n2 = 0.0;      n3 = 0.0;      }

		VnL = uL*n1+vL*n2+wL*n3;

		// Outer VOLUME
		rhoR = rho_r[n];
		uR   = u_r[n];
		vR   = ( DIM > 1 ? v_r[n] : 0.0 );
		wR   = ( DIM > 2 ? w_r[n] : 0.0 );
		pR   = p_r[n];

		VnR = n1*uR+n2*vR+n3*wR;

		cL = sqrt_T(GAMMA*pL/rhoL);
		cR = sqrt_T(GAMMA*pR/rhoR);

		// Riemann invariants
		RL = VnL + 2.0/GM1*cL;
		RR = VnR - 2.0/GM1*cR;

		Vn = 0.5*(RL+RR);
		c  = 0.25*GM1*(RL-RR);

		if (abs_T(Vn) >= abs_T(c)) { // Supersonic
			if (real_T(Vn) < 0.0) { // Inlet
//printf("j: Sup Inlet\n");
				for (var = 0; var < NVAR; var++) {
				for (eq = 0; eq < NEQ; eq++) {
					*ds_ds_ptr[ind_ds_ds++] = 0.0;
				}}
			} else { // Outlet
//printf("j: Sup Outlet\n");
				for (var = 0; var < NVAR; var++) {
				for (eq = 0; eq < NEQ; eq++) {
					if (var != eq)
						*ds_ds_ptr[ind_ds_ds++] = 0.0;
					else
						*ds_ds_ptr[ind_ds_ds++] = 1.0;
				}}
			}
		} else { // Subsonic
			Type dcLdW, rho, u, v, w, V2, ut, vt, wt, un, vn, wn, cnst1, drhodW, dudW, dvdW, dwdW, dpdW,
			     drhoLdW[NVAR], duLdW[NVAR], dvLdW[NVAR], dwLdW[NVAR], dpLdW[NVAR],
			     dVnLdW[NVAR], dRLdW[NVAR], dcdW[NVAR];

			// silence
			un = vn = wn = 0.0;

			if (DIM == 3) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0; drhoLdW[4] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = -wL; dpLdW[4]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;      duLdW[4] = 0.0;
				dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;      dvLdW[4] = 0.0;
				dwLdW[0] = -wL*rhoL_inv; dwLdW[1] = 0.0;      dwLdW[2] = 0.0;      dwLdW[3] = rhoL_inv; dwLdW[4] = 0.0;

				for (var = 0; var < NVAR; var++) {
					dpLdW[var] *= GM1;
					dVnLdW[var] = duLdW[var]*n1 + dvLdW[var]*n2 + dwLdW[var]*n3;
				}

				un = Vn*n1;
				vn = Vn*n2;
				wn = Vn*n3;
			} else if (DIM == 2) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;
				dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;

				for (var = 0; var < NVAR; var++) {
					dpLdW[var] *= GM1;
					dVnLdW[var] = duLdW[var]*n1 + dvLdW[var]*n2;
				}

				un = Vn*n1;
				vn = Vn*n2;
			} else if (DIM == 1) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;

				for (var = 0; var < NVAR; var++) {
					dpLdW[var] *= GM1;
					dVnLdW[var] = duLdW[var]*n1;
				}

				un = Vn*n1;
			}

			for (var = 0; var < NVAR; var++) {
				dcLdW      = 0.5*GAMMA/(cL*rhoL*rhoL)*(dpLdW[var]*rhoL-pL*drhoLdW[var]);
				dRLdW[var] = dVnLdW[var] + 2.0/GM1*dcLdW;
				dcdW[var]  = 0.25*GM1*dRLdW[var];
			}

			if (real_T(Vn) < 0.0) { // Inlet
//printf("j: Sub Inlet\n");
				Type sR;

				sR  = sqrt_T(pR/pow_T(rhoR,GAMMA));
				if (DIM == 3) {
					for (var = 0; var < NVAR; var++) {
						drhodW = pow_T(GAMMA,-1.0/GM1)*2.0/(GM1*sR)*pow_T(c/sR,-GM3/GM1)*dcdW[var];

						rho = pow_T(c*c/(GAMMA*sR*sR),1.0/GM1);

						ut = uR - VnR*n1;
						vt = vR - VnR*n2;
						wt = wR - VnR*n3;

						u = un + ut;
						v = vn + vt;
						w = wn + wt;
						V2 = u*u+v*v+w*w;

						dudW = 0.5*dRLdW[var]*n1;
						dvdW = 0.5*dRLdW[var]*n2;
						dwdW = 0.5*dRLdW[var]*n3;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*ds_ds_ptr[ind_ds_ds++] = drhodW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*u + rho*dudW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*v + rho*dvdW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*w + rho*dwdW;
						*ds_ds_ptr[ind_ds_ds++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW+w*dwdW));
					}
				} else if (DIM == 2) {
					for (var = 0; var < NVAR; var++) {
						drhodW = pow_T(GAMMA,-1.0/GM1)*2.0/(GM1*sR)*pow_T(c/sR,-GM3/GM1)*dcdW[var];

						rho = pow_T(c*c/(GAMMA*sR*sR),1.0/GM1);

						ut = uR - VnR*n1;
						vt = vR - VnR*n2;

						u = un + ut;
						v = vn + vt;
						V2 = u*u+v*v;

						dudW = 0.5*dRLdW[var]*n1;
						dvdW = 0.5*dRLdW[var]*n2;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*ds_ds_ptr[ind_ds_ds++] = drhodW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*u + rho*dudW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*v + rho*dvdW;
						*ds_ds_ptr[ind_ds_ds++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW));
					}
				} else if (DIM == 1) {
					for (var = 0; var < NVAR; var++) {
						drhodW = pow_T(GAMMA,-1.0/GM1)*2.0/(GM1*sR)*pow_T(c/sR,-GM3/GM1)*dcdW[var];

						rho = pow_T(c*c/(GAMMA*sR*sR),1.0/GM1);

						ut = uR - VnR*n1;
						u  = un + ut;
						V2 = u*u;

						dudW = 0.5*dRLdW[var]*n1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*ds_ds_ptr[ind_ds_ds++] = drhodW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*u + rho*dudW;
						*ds_ds_ptr[ind_ds_ds++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW));
					}
				}
			} else { // Outlet
//printf("j: Sub Outlet\n");
				Type sL, dsLdW;

				sL  = sqrt_T(pL/pow_T(rhoL,GAMMA));
				if (DIM == 3) {
					for (var = 0; var < NVAR; var++) {
						dsLdW = 0.5*sqrt_T(pow_T(rhoL,GAMMA)/pL)/pow_T(rhoL,2.0*GAMMA)*
								(dpLdW[var]*pow_T(rhoL,GAMMA)-GAMMA*pL*pow_T(rhoL,GM1)*drhoLdW[var]);

						drhodW = pow_T(GAMMA,-1.0/GM1)*2.0/GM1*pow_T(c,-GM3/GM1)*
								 pow_T(sL,-(GAMMA+1.0)/GM1)*(dcdW[var]*sL-c*dsLdW);

						rho = pow_T(c*c/(GAMMA*sL*sL),1.0/GM1);

						ut = uL - VnL*n1;
						vt = vL - VnL*n2;
						wt = wL - VnL*n3;

						u = un + ut;
						v = vn + vt;
						w = wn + wt;
						V2 = u*u+v*v+w*w;

						cnst1 = 0.5*dRLdW[var]-dVnLdW[var];
						dudW = duLdW[var]+n1*cnst1;
						dvdW = dvLdW[var]+n2*cnst1;
						dwdW = dwLdW[var]+n3*cnst1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*ds_ds_ptr[ind_ds_ds++] = drhodW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*u + rho*dudW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*v + rho*dvdW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*w + rho*dwdW;
						*ds_ds_ptr[ind_ds_ds++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW+w*dwdW));
					}
				} else if (DIM == 2) {
					for (var = 0; var < NVAR; var++) {
						dsLdW = 0.5*sqrt_T(pow_T(rhoL,GAMMA)/pL)/pow_T(rhoL,2.0*GAMMA)*
								(dpLdW[var]*pow_T(rhoL,GAMMA)-GAMMA*pL*pow_T(rhoL,GM1)*drhoLdW[var]);

						drhodW = pow_T(GAMMA,-1.0/GM1)*2.0/GM1*pow_T(c,-GM3/GM1)*
								 pow_T(sL,-(GAMMA+1.0)/GM1)*(dcdW[var]*sL-c*dsLdW);

						rho = pow_T(c*c/(GAMMA*sL*sL),1.0/GM1);

						ut = uL - VnL*n1;
						vt = vL - VnL*n2;

						u = un + ut;
						v = vn + vt;
						V2 = u*u+v*v;

						cnst1 = 0.5*dRLdW[var]-dVnLdW[var];
						dudW = duLdW[var]+n1*cnst1;
						dvdW = dvLdW[var]+n2*cnst1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*ds_ds_ptr[ind_ds_ds++] = drhodW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*u + rho*dudW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*v + rho*dvdW;
						*ds_ds_ptr[ind_ds_ds++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW));
					}
				} else if (DIM == 1) {
					for (var = 0; var < NVAR; var++) {
						dsLdW = 0.5*sqrt_T(pow_T(rhoL,GAMMA)/pL)/pow_T(rhoL,2.0*GAMMA)*
								(dpLdW[var]*pow_T(rhoL,GAMMA)-GAMMA*pL*pow_T(rhoL,GM1)*drhoLdW[var]);

						drhodW = pow_T(GAMMA,-1.0/GM1)*2.0/GM1*pow_T(c,-GM3/GM1)*
								 pow_T(sL,-(GAMMA+1.0)/GM1)*(dcdW[var]*sL-c*dsLdW);

						rho = pow_T(c*c/(GAMMA*sL*sL),1.0/GM1);

						ut = uL - VnL*n1;
						u = un + ut;
						V2 = u*u;

						cnst1 = 0.5*dRLdW[var]-dVnLdW[var];
						dudW = duLdW[var]+n1*cnst1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*ds_ds_ptr[ind_ds_ds++] = drhodW;
						*ds_ds_ptr[ind_ds_ds++] = drhodW*u + rho*dudW;
						*ds_ds_ptr[ind_ds_ds++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW));
					}
				}
			}
		}
		for (i = 0, iMax = NEQ*NVAR; i < iMax; i++)
			ds_ds_ptr[i]++;
	}
	bv->ds_ds = (const struct const_Multiarray_T*) ds_ds;
	}
	convert_variables_T(sol,'p','c');
	bv->s = (struct const_Multiarray_T*)sol;

	convert_variables_T((struct Multiarray_T*)sol_l,'p','c');
	destructor_const_Multiarray_T(sol_r);

	constructor_Boundary_Value_T_grad_from_internal(bv,bv_i,NVAR);
}

void constructor_Boundary_Value_T_euler_slipwall
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
	 const struct Simulation* sim)
{
	UNUSED(face);
	UNUSED(sim);
	const bool* c_m = bv_i->compute_member;
	assert(c_m[0] == true);

	const struct const_Multiarray_T* sol_l = bv_i->s;

	const Type*const rho_l  = get_col_const_Multiarray_T(0,sol_l),
	          *const rhou_l = get_col_const_Multiarray_T(1,sol_l),
	          *const rhov_l = (DIM > 1 ? get_col_const_Multiarray_T(2,sol_l) : NULL),
	          *const rhow_l = (DIM > 2 ? get_col_const_Multiarray_T(3,sol_l) : NULL),
	          *const E_l    = get_col_const_Multiarray_T(NVAR-1,sol_l);

	const ptrdiff_t n_n = sol_l->extents[0];
	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,NVAR}); // keep

	Type*const rho  = get_col_Multiarray_T(0,sol),
	    *const rhou = get_col_Multiarray_T(1,sol),
	    *const rhov = (DIM > 1 ? get_col_Multiarray_T(2,sol) : NULL),
	    *const rhow = (DIM > 2 ? get_col_Multiarray_T(3,sol) : NULL),
	    *const E    = get_col_Multiarray_T(NVAR-1,sol);

	const struct const_Multiarray_d* normals = bv_i->normals;
	assert(normals->layout == 'R');

	for (int n = 0; n < n_n; n++) {
		rho[n] = rho_l[n];
		E[n]   = E_l[n];   // Equivalent to setting p[n] = p_l[n] (See comments).

		const double* data_n = get_row_const_Multiarray_d(n,normals);
		const Type rhouvw_l[] = { rhou_l[n], (DIM > 1 ? rhov_l[n] : 0.0), (DIM > 2 ? rhow_l[n] : 0.0), };
		const Type rhoVn_l = compute_Vn(data_n,rhouvw_l);

		Type* rhouvw[] = { &rhou[n], (DIM > 1 ? &rhov[n] : NULL), (DIM > 2 ? &rhow[n] : NULL), };
		compute_opposite_normal_uvw(data_n,rhoVn_l,rhouvw_l,rhouvw);
	}
	bv->s = (struct const_Multiarray_T*)sol;

	if (c_m[1] == true) {
		struct Multiarray_T* ds_ds = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_n,NVAR,NVAR}); // keep

		Type *ds_ds_ptr[NVAR*NVAR];
		for (int ind = 0, vr_l = 0; vr_l < NVAR; vr_l++) {
		for (int vr_r = 0; vr_r < NVAR; vr_r++) {
			ds_ds_ptr[ind] = &ds_ds->data[n_n*(ind)];
			++ind;
		}}
		const double* n_ptr = get_row_const_Multiarray_d(0,normals);

		if (DIM == 3) {
			for (int n = 0; n < n_n; n++) {
				const Type n1 = *n_ptr++,
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
				const Type n1 = *n_ptr++,
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
				const Type n1 = *n_ptr++;

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
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds;
	}
	assert(c_m[2] == false);
}

void constructor_Boundary_Value_T_euler_supersonic_inflow
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
	 const struct Simulation* sim)
{
	UNUSED(face);
	const bool* c_m = bv_i->compute_member;
	assert(c_m[0] == true);

	const struct const_Multiarray_d* xyz = bv_i->xyz;
	const struct const_Multiarray_T* sol_r = constructor_sol_bv(xyz,sim); // destructed

	const ptrdiff_t n_n = xyz->extents[0];

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,NVAR}); // keep

	const Type*const rho_r  = get_col_const_Multiarray_T(0,sol_r),
	          *const rhou_r = get_col_const_Multiarray_T(1,sol_r),
	          *const E_r    = get_col_const_Multiarray_T(NVAR-1,sol_r);

	IF_DIM_GE_2( const Type*const rhov_r = (DIM > 1 ? get_col_const_Multiarray_T(2,sol_r) : NULL); )
	IF_DIM_GE_3( const Type*const rhow_r = (DIM > 2 ? get_col_const_Multiarray_T(3,sol_r) : NULL); )

	Type*const rho  = get_col_Multiarray_T(0,sol),
	    *const rhou = get_col_Multiarray_T(1,sol),
	    *const E    = get_col_Multiarray_T(NVAR-1,sol);

	IF_DIM_GE_2( Type*const rhov = (DIM > 1 ? get_col_Multiarray_T(2,sol) : NULL); )
	IF_DIM_GE_3( Type*const rhow = (DIM > 2 ? get_col_Multiarray_T(3,sol) : NULL); )

	for (int n = 0; n < n_n; n++) {
		rho[n] = rho_r[n];
		IF_DIM_GE_1( rhou[n] = rhou_r[n]; )
		IF_DIM_GE_2( rhov[n] = rhov_r[n]; )
		IF_DIM_GE_3( rhow[n] = rhow_r[n]; )
		E[n]   = E_r[n];
	}

	bv->s = (struct const_Multiarray_T*)sol;

	if (c_m[1] == true) {
		struct Multiarray_T* ds_ds = constructor_zero_Multiarray_T('C',3,(ptrdiff_t[]){n_n,NVAR,NVAR}); // keep
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds;
	}
	destructor_const_Multiarray_T(sol_r);

	constructor_Boundary_Value_T_grad_from_internal(bv,bv_i,NVAR);
}

void constructor_Boundary_Value_T_euler_supersonic_outflow
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
	 const struct Simulation* sim)
{
	UNUSED(face);
	UNUSED(sim);
	const bool* c_m = bv_i->compute_member;
	assert(c_m[0] == true);

	const struct const_Multiarray_d* xyz = bv_i->xyz;
	const struct const_Multiarray_T* sol_l = bv_i->s;

	const ptrdiff_t n_n = xyz->extents[0];

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,NVAR}); // keep

	const Type*const rho_l  = get_col_const_Multiarray_T(0,sol_l),
	          *const rhou_l = get_col_const_Multiarray_T(1,sol_l),
	          *const E_l    = get_col_const_Multiarray_T(NVAR-1,sol_l);

	IF_DIM_GE_2( const Type*const rhov_l = (DIM > 1 ? get_col_const_Multiarray_T(2,sol_l) : NULL); )
	IF_DIM_GE_3( const Type*const rhow_l = (DIM > 2 ? get_col_const_Multiarray_T(3,sol_l) : NULL); )

	Type*const rho  = get_col_Multiarray_T(0,sol),
	    *const rhou = get_col_Multiarray_T(1,sol),
	    *const E    = get_col_Multiarray_T(NVAR-1,sol);

	IF_DIM_GE_2( Type*const rhov = (DIM > 1 ? get_col_Multiarray_T(2,sol) : NULL); )
	IF_DIM_GE_3( Type*const rhow = (DIM > 2 ? get_col_Multiarray_T(3,sol) : NULL); )

	for (int n = 0; n < n_n; n++) {
		rho[n] = rho_l[n];
		IF_DIM_GE_1( rhou[n] = rhou_l[n]; )
		IF_DIM_GE_2( rhov[n] = rhov_l[n]; )
		IF_DIM_GE_3( rhow[n] = rhow_l[n]; )
		E[n]   = E_l[n];
	}
	bv->s = (struct const_Multiarray_T*)sol;

	if (c_m[1] == true) {
		struct Multiarray_T* ds_ds = constructor_zero_Multiarray_T('C',3,(ptrdiff_t[]){n_n,NVAR,NVAR}); // keep

		Type *ds_ds_ptr[NVAR*NVAR];
		for (int ind = 0, vr_l = 0; vr_l < NVAR; vr_l++) {
		for (int vr_r = 0; vr_r < NVAR; vr_r++) {
			ds_ds_ptr[ind] = get_col_Multiarray_T(ind,ds_ds);
			++ind;
		}}

		for (int n = 0; n < n_n; n++) {
			int ind_ds_ds = 0;
			for (int var = 0; var < NVAR; var++) {
			for (int eq = 0; eq < NEQ; eq++) {
				if (var != eq)
					*ds_ds_ptr[ind_ds_ds++] = 0.0;
				else
					*ds_ds_ptr[ind_ds_ds++] = 1.0;
			}}
			for (int i = 0; i < NVAR*NVAR; i++)
				ds_ds_ptr[i]++;
		}
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds;
	}

	constructor_Boundary_Value_T_grad_from_internal(bv,bv_i,NVAR);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static Type compute_Vn (const double*const n, const Type*const uvw)
{
	Type Vn = 0.0;
	for (int d = 0; d < DIM; ++d)
		Vn += n[d]*uvw[d];
	return Vn;
}

static void set_uvw (const Type*const uvw_i, Type*const*const uvw)
{
	for (int d = 0; d < DIM; ++d)
		*uvw[d] = uvw_i[d];
}

static void compute_Vt
	(const double*const n, const Type Vn, const Type*const uvw, Type*const uvw_t)
{
	for (int d = 0; d < DIM; ++d)
		uvw_t[d] = uvw[d] - Vn*n[d];
}

static void compute_uvw
	(const double*const n, const Type Vn, const Type*const uvw_t, Type*const*const uvw)
{
	for (int d = 0; d < DIM; ++d)
		*uvw[d] = Vn*n[d] + uvw_t[d];
}

static void compute_opposite_normal_uvw
	(const double*const n, const Type Vn, const Type*const uvw_i, Type*const*const uvw)
{
	for (int d = 0; d < DIM; ++d)
		*uvw[d] = uvw_i[d] - 2.0*Vn*n[d];
}
