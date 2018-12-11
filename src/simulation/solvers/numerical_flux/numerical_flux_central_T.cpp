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
 *  \brief Provides the templated diffusion numerical flux functions.
 */

#include <assert.h>
#include <stddef.h>

#include "macros.h"
#include "definitions_core.h"


#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_boundary.h"
#include "def_templates_flux.h"

// Static function declarations ************************************************************************************* //

/** \brief Container for \ref Flux_T containers from the 'L'eft and 'R'ight sides which used to compute the central
 *         numerical flux members. */
struct Fluxes_LR {
	struct Flux_T* flux_l, ///< Flux members from the 'l'eft side.
	             * flux_r; ///< Flux members from the 'r'ight side.
};

/** \brief Constructor for a statically allocated \ref Fluxes_LR container.
 *  \return See brief. */
static struct Fluxes_LR constructor_Fluxes_LR
	(const struct Numerical_Flux_Input_T*const num_flux_i ///< \ref Numerical_Flux_Input_T.
	);

/// \brief Destructor for a statically allocated \ref Fluxes_LR container.
static void destructor_Fluxes_LR
	(struct Fluxes_LR*const fluxes_lr ///< Standard.
	);

/// \brief Sets `nnf_a` = scale* col_sum(dot_mult(`n_l`,`f_a`)) for each slice of `f_a`.
static void set_Numerical_Flux_member
	(const struct const_Multiarray_T*const n_l, ///< The 'n'ormal vectors as seen from the 'l'eft volume.
	 struct Multiarray_T*const f_a,             ///< Pointer to the 'a'rbitrary 'f'lux member.
	 struct Multiarray_T*const nnf_a,           ///< Pointer to the 'a'rbitrary 'n'ormal 'n'umerical 'f'lux member.
	 const Real scale,                          ///< Scaling constant.
	 const int loop_max[3]                      ///< The maximum values for the loop counters.
	);

/// \brief Set any normal numerical flux members which are provided directly from the boundary conditions.
static void set_provided_Numerical_Flux_members
	(const struct Numerical_Flux_Input_T*const num_flux_i, ///< \ref Numerical_Flux_Input_T.
	 struct mutable_Numerical_Flux_T*const num_flux        ///< \ref mutable_Numerical_Flux_T.
	);

/// \brief Set any normal numerical flux Jacobian members which are provided directly from the boundary conditions.
static void set_provided_Numerical_Flux_jacobian_members
	(const struct Numerical_Flux_Input_T*const num_flux_i, ///< \ref Numerical_Flux_Input_T.
	 struct mutable_Numerical_Flux_T*const num_flux        ///< \ref mutable_Numerical_Flux_T.
	);

// Interface functions ********************************************************************************************** //

/** \brief Version of \ref compute_Numerical_Flux_fptr_T computing the numerical fluxes as the average of the left and
 *         right fluxes. */
static void compute_Numerical_Flux_T_central
	(const struct Numerical_Flux_Input_T* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_T* num_flux        ///< See brief.
	)
{
	struct Fluxes_LR fluxes_lr = constructor_Fluxes_LR(num_flux_i); // destructed
	struct Flux_T*const flux_l = fluxes_lr.flux_l,
	             *const flux_r = fluxes_lr.flux_r;

	const struct const_Multiarray_T*const n_l = num_flux_i->bv_l.normals;
	transpose_Multiarray_T((struct Multiarray_T*)n_l,true);

	struct Multiarray_T*const f_avg = (struct Multiarray_T*)
		constructor_sum_Multiarrays_const_Multiarray_T(0.5,flux_l->f,0.5,flux_r->f); // destructed

	const int n_eq_set = ( !num_flux_i->bv_r.nf_E_provided ? NEQ : NEQ-1 );
	set_Numerical_Flux_member(n_l,f_avg,num_flux->nnf,1.0,(int[]){n_eq_set,1,1});
	if (num_flux_i->bv_r.nf_E_provided)
		set_provided_Numerical_Flux_members(num_flux_i,num_flux);

	transpose_Multiarray_T((struct Multiarray_T*)n_l,true);

	destructor_Fluxes_LR(&fluxes_lr);
/// \todo Ensure that all is working correctly for nonlinear BCs.
// Different treatment for boundary faces in previous implementation (evaluate flux at common states). Result is
// identical for linear equations.
}

/** \brief Version of \ref compute_Numerical_Flux_fptr_T computing the numerical fluxes and Jacobians for the flux of
 *         \ref compute_Numerical_Flux_T_central. */
static void compute_Numerical_Flux_T_central_jacobian
	(const struct Numerical_Flux_Input_T* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_T* num_flux        ///< See brief.
	)
{
	struct Fluxes_LR fluxes_lr = constructor_Fluxes_LR(num_flux_i); // destructed
	struct Flux_T*const flux_l = fluxes_lr.flux_l,
	             *const flux_r = fluxes_lr.flux_r;

	struct Flux_Input_T*const flux_i = num_flux_i->flux_i;

	const struct const_Multiarray_T*const n_l = num_flux_i->bv_l.normals;
	transpose_Multiarray_T((struct Multiarray_T*)n_l,true);

	const bool*const c_m = flux_i->compute_member;

	const int n_eq_set = ( !num_flux_i->bv_r.nf_E_provided ? NEQ : NEQ-1 );
	assert(c_m[0]);
	{
	struct Multiarray_T*const f_avg = (struct Multiarray_T*)
		constructor_sum_Multiarrays_const_Multiarray_T(0.5,flux_l->f,0.5,flux_r->f); // destructed
	set_Numerical_Flux_member(n_l,f_avg,num_flux->nnf,1.0,(int[]){n_eq_set,1,1});
	destructor_Multiarray_T(f_avg);
	}

	if (c_m[1])
	{
	struct Multiarray_T*const df_ds_l = constructor_copy_Multiarray_T((struct Multiarray_T*)flux_l->df_ds); // dest.
	set_Numerical_Flux_member(n_l,df_ds_l,num_flux->neigh_info[0].dnnf_ds,0.5,(int[]){n_eq_set,NVR,1});
	destructor_Multiarray_T(df_ds_l);

	struct Multiarray_T*const df_ds_r = constructor_copy_Multiarray_T((struct Multiarray_T*)flux_r->df_ds); // dest.
	set_Numerical_Flux_member(n_l,df_ds_r,num_flux->neigh_info[1].dnnf_ds,0.5,(int[]){n_eq_set,NVR,1});
	destructor_Multiarray_T(df_ds_r);
	}

	assert(c_m[2]);
	{
	struct Multiarray_T*const df_dg_l = constructor_copy_Multiarray_T((struct Multiarray_T*)flux_l->df_dg); // dest.
	set_Numerical_Flux_member(n_l,df_dg_l,num_flux->neigh_info[0].dnnf_dg,0.5,(int[]){n_eq_set,NVR,DIM});
	destructor_Multiarray_T(df_dg_l);

	struct Multiarray_T*const df_dg_r = constructor_copy_Multiarray_T((struct Multiarray_T*)flux_r->df_dg); // dest.
	set_Numerical_Flux_member(n_l,df_dg_r,num_flux->neigh_info[1].dnnf_dg,0.5,(int[]){n_eq_set,NVR,DIM});
	destructor_Multiarray_T(df_dg_r);
	}

	if (num_flux_i->bv_r.nf_E_provided)
		set_provided_Numerical_Flux_jacobian_members(num_flux_i,num_flux);

	assert(c_m[3] == 0); // Add support.
	assert(c_m[4] == 0); // Add support.
	assert(c_m[5] == 0); // Add support.

	transpose_Multiarray_T((struct Multiarray_T*)n_l,true);

	destructor_Fluxes_LR(&fluxes_lr);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Sets each entry of `nnf_a` = val for each column in the multiarray corresponding to the Energy equation.
static void set_Numerical_Flux_Energy_member
	(struct Multiarray_T*const nnf_a, ///< Pointer to the 'a'rbitrary 'n'ormal 'n'umerical 'f'lux member.
	 const Type val,                  ///< The desired value.
	 const int loop_max[2]            ///< The maximum values for the loop counters.
	);

static struct Fluxes_LR constructor_Fluxes_LR (const struct Numerical_Flux_Input_T*const num_flux_i)
{
	static struct Fluxes_LR fluxes_lr;

	// Note that only the 2nd order flux contributions should be included in the viscous numerical flux.
	struct Flux_Input_T*const flux_i = num_flux_i->flux_i;
	compute_Flux_fptr_T compute_Flux_initial = flux_i->compute_Flux;
	flux_i->compute_Flux = compute_Flux_2_T;

	flux_i->s = num_flux_i->bv_l.s;
	flux_i->g = num_flux_i->bv_l.g;
	fluxes_lr.flux_l = constructor_Flux_T(flux_i); // destructed

	flux_i->s = num_flux_i->bv_r.s;
	flux_i->g = num_flux_i->bv_r.g;
	fluxes_lr.flux_r = constructor_Flux_T(flux_i); // destructed

	flux_i->compute_Flux = compute_Flux_initial;

	return fluxes_lr;
}

static void destructor_Fluxes_LR (struct Fluxes_LR*const fluxes_lr)
{
	destructor_Flux_T(fluxes_lr->flux_l);
	destructor_Flux_T(fluxes_lr->flux_r);
}

static void set_Numerical_Flux_member
	(const struct const_Multiarray_T*const n_l, struct Multiarray_T*const f_a, struct Multiarray_T*const nnf_a,
	 const Real scale, const int loop_max[3])
{
	for (int d  = 0; d  < loop_max[2]; ++d)  {
	for (int vr = 0; vr < loop_max[1]; ++vr) {
	for (int eq = 0; eq < loop_max[0]; ++eq) {
		struct Multiarray_T f_s = interpret_Multiarray_as_slice_T(f_a,2,(ptrdiff_t[]){eq,vr,d});
		multiply_in_place_Multiarray_T(scale,&f_s,n_l);

		const struct const_Matrix_T f_M = interpret_const_Multiarray_as_Matrix_T((struct const_Multiarray_T*)&f_s);
		struct Vector_T nnf_V = interpret_Multiarray_slice_as_Vector_T(nnf_a,(ptrdiff_t[]){eq,vr,d});

		set_to_sum_Vector_T('C',&f_M,&nnf_V,false);
	}}}
}

static void set_provided_Numerical_Flux_members
	(const struct Numerical_Flux_Input_T*const num_flux_i, struct mutable_Numerical_Flux_T*const num_flux)
{
	const bool*const c_m = num_flux_i->flux_i->compute_member;

	assert(c_m[0]);
	set_Numerical_Flux_Energy_member(num_flux->nnf,num_flux_i->bv_r.nf_E,(int[]){1,1});

	const_cast_b(&num_flux_i->bv_r.nf_E_provided,false);
}

static void set_provided_Numerical_Flux_jacobian_members
	(const struct Numerical_Flux_Input_T*const num_flux_i, struct mutable_Numerical_Flux_T*const num_flux)
{
	set_provided_Numerical_Flux_members(num_flux_i,num_flux);

	const bool*const c_m = num_flux_i->flux_i->compute_member;

	if (c_m[1]) {
		set_Numerical_Flux_Energy_member(num_flux->neigh_info[0].dnnf_ds,0.0,(int[]){NVR,1});
		set_Numerical_Flux_Energy_member(num_flux->neigh_info[1].dnnf_ds,0.0,(int[]){NVR,1});
	}

	assert(c_m[2]);
	set_Numerical_Flux_Energy_member(num_flux->neigh_info[0].dnnf_dg,0.0,(int[]){NVR,DIM});
	set_Numerical_Flux_Energy_member(num_flux->neigh_info[1].dnnf_dg,0.0,(int[]){NVR,DIM});

	const_cast_b(&num_flux_i->bv_r.nf_E_provided,false);
}

// Level 1 ********************************************************************************************************** //

static void set_Numerical_Flux_Energy_member
	(struct Multiarray_T*const nnf_a, const Type val, const int loop_max[2])
{
	const int eq = NEQ-1;
	for (int d  = 0; d  < loop_max[1]; ++d)  {
	for (int vr = 0; vr < loop_max[0]; ++vr) {
		struct Vector_T nnf_V = interpret_Multiarray_slice_as_Vector_T(nnf_a,(ptrdiff_t[]){eq,vr,d});
		set_to_value_Vector_T(&nnf_V,val);
	}}
}

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_boundary.h"
#include "undef_templates_flux.h"
