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


#include "def_templates_numerical_flux.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_boundary.h"
#include "def_templates_flux.h"

// Static function declarations ************************************************************************************* //

#define NEQ 1 ///< Number of equations.
#define NVR 1 ///< Number of variables.

/// \brief Sets `nnf_a` = scale* col_sum(dot_mult(`n_l`,`f_a`)) for each slice of `f_a`.
static void set_Numerical_Flux_member
	(const struct const_Multiarray_R*const n_l, ///< The 'n'ormal vectors as seen from the 'l'eft volume.
	 struct Multiarray_T*const f_a,             ///< Pointer to the 'a'rbitrary 'f'lux member.
	 struct Multiarray_T*const nnf_a,           ///< Pointer to the 'a'rbitrary 'n'ormal 'n'umerical 'f'lux member.
	 const Real scale,                          ///< Scaling constant.
	 const int loop_max[3]                      ///< The maximum values for the loop counters.
	);

// Interface functions ********************************************************************************************** //

void compute_Numerical_Flux_T_diffusion_central
	(const struct Numerical_Flux_Input_T* num_flux_i, struct mutable_Numerical_Flux_T* num_flux)
{
	struct Flux_Input_T*const flux_i = num_flux_i->flux_i;

	flux_i->g = num_flux_i->bv_l.g;
	struct Flux_T*const flux_l = constructor_Flux_T(flux_i); // destructed

	flux_i->g = num_flux_i->bv_r.g;
	struct Flux_T*const flux_r = constructor_Flux_T(flux_i); // destructed

	const struct const_Multiarray_R*const n_l = num_flux_i->bv_l.normals;
	transpose_Multiarray_R((struct Multiarray_R*)n_l,true);

	struct Multiarray_T*const f_avg = (struct Multiarray_T*)
		constructor_sum_Multiarrays_const_Multiarray_T(0.5,flux_l->f,0.5,flux_r->f); // destructed
	set_Numerical_Flux_member(n_l,f_avg,num_flux->nnf,1.0,(int[]){NEQ,1,1});

	transpose_Multiarray_R((struct Multiarray_R*)n_l,true);

	destructor_Flux_T(flux_l);
	destructor_Flux_T(flux_r);

/// \todo Ensure that all is working correctly for nonlinear BCs.
// Different treatment for boundary faces in previous implementation (evaluate flux at common states). Result is
// identical for linear equations.
}

void compute_Numerical_Flux_T_diffusion_central_jacobian
	(const struct Numerical_Flux_Input_T* num_flux_i, struct mutable_Numerical_Flux_T* num_flux)
{
	struct Flux_Input_T*const flux_i = num_flux_i->flux_i;

	flux_i->g = num_flux_i->bv_l.g;
	struct Flux_T*const flux_l = constructor_Flux_T(flux_i); // destructed

	flux_i->g = num_flux_i->bv_r.g;
	struct Flux_T*const flux_r = constructor_Flux_T(flux_i); // destructed

	const struct const_Multiarray_R*const n_l = num_flux_i->bv_l.normals;
	transpose_Multiarray_R((struct Multiarray_R*)n_l,true);


	struct Multiarray_T*const f_avg = (struct Multiarray_T*)
		constructor_sum_Multiarrays_const_Multiarray_T(0.5,flux_l->f,0.5,flux_r->f); // destructed
	set_Numerical_Flux_member(n_l,f_avg,num_flux->nnf,1.0,(int[]){NEQ,1,1});
	destructor_Multiarray_T(f_avg);

	struct Multiarray_T*const df_dg_l = constructor_copy_Multiarray_T((struct Multiarray_T*)flux_l->df_dg); // dest.
	set_Numerical_Flux_member(n_l,df_dg_l,num_flux->neigh_info[0].dnnf_dg,0.5,(int[]){NEQ,NVR,DIM});
	destructor_Multiarray_T(df_dg_l);

	struct Multiarray_T*const df_dg_r = constructor_copy_Multiarray_T((struct Multiarray_T*)flux_r->df_dg); // dest.
	set_Numerical_Flux_member(n_l,df_dg_r,num_flux->neigh_info[1].dnnf_dg,0.5,(int[]){NEQ,NVR,DIM});
	destructor_Multiarray_T(df_dg_r);


	transpose_Multiarray_R((struct Multiarray_R*)n_l,true);

	destructor_Flux_T(flux_l);
	destructor_Flux_T(flux_r);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_Numerical_Flux_member
	(const struct const_Multiarray_R*const n_l, struct Multiarray_T*const f_a, struct Multiarray_T*const nnf_a,
	 const Real scale, const int loop_max[3])
{
	for (int d  = 0; d  < loop_max[2]; ++d)  {
	for (int vr = 0; vr < loop_max[1]; ++vr) {
	for (int eq = 0; eq < loop_max[0]; ++eq) {
		struct Multiarray_T f_s = interpret_Multiarray_as_slice_T(f_a,2,(ptrdiff_t[]){eq,vr,d});
		multiply_in_place_Multiarray_TR(scale,&f_s,n_l);

		const struct const_Matrix_T f_M = interpret_const_Multiarray_as_Matrix_T((struct const_Multiarray_T*)&f_s);
		struct Vector_T nnf_V = interpret_Multiarray_slice_as_Vector_T(nnf_a,(ptrdiff_t[]){eq,vr,d});

		set_to_sum_Vector_T('C',&f_M,&nnf_V);
	}}}
}
