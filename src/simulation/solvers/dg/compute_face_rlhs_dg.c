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

#include "compute_face_rlhs_dg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "face_solver_dg.h"
#include "element_solver.h"
#include "volume.h"
#include "volume_solver_dg.h"

#include "computational_elements.h"
#include "compute_rlhs.h"
#include "compute_face_rlhs.h"
#include "multiarray_operator.h"
#include "numerical_flux.h"
#include "operator.h"
#include "simulation.h"
#include "solve_dg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Scale \ref Numerical_Flux_T::nnf and \ref Numerical_Flux_T::Neigh_Info_NF_T::dnnf_ds by the face Jacobian.
static void scale_by_Jacobian_i1
	(const struct Numerical_Flux* num_flux, ///< Defined for \ref scale_by_Jacobian_fptr_T.
	 struct Face* face,                     ///< Defined for \ref scale_by_Jacobian_fptr_T.
	 const struct Simulation* sim           ///< Defined for \ref scale_by_Jacobian_fptr_T.
	);

/// \brief Compute the rhs and first order lhs terms.
static void compute_rlhs_1
	(const struct Numerical_Flux* num_flux,     ///< Defined for \ref compute_rlhs_fptr.
	 struct DG_Solver_Face* dg_s_face,          ///< Defined for \ref compute_rlhs_fptr.
	 struct Solver_Storage_Implicit* s_store_i, ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim               ///< Defined for \ref compute_rlhs_fptr.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_face_rlhs_dg_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Finalize the lhs term contribution from the \ref Face for the dg scheme.
static void finalize_lhs_f_dg
	(const int side_index[2],                  /**< The indices of the affectee, affector, respectively. See the
	                                            *   comments in \ref compute_face_rlhs_dg.h for the convention. */
	 const struct Numerical_Flux* num_flux,    ///< Defined for \ref compute_rlhs_fptr.
	 struct DG_Solver_Face* dg_s_face,         ///< Defined for \ref compute_rlhs_fptr.
	 struct Solver_Storage_Implicit* s_store_i ///< Defined for \ref compute_rlhs_fptr.
	);

static void scale_by_Jacobian_i1
	(const struct Numerical_Flux* num_flux, struct Face* face, const struct Simulation* sim)
{
UNUSED(sim);
	assert((!face->boundary && num_flux->neigh_info[1].dnnf_ds != NULL) ||
	       ( face->boundary && num_flux->neigh_info[1].dnnf_ds == NULL));

	struct Solver_Face* s_face = (struct Solver_Face*)face;

	const struct const_Vector_d jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->nnf,&jacobian_det_fc,false);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->neigh_info[0].dnnf_ds,&jacobian_det_fc,false);
	if (!face->boundary)
		scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->neigh_info[1].dnnf_ds,&jacobian_det_fc,false);
}

static void compute_rlhs_1
	(const struct Numerical_Flux* num_flux, struct DG_Solver_Face* dg_s_face,
	 struct Solver_Storage_Implicit* s_store_i, const struct Simulation* sim)
{
	const struct Face* face = (struct Face*) dg_s_face;

	/// See \ref compute_face_rlhs_dg.h for the `lhs_**` notation.
	compute_rhs_f_dg(num_flux,dg_s_face,s_store_i,sim);

	finalize_lhs_f_dg((int[]){0,0},num_flux,dg_s_face,s_store_i); // lhs_ll
//EXIT_UNSUPPORTED;
	if (!face->boundary) {
//printf("\n\n\n");
		finalize_lhs_f_dg((int[]){0,1},num_flux,dg_s_face,s_store_i); // lhs_lr

		for (int i = 0; i < 2; ++i) {
			const struct Neigh_Info_NF* n_i = &num_flux->neigh_info[i];
			permute_Multiarray_d_fc((struct Multiarray_d*)n_i->dnnf_ds,'R',1,(struct Solver_Face*)face);
			scale_Multiarray_d((struct Multiarray_d*)n_i->dnnf_ds,-1.0); // Use "-ve" normal.
		}

		finalize_lhs_f_dg((int[]){1,0},num_flux,dg_s_face,s_store_i); // lhs_rl
		finalize_lhs_f_dg((int[]){1,1},num_flux,dg_s_face,s_store_i); // lhs_rr
//EXIT_UNSUPPORTED;
	}
}

// Level 1 ********************************************************************************************************** //

static void finalize_lhs_f_dg
	(const int side_index[2], const struct Numerical_Flux* num_flux, struct DG_Solver_Face* dg_s_face,
	 struct Solver_Storage_Implicit* s_store_i)
{
	struct Face* face              = (struct Face*) dg_s_face;
	struct Solver_Face* s_face     = (struct Solver_Face*) face;

	struct Matrix_d* lhs = constructor_lhs_f_1(side_index,num_flux,s_face); // destructed

	struct Solver_Volume* s_vol[2] = { (struct Solver_Volume*) face->neigh_info[0].volume,
	                                   (struct Solver_Volume*) face->neigh_info[1].volume, };
//printf("%td %td\n",s_vol[side_index[0]]->ind_dof,s_vol[side_index[1]]->ind_dof);
	set_petsc_Mat_row_col(s_store_i,s_vol[side_index[0]],0,s_vol[side_index[1]],0);
	add_to_petsc_Mat(s_store_i,(struct const_Matrix_d*)lhs);

	destructor_Matrix_d(lhs);
}
