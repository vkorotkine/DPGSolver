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

#include "definitions_test_case.h"

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
	(const struct Numerical_Flux*const num_flux, ///< Defined for \ref scale_by_Jacobian_fptr_T.
	 struct Face*const face,                     ///< Defined for \ref scale_by_Jacobian_fptr_T.
	 const struct Simulation*const sim           ///< Defined for \ref scale_by_Jacobian_fptr_T.
	);

/// \brief Scale \ref Numerical_Flux_T::nnf and \ref Numerical_Flux_T::Neigh_Info_NF_T::dnnf_dg by the face Jacobian.
static void scale_by_Jacobian_i2
	(const struct Numerical_Flux*const num_flux, ///< Defined for \ref scale_by_Jacobian_fptr_T.
	 struct Face*const face,                     ///< Defined for \ref scale_by_Jacobian_fptr_T.
	 const struct Simulation*const sim           ///< Defined for \ref scale_by_Jacobian_fptr_T.
	);

/// \brief Compute the rhs and lhs terms for 1st order fluxes.
static void compute_rlhs_1
	(const struct Numerical_Flux* num_flux, ///< Defined for \ref compute_rlhs_fptr.
	 struct DG_Solver_Face* dg_s_face,      ///< Defined for \ref compute_rlhs_fptr.
	 struct Solver_Storage_Implicit* ssi,   ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim           ///< Defined for \ref compute_rlhs_fptr.
	);

/// \brief Compute the rhs and lhs terms for fluxes depending only on solution gradients.
static void compute_rlhs_2
	(const struct Numerical_Flux* num_flux, ///< Defined for \ref compute_rlhs_fptr.
	 struct DG_Solver_Face* dg_s_face,      ///< Defined for \ref compute_rlhs_fptr.
	 struct Solver_Storage_Implicit* ssi,   ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim           ///< Defined for \ref compute_rlhs_fptr.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_face_rlhs_dg_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Finalize the 1st order lhs term contribution from the \ref Face for the dg scheme.
static void finalize_lhs_1_f_dg
	(const int side_index[2],               /**< The indices of the affectee, affector, respectively. See the
	                                         *   comments in \ref compute_face_rlhs_dg.h for the convention. */
	 const struct Numerical_Flux* num_flux, ///< Defined for \ref compute_rlhs_fptr.
	 struct DG_Solver_Face* dg_s_face,      ///< Defined for \ref compute_rlhs_fptr.
	 struct Solver_Storage_Implicit* ssi    ///< Defined for \ref compute_rlhs_fptr.
	);

/// \brief Finalize the 2nd order lhs term contribution from the \ref Face for the dg scheme.
static void finalize_lhs_2_f_dg
	(const int side_index[2],                     ///< Defined for \ref finalize_lhs_1_f_dg.
	 const struct Numerical_Flux*const num_flux,  ///< Defined for \ref compute_rlhs_fptr.
	 const struct DG_Solver_Face*const dg_s_face, ///< Defined for \ref compute_rlhs_fptr.
	 struct Solver_Storage_Implicit*const ssi,    ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation*const sim            ///< \ref Simulation.
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

static void scale_by_Jacobian_i2
	(const struct Numerical_Flux*const num_flux, struct Face*const face, const struct Simulation*const sim)
{
UNUSED(sim);
	assert((!face->boundary && num_flux->neigh_info[1].dnnf_dg != NULL) ||
	       ( face->boundary && num_flux->neigh_info[1].dnnf_dg == NULL));

	struct Solver_Face* s_face = (struct Solver_Face*)face;

	const struct const_Vector_d jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->nnf,&jacobian_det_fc,false);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->neigh_info[0].dnnf_dg,&jacobian_det_fc,false);
	if (!face->boundary)
		scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->neigh_info[1].dnnf_dg,&jacobian_det_fc,false);
}

static void compute_rlhs_1
	(const struct Numerical_Flux* num_flux, struct DG_Solver_Face* dg_s_face,
	 struct Solver_Storage_Implicit* ssi, const struct Simulation* sim)
{
	const struct Face* face = (struct Face*) dg_s_face;

	/// See \ref compute_face_rlhs_dg.h for the `lhs_**` notation.
	compute_rhs_f_dg(num_flux,dg_s_face,ssi,sim);

	/** \note The procedure use below is inefficient in the sense that some of the computed operators in \ref
	 *        finalize_lhs_1_f_dg are identical for lhs_ll/lhs_lr to those of lhs_rl/lhs_rr after permutation.
	 *        However, this operators are computing using matrix-vector operations and operators resulting from the
	 *        matrix-matrix products are not duplicated. */
	finalize_lhs_1_f_dg((int[]){0,0},num_flux,dg_s_face,ssi); // lhs_ll
	if (!face->boundary) {
		finalize_lhs_1_f_dg((int[]){0,1},num_flux,dg_s_face,ssi); // lhs_lr

		for (int i = 0; i < 2; ++i) {
			const struct Neigh_Info_NF* n_i = &num_flux->neigh_info[i];
			permute_Multiarray_d_fc((struct Multiarray_d*)n_i->dnnf_ds,'R',1,(struct Solver_Face*)face);
			scale_Multiarray_d((struct Multiarray_d*)n_i->dnnf_ds,-1.0); // Use "-ve" normal.
		}

		finalize_lhs_1_f_dg((int[]){1,0},num_flux,dg_s_face,ssi); // lhs_rl
		finalize_lhs_1_f_dg((int[]){1,1},num_flux,dg_s_face,ssi); // lhs_rr
	}
}

static void compute_rlhs_2
	(const struct Numerical_Flux* num_flux, struct DG_Solver_Face* dg_s_face,
	 struct Solver_Storage_Implicit* ssi, const struct Simulation* sim)
{
	const struct Face*const face = (struct Face*) dg_s_face;

	compute_rhs_f_dg(num_flux,dg_s_face,ssi,sim);

	finalize_lhs_2_f_dg((int[]){0,0},num_flux,dg_s_face,ssi,sim); // lhs_ll (and lhs_lr if not on a boundary).
	if (!face->boundary) {
		finalize_lhs_2_f_dg((int[]){0,1},num_flux,dg_s_face,ssi,sim); // lhs_lr and lhs_ll

		for (int i = 0; i < 2; ++i) {
			const struct Neigh_Info_NF* n_i = &num_flux->neigh_info[i];
			permute_Multiarray_d_fc((struct Multiarray_d*)n_i->dnnf_dg,'R',1,(struct Solver_Face*)face);
			scale_Multiarray_d((struct Multiarray_d*)n_i->dnnf_dg,-1.0); // Use "-ve" normal.
		}

		finalize_lhs_2_f_dg((int[]){1,0},num_flux,dg_s_face,ssi,sim); // lhs_rl and lhs_rr
		finalize_lhs_2_f_dg((int[]){1,1},num_flux,dg_s_face,ssi,sim); // lhs_rl and lhs_rr
//EXIT_UNSUPPORTED;
	}
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for the 'r'ight 'p'artial contribution (d_g_coef__d_s_coef) to the lhs matrix for the gradient
 *         and solution coefficient components corresponding to the input side indices.
 *  \return See brief. */
static const struct const_Matrix_d* constructor_lhs_p_r_gs
	(const int side_index_g,                      ///< Side index of the solution gradients.
	 const int side_index_s,                      ///< Side index of the solution.
	 const struct DG_Solver_Face*const dg_s_face, ///< \ref DG_Solver_Face_T.
	 const struct Simulation*const sim            ///< \ref Simulation.
	);

static void finalize_lhs_1_f_dg
	(const int side_index[2], const struct Numerical_Flux* num_flux, struct DG_Solver_Face* dg_s_face,
	 struct Solver_Storage_Implicit* ssi)
{
	struct Face* face              = (struct Face*) dg_s_face;
	struct Solver_Face* s_face     = (struct Solver_Face*) face;

	struct Matrix_d* lhs = constructor_lhs_f_1(side_index,num_flux,s_face); // destructed

	struct Solver_Volume* s_vol[2] = { (struct Solver_Volume*) face->neigh_info[0].volume,
	                                   (struct Solver_Volume*) face->neigh_info[1].volume, };
	set_petsc_Mat_row_col(ssi,s_vol[side_index[0]],0,s_vol[side_index[1]],0);
	add_to_petsc_Mat(ssi,(struct const_Matrix_d*)lhs);

	destructor_Matrix_d(lhs);
}

static void finalize_lhs_2_f_dg
	(const int side_index[2], const struct Numerical_Flux*const num_flux, const struct DG_Solver_Face*const dg_s_face,
	 struct Solver_Storage_Implicit*const ssi, const struct Simulation*const sim)
{
	struct Face* face                  = (struct Face*) dg_s_face;
	struct Solver_Face* s_face         = (struct Solver_Face*) face;
	struct Solver_Volume*const s_vol[] = { (struct Solver_Volume*) face->neigh_info[0].volume,
	                                       (struct Solver_Volume*) face->neigh_info[1].volume, };

	const struct const_Matrix_d*const lhs_p_l =
		(struct const_Matrix_d*)constructor_lhs_p_f_2(side_index,num_flux,s_face); // destructed
	const struct const_Matrix_d* lhs_p_r_i = NULL;
	const struct const_Matrix_d* lhs_i     = NULL;

	lhs_p_r_i = constructor_lhs_p_r_gs(side_index[1],0,dg_s_face,sim);            // destructed
	lhs_i     = constructor_mm_const_Matrix_d('N','N',1.0,lhs_p_l,lhs_p_r_i,'R'); // destructed
	destructor_const_Matrix_d(lhs_p_r_i);

#if 0
printf("%d %d\n",((struct Volume*)s_vol[side_index[0]])->index,((struct Volume*)s_vol[0])->index);
//print_const_Matrix_d(lhs_p_l);
//print_const_Matrix_d(lhs_p_r_i);
print_const_Matrix_d(lhs_i);
//EXIT_UNSUPPORTED;
#endif

	set_petsc_Mat_row_col(ssi,s_vol[side_index[0]],0,s_vol[0],0);
	add_to_petsc_Mat(ssi,lhs_i);
	destructor_const_Matrix_d(lhs_i);

	if (!face->boundary) {
		lhs_p_r_i = constructor_lhs_p_r_gs(side_index[1],1,dg_s_face,sim);            // destructed
		lhs_i     = constructor_mm_const_Matrix_d('N','N',1.0,lhs_p_l,lhs_p_r_i,'R'); // destructed
		destructor_const_Matrix_d(lhs_p_r_i);

#if 0
printf("%d %d\n",((struct Volume*)s_vol[side_index[0]])->index,((struct Volume*)s_vol[1])->index);
print_const_Matrix_d(lhs_i);
#endif
		set_petsc_Mat_row_col(ssi,s_vol[side_index[0]],0,s_vol[1],0);
		add_to_petsc_Mat(ssi,lhs_i);
		destructor_const_Matrix_d(lhs_i);
	}

	destructor_const_Matrix_d(lhs_p_l);
}

// Level 2 ********************************************************************************************************** //

static const struct const_Matrix_d* constructor_lhs_p_r_gs
	(const int side_index_g, const int side_index_s, const struct DG_Solver_Face*const dg_s_face,
	 const struct Simulation*const sim)
{
	const struct Face*const face = (struct Face*) dg_s_face;
	const struct Solver_Volume*const s_vol_g  = (struct Solver_Volume*) face->neigh_info[side_index_g].volume,
	                          *const s_vol_s  = (struct Solver_Volume*) face->neigh_info[side_index_s].volume;
	const struct DG_Solver_Volume* dg_s_vol_g = (struct DG_Solver_Volume*) s_vol_g;

	const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
	const int n_vr = test_case->n_var;

	const ptrdiff_t n_dof_s = s_vol_s->sol_coef->extents[0],
	                n_dof_g = s_vol_g->grad_coef->extents[0];

	struct Matrix_d*const lhs_p_r = constructor_zero_Matrix_d('R',n_dof_g*n_vr*DIM,n_vr*n_dof_s); // returned

	if (side_index_g == side_index_s) {
		assert(dg_s_vol_g->d_g_coef_v__d_s_coef[0]->ext_0 == n_dof_g);
		assert(dg_s_vol_g->d_g_coef_v__d_s_coef[0]->ext_1 == n_dof_s);
		add_to_lhs_p_r(1.0,dg_s_vol_g->d_g_coef_v__d_s_coef,lhs_p_r,false,sim);
	}
//print_Matrix_d(lhs_p_r);

	const double s = compute_scaling_weak_gradient(dg_s_face,test_case);
	const int s_ind_g = side_index_g,
	          s_ind_s = side_index_s;

	assert(dg_s_face->neigh_info[s_ind_g].d_g_coef_f__d_s_coef[s_ind_s][0]->ext_0 == n_dof_g);
	assert(dg_s_face->neigh_info[s_ind_g].d_g_coef_f__d_s_coef[s_ind_s][0]->ext_1 == n_dof_s);
	add_to_lhs_p_r(s,dg_s_face->neigh_info[s_ind_g].d_g_coef_f__d_s_coef[s_ind_s],lhs_p_r,face->boundary,sim);
//printf("%f\n",s);
//print_Matrix_d(lhs_p_r);
//EXIT_UNSUPPORTED;

	return (struct const_Matrix_d*) lhs_p_r;
}
