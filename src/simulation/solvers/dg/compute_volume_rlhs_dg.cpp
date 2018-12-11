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

#include "compute_volume_rlhs_dg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "face_solver_dg.h"
#include "volume_solver_dg.h"
#include "element_solver_dg.h"

#include "compute_rlhs.h"
#include "compute_volume_rlhs.h"
#include "flux.h"
#include "intrusive.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "solve_dg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Version of \ref compute_rlhs_v_fptr_T computing the rhs and lhs terms for 1st order equations only.
static void compute_rlhs_1
	(const struct Flux_Ref*const flux_r,      ///< See brief.
	 struct Solver_Volume*const s_vol,        ///< See brief.
	 struct Solver_Storage_Implicit*const ssi ///< See brief.
	);

/// \brief Version of \ref compute_rlhs_v_fptr_T computing the rhs and lhs terms for 2nd order equations only.
static void compute_rlhs_2
	(const struct Flux_Ref*const flux_r,      ///< See brief.
	 struct Solver_Volume*const s_vol,        ///< See brief.
	 struct Solver_Storage_Implicit*const ssi ///< See brief.
	);

/// \brief Version of \ref compute_rlhs_v_fptr_T computing the rhs and lhs terms for both 1st and 2nd order equations.
static void compute_rlhs_12
	(const struct Flux_Ref*const flux_r,      ///< See brief.
	 struct Solver_Volume*const s_vol,        ///< See brief.
	 struct Solver_Storage_Implicit*const ssi ///< See brief.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_volume_rlhs_dg_T.cpp"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "compute_volume_rlhs_dg_T.cpp"
#include "undef_templates_type.h"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_rlhs_v_fptr_T computing the lhs term for linearization wrt the solution.
static void compute_lhs_1
	(const struct Flux_Ref*const flux_r,      ///< See brief.
	 struct DG_Solver_Volume*const dg_s_vol,  ///< See brief.
	 struct Solver_Storage_Implicit*const ssi ///< See brief.
	);

/** \brief Version of \ref compute_rlhs_v_fptr_T computing the lhs term for linearization wrt the solution through the
 *         solution gradients. */
static void compute_lhs_2
	(const struct Flux_Ref*const flux_r,      ///< See brief.
	 struct DG_Solver_Volume*const dg_s_vol,  ///< See brief.
	 struct Solver_Storage_Implicit*const ssi ///< See brief.
	);

static void compute_rlhs_1
	(const struct Flux_Ref*const flux_r, struct Solver_Volume*const s_vol, struct Solver_Storage_Implicit*const ssi)
{
	struct DG_Solver_Volume*const dg_s_vol = (struct DG_Solver_Volume*) s_vol;
	compute_rhs_v_dg_like(flux_r,s_vol,ssi);
	compute_lhs_1(flux_r,dg_s_vol,ssi);
}

static void compute_rlhs_2
	(const struct Flux_Ref*const flux_r, struct Solver_Volume*const s_vol, struct Solver_Storage_Implicit*const ssi)
{
	struct DG_Solver_Volume*const dg_s_vol = (struct DG_Solver_Volume*) s_vol;
	compute_rhs_v_dg_like(flux_r,s_vol,ssi);
	compute_lhs_2(flux_r,dg_s_vol,ssi);
}

static void compute_rlhs_12
	(const struct Flux_Ref*const flux_r, struct Solver_Volume*const s_vol, struct Solver_Storage_Implicit*const ssi)
{
	struct DG_Solver_Volume*const dg_s_vol = (struct DG_Solver_Volume*) s_vol;
	compute_rhs_v_dg_like(flux_r,s_vol,ssi);
	compute_lhs_1(flux_r,dg_s_vol,ssi);
	compute_lhs_2(flux_r,dg_s_vol,ssi);
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for the 'r'ight 'p'artial contribution (d_g_coef__d_s_coef) to the lhs matrix for the 'i'nternal
 *         \ref Volume.
 *  \return See brief.
 *  \note 'i'nternal here refers to being the internal volume wrt to the faces.
 */
static const struct const_Matrix_d* constructor_lhs_p_r_internal
	(const struct DG_Solver_Volume*const dg_s_vol ///< \ref DG_Solver_Volume_T.
	);

/** \brief Add off-diagonal contributions from the linearization of the volume term wrt to the solution through the
 *         gradients to the Petsc Mat. */
static void add_to_petsc_Mat_offdiagonal_volume_2
	(const struct const_Matrix_d*const lhs_p_l, ///< The return of \ref constructor_lhs_p_r_internal.
	 struct DG_Solver_Volume*const dg_s_vol_i,  ///< The 'i'nternal volume.
	 struct Solver_Storage_Implicit*const ssi   ///< Standard.
	);

static void compute_lhs_1
	(const struct Flux_Ref*const flux_r, struct DG_Solver_Volume*const dg_s_vol,
	 struct Solver_Storage_Implicit*const ssi)
{
/// \todo Add special case for collocated.
// If collocation is enabled, note that the diagonal weight scaling must be added back in to recover the symmetry of the
// residual Jacobian. Add it just before adding the contribution to the petsc mat. Also add for face terms and RHS
// terms (volume, face, source or simply the complete rhs).
assert(get_set_collocated(NULL) == false);

	struct Solver_Volume* s_vol = (struct Solver_Volume*) dg_s_vol;
	struct Matrix_d*const lhs = constructor_lhs_v_1(flux_r,s_vol); // destructed
	set_petsc_Mat_row_col_dg(ssi,s_vol,0,s_vol,0);
	add_to_petsc_Mat(ssi,(struct const_Matrix_d*)lhs);
	destructor_Matrix_d(lhs);
}

static void compute_lhs_2
	(const struct Flux_Ref*const flux_r, struct DG_Solver_Volume*const dg_s_vol,
	 struct Solver_Storage_Implicit*const ssi)
{
	// lhs (dependence on sol through grad through the flux)
	struct Solver_Volume* s_vol = (struct Solver_Volume*) dg_s_vol;
	const struct const_Matrix_d*const lhs_p_l =
		(struct const_Matrix_d*) constructor_lhs_p_v_2(flux_r,s_vol); // destructed

	const struct const_Matrix_d*const lhs_p_r_i = constructor_lhs_p_r_internal(dg_s_vol); // destructed

	const struct const_Matrix_d*const lhs_i =
		constructor_mm_const_Matrix_d('N','N',1.0,lhs_p_l,lhs_p_r_i,'R'); // destructed
	destructor_const_Matrix_d(lhs_p_r_i);

#if 0 // OK
print_const_Matrix_d(lhs_i);
#endif

	set_petsc_Mat_row_col_dg(ssi,s_vol,0,s_vol,0);
	add_to_petsc_Mat(ssi,lhs_i);
	destructor_const_Matrix_d(lhs_i);

	/* Note: As the local gradient depends on the solution in neighbouring elements, off-diagonal terms are also
	 *       present when the volume has a face which is not on a boundary (always the case for a mesh having more
	 *       than a single element). */
	add_to_petsc_Mat_offdiagonal_volume_2(lhs_p_l,dg_s_vol,ssi);

	destructor_const_Matrix_d(lhs_p_l);
}

// Level 2 ********************************************************************************************************** //

static const struct const_Matrix_d* constructor_lhs_p_r_internal (const struct DG_Solver_Volume*const dg_s_vol)
{
	struct Solver_Volume* s_vol = (struct Solver_Volume*) dg_s_vol;

	const int n_vr = get_set_n_var_eq(NULL)[0];

	const ptrdiff_t n_dof_s = s_vol->sol_coef->extents[0],
	                n_dof_g = s_vol->grad_coef->extents[0];

	struct Matrix_d*const lhs_p_r = constructor_zero_Matrix_d('R',n_dof_g*n_vr*DIM,n_vr*n_dof_s); // returned

	assert(dg_s_vol->d_g_coef_v__d_s_coef[0]->ext_0 == n_dof_g);
	assert(dg_s_vol->d_g_coef_v__d_s_coef[0]->ext_1 == n_dof_s);

	add_to_lhs_p_r(1.0,dg_s_vol->d_g_coef_v__d_s_coef,lhs_p_r,false);

	const struct Volume*const vol = (struct Volume*) dg_s_vol;
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face = vol->faces[i][j];
		if (!face)
			continue;

		const struct DG_Solver_Face*const dg_s_face = (struct DG_Solver_Face*) face;
		const int s_ind = compute_side_index_face(face,vol);

		const int mult = ( !face->boundary ? 1 : n_vr );
		assert(dg_s_face->neigh_info[s_ind].d_g_coef_f__d_s_coef[s_ind][0]->ext_0 == n_dof_g*mult);
		assert(dg_s_face->neigh_info[s_ind].d_g_coef_f__d_s_coef[s_ind][0]->ext_1 == n_dof_s*mult);
		add_to_lhs_p_r(1.0,dg_s_face->neigh_info[s_ind].d_g_coef_f__d_s_coef[s_ind],lhs_p_r,face->boundary);
	}}

	return (struct const_Matrix_d*) lhs_p_r;
}

static void add_to_petsc_Mat_offdiagonal_volume_2
	(const struct const_Matrix_d*const lhs_p_l, struct DG_Solver_Volume*const dg_s_vol_i,
	 struct Solver_Storage_Implicit*const ssi)
{
	const int n_vr = get_set_n_var_eq(NULL)[0];

	const struct Volume*const vol_i          = (struct Volume*) dg_s_vol_i;
	const struct Solver_Volume*const s_vol_i = (struct Solver_Volume*) dg_s_vol_i;
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face = vol_i->faces[i][j];
		if (!face || face->boundary)
			continue;

		const int s_ind_i = compute_side_index_face(face,vol_i),
		          s_ind_o = ( s_ind_i == 0 ? 1 : 0 );
		const struct DG_Solver_Volume*const dg_s_vol_o = (struct DG_Solver_Volume*)face->neigh_info[s_ind_o].volume;
		const struct Solver_Volume*const s_vol_o = (struct Solver_Volume*) dg_s_vol_o;

		const ptrdiff_t n_dof_s = s_vol_o->sol_coef->extents[0],
				    n_dof_g = s_vol_i->grad_coef->extents[0];

		struct Matrix_d*const lhs_p_r = constructor_zero_Matrix_d('R',n_dof_g*n_vr*DIM,n_vr*n_dof_s); // destructed

		const struct DG_Solver_Face*const dg_s_face = (struct DG_Solver_Face*) face;
		assert(dg_s_face->neigh_info[s_ind_i].d_g_coef_f__d_s_coef[s_ind_o][0]->ext_0 == n_dof_g);
		assert(dg_s_face->neigh_info[s_ind_i].d_g_coef_f__d_s_coef[s_ind_o][0]->ext_1 == n_dof_s);
		add_to_lhs_p_r(1.0,dg_s_face->neigh_info[s_ind_i].d_g_coef_f__d_s_coef[s_ind_o],lhs_p_r,face->boundary);

		const struct const_Matrix_d*const lhs_o =
			constructor_mm_const_Matrix_d('N','N',1.0,lhs_p_l,(struct const_Matrix_d*)lhs_p_r,'R'); // destructed
		destructor_Matrix_d(lhs_p_r);
#if 0 // OK
printf("face: %d %d %d\n",((struct Volume*)s_vol_i)->index,((struct Volume*)s_vol_o)->index,s_ind_i);
print_const_Matrix_d(lhs_o);
#endif

		set_petsc_Mat_row_col_dg(ssi,s_vol_i,0,s_vol_o,0);
		add_to_petsc_Mat(ssi,lhs_o);
		destructor_const_Matrix_d(lhs_o);
	}}
}
