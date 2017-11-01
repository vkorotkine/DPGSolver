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

#include "solve_dg.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "petscvec.h"

#include "macros.h"
#include "definitions_test_case.h"

#include "face.h"
#include "face_solver.h"
#include "element.h"
#include "element_solver_dg.h"
#include "volume.h"
#include "volume_solver_dg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "compute_grad_coef_dg.h"
#include "compute_volume_rlhs_dg.h"
#include "compute_face_rlhs_dg.h"
#include "compute_source_rlhs_dg.h"
#include "intrusive.h"
#include "math_functions.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Call the common functions when computing the rlhs values for explicit and implicit dg schemes.
static void compute_rlhs_common_dg
	(const struct Simulation* sim,             ///< \ref Simulation.
	 struct Solver_Storage_Implicit* s_store_i ///< \ref Solver_Storage_Implicit.
	);

/// \brief Scale the rhs terms by the the inverse mass matrices (for explicit schemes).
static void scale_rhs_by_m_inv
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Compute the maximum value of the rhs term for the first variable.
 *  \return See brief. */
static double compute_max_rhs
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Fill \ref Solver_Storage_Implicit::b with the negative of rhs values.
static void fill_petsc_Vec_b_dg
	(const struct Simulation* sim,             ///< Defined for \ref compute_rlhs_dg.
	 struct Solver_Storage_Implicit* s_store_i ///< Defined for \ref compute_rlhs_dg.
	);

// Interface functions ********************************************************************************************** //

double compute_rhs_dg (const struct Simulation* sim)
{
	compute_rlhs_common_dg(sim,NULL);
	scale_rhs_by_m_inv(sim);

	return compute_max_rhs(sim);
}

double compute_rlhs_dg (const struct Simulation* sim, struct Solver_Storage_Implicit* s_store_i)
{
	compute_rlhs_common_dg(sim,s_store_i);
	fill_petsc_Vec_b_dg(sim,s_store_i);

	return compute_max_rhs(sim);
}

void permute_Multiarray_d_fc
	(struct Multiarray_d* data, const char perm_layout, const int side_index_dest, const struct Solver_Face* s_face)
{
	const struct Neigh_Info* neigh_info = &((struct Face*)s_face)->neigh_info[side_index_dest];

	struct Volume* vol = neigh_info->volume;
	const struct const_DG_Solver_Element* f_e =
		(const struct const_DG_Solver_Element*) get_element_by_face(vol->element,neigh_info->ind_lf);

	const int ind_ord = neigh_info->ind_ord,
	          p_f     = s_face->p_ref;
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	const struct const_Vector_i* nc_fc =
		get_const_Multiarray_Vector_i(f_e->nc_fc[curved],(ptrdiff_t[]){ind_ord,0,0,p_f,p_f});

	permute_Multiarray_d_V(data,nc_fc,perm_layout);
}

void permute_Matrix_d_fc
	(struct Matrix_d* data, const char perm_layout, const int side_index_dest, const struct Solver_Face* s_face)
{
	assert(perm_layout == 'R');
	assert(data->layout == 'R');

	const struct Neigh_Info* neigh_info = &((struct Face*)s_face)->neigh_info[side_index_dest];

	struct Volume* vol = neigh_info->volume;
	const struct const_DG_Solver_Element* f_e =
		(const struct const_DG_Solver_Element*) get_element_by_face(vol->element,neigh_info->ind_lf);

	const int ind_ord = neigh_info->ind_ord,
	          p_f     = s_face->p_ref;
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	const struct const_Vector_i* nc_fc =
		get_const_Multiarray_Vector_i(f_e->nc_fc[curved],(ptrdiff_t[]){ind_ord,0,0,p_f,p_f});

	permute_Matrix_d_V(data,nc_fc);
}

void set_petsc_Mat_row_col
	(struct Solver_Storage_Implicit*const s_store_i, const struct Solver_Volume* v_l, const int eq,
	 const struct Solver_Volume* v_r, const int vr)
{
	s_store_i->row = v_l->ind_dof+v_l->sol_coef->extents[0]*eq;
	s_store_i->col = v_r->ind_dof+v_r->sol_coef->extents[0]*vr;
}

void add_to_petsc_Mat (const struct Solver_Storage_Implicit*const s_store_i, const struct const_Matrix_d* lhs)
{
	const ptrdiff_t ext_0 = lhs->ext_0,
	                ext_1 = lhs->ext_1;

	PetscInt idxm[ext_0],
	         idxn[ext_1];

	for (int i = 0; i < ext_0; ++i)
		idxm[i] = s_store_i->row+i;

	for (int i = 0; i < ext_1; ++i)
		idxn[i] = s_store_i->col+i;

	const PetscScalar*const vv = lhs->data;
	MatSetValues(s_store_i->A,ext_0,idxm,ext_1,idxn,vv,ADD_VALUES);
}

// Static functions ************************************************************************************************* //
// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Set the memory of the rhs and lhs (if applicable) terms to zero for the volumes.
static void zero_memory_volumes
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Version of \ref scale_rhs_by_m_inv used for non-collocated runs.
static void scale_rhs_by_m_inv_std
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Version of \ref scale_rhs_by_m_inv used for collocated runs.
static void scale_rhs_by_m_inv_col
	(const struct Simulation* sim ///< \ref Simulation.
	);

static void compute_rlhs_common_dg (const struct Simulation* sim, struct Solver_Storage_Implicit* s_store_i)
{
	zero_memory_volumes(sim);
	compute_grad_coef_dg(sim);
	compute_volume_rlhs_dg(sim,s_store_i);
	compute_face_rlhs_dg(sim,s_store_i);
	compute_source_rhs_dg(sim);
}

static void scale_rhs_by_m_inv (const struct Simulation* sim)
{
	assert((sim->test_case->solver_proc == SOLVER_E) ||
	       (sim->test_case->solver_proc == SOLVER_EI));

	if (!sim->collocated)
		scale_rhs_by_m_inv_std(sim);
	else
		scale_rhs_by_m_inv_col(sim);
}

static double compute_max_rhs (const struct Simulation* sim)
{
	double max_rhs = 0.0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) curr;

		struct Multiarray_d* rhs = dg_s_vol->rhs;
		double max_rhs_curr = norm_d(rhs->extents[0],rhs->data,"Inf");
		if (max_rhs_curr > max_rhs)
			max_rhs = max_rhs_curr;
	}
	return max_rhs;
}

static void fill_petsc_Vec_b_dg (const struct Simulation* sim, struct Solver_Storage_Implicit* s_store_i)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const int ind_dof        = ((struct Solver_Volume*)curr)->ind_dof;
		struct Multiarray_d* rhs = ((struct DG_Solver_Volume*)curr)->rhs;

		const int ni = compute_size(rhs->order,rhs->extents);

		PetscInt    ix[ni];
		PetscScalar y[ni];

		for (int i = 0; i < ni; ++i) {
			ix[i] = ind_dof+i;
			y[i]  = -(rhs->data[i]);
		}

		VecSetValues(s_store_i->b,ni,ix,y,INSERT_VALUES);
	}
}

// Level 1 ********************************************************************************************************** //

static void zero_memory_volumes (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		set_to_value_Multiarray_d(((struct DG_Solver_Volume*)curr)->rhs,0.0);
	}
}

static void scale_rhs_by_m_inv_std (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) curr;
		mm_NN1C_overwrite_Multiarray_d(dg_s_vol->m_inv,&dg_s_vol->rhs);
	}
}

static void scale_rhs_by_m_inv_col (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol       = (struct Solver_Volume*) curr;
		struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) curr;

		const struct const_Vector_d jac_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
		scale_Multiarray_by_Vector_d('L',1.0,dg_s_vol->rhs,&jac_det_vc,true);
	}
}
