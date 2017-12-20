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

#include "compute_all_rlhs_dpg.h"

#include "face_solver_dpg.h"
#include "volume_solver_dpg.h"
#include "element_solver_dpg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "compute_rlhs.h"
#include "compute_face_rlhs.h"
#include "compute_volume_rlhs.h"
#include "flux.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "numerical_flux.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

#include "test_complex_computational_elements.h"
#include "test_complex_compute_all_rhs_dpg.h"
#include "test_complex_test_case.h"

// Static function declarations ************************************************************************************* //

struct Norm_DPG;

/** \brief Constructor for the Jacobian of the lhs volume term with respect to the solution coefficients of 1st order
 *         equations only.
 *  \return See additional comments.
 *
 *  The term returned from this function, henceforth denoted by dLds, will eventually be transposed and multiplied on
 *  the right by the solution coefficient related rows of the RHS vector. To perform the multiplication using a single
 *  BLAS 2 call (i.e. in the most efficient manner), the dLds term is stored in memory with the following constraints:
 *  - ext_0 = n_dof_t*n_eq (size of RHS(s_coef), where n_dof_t == number of test basis functions);
 *  - ext_1 = (n_dof_s*n_vr)^2 (number of columns of standard LHS matrix times number of components to linearize with
 *                              respect to);
 *  - layout = 'C' (Required to ensure that memory of each slice of the 3-tensor is contiguous).
 *
 *  This function is nearly identical to \ref constructor_lhs_v_1_T, except that the \ref Solver_Element::cv0_vs_vc
 *  operator is replaced by the \ref DPG_Solver_Element::cvcv0_vs_vc operator.
 */
static struct Matrix_d* constructor_dlhs_ds_v_1
	(const struct Flux_Ref* flux_r,             ///< \ref Flux_Ref_T.
	 const struct DPG_Solver_Volume* dpg_s_vol, ///< \ref DPG_Solver_Volume_T.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

/** \brief Adds the linearization of \ref Norm_DPG::N to the Hessian of the residual (rhs) wrt the solution coefficients
 *         if the norm has a dependence on the solution. */
static void add_to_dL_ds__norm
	(struct Matrix_d*const dL_ds,              /**< The term storing the linearization of the LHS wrt to the
	                                            *   solution coefficients. */
	 const struct Norm_DPG* norm,              ///< \ref Norm_DPG.
	 const struct const_Matrix_d* optimal_test ///< The optimal test functions.
	);

/** \brief Adds the linearization of the face boundary terms to the Hessian of the residual (rhs) wrt the solution
 *         coefficients. */
static void add_to_dL_ds__face_boundary
	(struct Matrix_d*const dL_ds,               ///< Defined for \ref add_to_dL_ds__norm.
	 const struct DPG_Solver_Volume* dpg_s_vol, ///< Defined for \ref add_to_rlhs__face_boundary.
	 const struct Simulation* sim,              ///< Defined for \ref add_to_rlhs__face_boundary.
	 const char lin_method                      ///< Linearization method. Options: 'a'nalytical, 'c'omplex step.
	);

/** \brief Add entries from the current volume to Solver_Storage_Implicit::A and Solver_Storage_Implicit::b.
 *  \attention **When using the schur complement method to solve the global system, the block diagonal contributions are
 *             inverted before being added to global system matrix.**
 */
static void add_to_petsc_Mat_Vec_dpg
	(const struct Solver_Volume* s_vol,    ///< The current volume.
	 const struct const_Vector_d* rhs_neg, ///< The 'neg'ated local 'r'ight-'h'and 's'ide vector.
	 const struct const_Matrix_d* lhs,     ///< The local 'l'eft-'h'and 's'ide matrix.
	 struct Solver_Storage_Implicit* ssi,  ///< \ref Solver_Storage_Implicit.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_all_rlhs_dpg_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Get the pointer to the appropriate \ref DPG_Solver_Element::cvcv0_vs_vc operator.
 *  \return See brief. */
const struct Operator* get_operator__cvcv0_vs_vc
	(const struct DPG_Solver_Volume* dpg_s_vol ///< The current volume.
	);

/// \brief Version of \ref add_to_dL_ds__face_boundary using the complex step method.
static void add_to_dL_ds__face_boundary_cmplx_step
	(struct Matrix_d*const dL_ds,               ///< See brief.
	 const struct DPG_Solver_Volume* dpg_s_vol, ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

static struct Matrix_d* constructor_dlhs_ds_v_1
	(const struct Flux_Ref* flux_r, const struct DPG_Solver_Volume* dpg_s_vol, const struct Simulation* sim)
{
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) dpg_s_vol;

	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc(s_vol);
	const struct Operator* cvcv0_vs_vc         = get_operator__cvcv0_vs_vc(dpg_s_vol);

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	const ptrdiff_t n_dof_t  = tw1_vt_vc.data[0]->op_std->ext_0,
	                n_dof_s  = s_vol->sol_coef->extents[0],
	                n_dof_s2 = n_dof_s*n_dof_s,
	                n_eq     = test_case->n_eq,
	                n_vr     = test_case->n_var,
	                n_vr2    = n_vr*n_vr;

	const ptrdiff_t ext_1 = tw1_vt_vc.data[0]->op_std->ext_1;
	struct Matrix_d* tw1_r  = constructor_empty_Matrix_d('R',n_dof_t,ext_1);               // destructed
	struct Matrix_d* lhs_l  = constructor_empty_Matrix_d('C',n_dof_t,n_dof_s2);            // destructed
	struct Matrix_d* lhs    = constructor_empty_Matrix_d('C',n_eq*n_dof_t,n_vr2*n_dof_s2); // returned

	const struct const_Multiarray_d* d2fr_ds2_Ma = flux_r->d2fr_ds2;
	struct Vector_d d2fr_ds2 = { .ext_0 = d2fr_ds2_Ma->extents[0], .owns_data = false, .data = NULL, };

	/** \note: Many entries in the Hessian are 0 and need not be included in the loop below.
	 *  \todo Include a list of non-zero entries in the Hessian term (as part of \ref Flux_Ref_T) and skip
	 *        unnecessary terms below. Change to constructor_zero_\*. Also do this for the Jacobian terms. */
	for (int vr2 = 0; vr2 < n_vr; ++vr2) {
	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_d(tw1_r,0.0);
		for (int dim = 0; dim < DIM; ++dim) {
			const ptrdiff_t ind = compute_index_sub_container(
				d2fr_ds2_Ma->order,1,d2fr_ds2_Ma->extents,(ptrdiff_t[]){eq,vr,vr2,dim});
			d2fr_ds2.data = (double*)&d2fr_ds2_Ma->data[ind];
			mm_diag_d('R',1.0,1.0,tw1_vt_vc.data[dim]->op_std,(struct const_Vector_d*)&d2fr_ds2,tw1_r,false);
		}

		mm_d('N','N',1.0,0.0,(struct const_Matrix_d*)tw1_r,cvcv0_vs_vc->op_std,lhs_l);
//print_Matrix_d(lhs_l);

		for (int dof_s = 0; dof_s < n_dof_s; ++dof_s) {
			const ptrdiff_t row = eq*n_dof_t,
			                col = (n_dof_s*n_vr*vr2+dof_s*n_vr+vr)*n_dof_s;
			const double*const data = (const double*) get_col_Matrix_d(dof_s*n_dof_s,lhs_l);
			const struct const_Matrix_d lhs_ll =
				{ .layout = lhs_l->layout, .ext_0 = n_dof_t, .ext_1 = n_dof_s, .data = data, };
			set_block_Matrix_d(lhs,&lhs_ll,row,col,'i');
		}
	}}}
	destructor_Matrix_d(tw1_r);
	destructor_Matrix_d(lhs_l);

	return lhs;
}

static void add_to_petsc_Mat_Vec_dpg
	(const struct Solver_Volume* s_vol, const struct const_Vector_d* rhs_neg, const struct const_Matrix_d* lhs,
	 struct Solver_Storage_Implicit* ssi, const struct Simulation* sim)
{
	assert(sizeof(int) == sizeof(PetscInt));
	assert(sizeof(double) == sizeof(PetscScalar));

	const ptrdiff_t ext_0 = rhs_neg->ext_0;

	const struct const_Vector_i* idxm = constructor_petsc_idxm_dpg(ext_0,s_vol); // destructed.

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	if (test_case->use_schur_complement) {
		const ptrdiff_t dof_s = compute_size(s_vol->sol_coef->order,s_vol->sol_coef->extents),
		                dof_g = compute_size(s_vol->grad_coef->order,s_vol->grad_coef->extents);
		invert_sub_block_Matrix_d((struct Matrix_d*)lhs,0,0,dof_s);
		assert(dof_g == 0); // Add support.
	}

	MatSetValues(ssi->A,(PetscInt)ext_0,idxm->data,(PetscInt)ext_0,idxm->data,lhs->data,ADD_VALUES);
	VecSetValues(ssi->b,(PetscInt)ext_0,idxm->data,rhs_neg->data,ADD_VALUES);

	destructor_const_Vector_i(idxm);
}

static void add_to_dL_ds__norm
	(struct Matrix_d*const dL_ds, const struct Norm_DPG* norm, const struct const_Matrix_d* optimal_test)
{
	if (!norm->dN_ds)
		return;

UNUSED(dL_ds);
UNUSED(optimal_test);
EXIT_ADD_SUPPORT;
}

static void add_to_dL_ds__face_boundary
	(struct Matrix_d*const dL_ds, const struct DPG_Solver_Volume* dpg_s_vol, const struct Simulation* sim,
	 const char lin_method)
{
	if (lin_method == 'a')
		EXIT_ADD_SUPPORT; // Will require linearization of boundary condition/numerical flux Jacobians.
	else if (lin_method == 'c')
		add_to_dL_ds__face_boundary_cmplx_step(dL_ds,dpg_s_vol,sim);
	else
		EXIT_ERROR("Unsupported: %c\n",lin_method);
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for a list of volumes including only a copy of the current volume.
 *  \return See brief. */
static struct Intrusive_List* constructor_Volumes_dpg_local
	(const struct DPG_Solver_Volume*const dpg_s_vol, ///< The current \ref DPG_Solver_Volume_T.
	 struct Simulation*const sim                     ///< The complex \ref Simulation.
	);

/** \brief Constructor for a list of faces including only copies of boundary faces neighbouring the current volume.
 *  \return See brief. */
static struct Intrusive_List* constructor_Faces_dpg_local
	(const struct DPG_Solver_Volume*const dpg_s_vol, ///< The current \ref DPG_Solver_Volume_T.
	 struct Simulation*const sim                     ///< The complex \ref Simulation.
	);

/** \brief Copy members of real to complex computational elements of the current type as determined by the name of the
 *         \ref Simulation::volumes list. */
static void copy_members_computational_elements_dpg
	(const struct DPG_Solver_Volume*const dpg_s_vol, ///< The current \ref DPG_Solver_Volume_T.
	 const struct Simulation*const sim               ///< The complex \ref Simulation.
	);

const struct Operator* get_operator__cvcv0_vs_vc (const struct DPG_Solver_Volume* dpg_s_vol)
{
	const struct Volume* vol           = (struct Volume*) dpg_s_vol;
	const struct Solver_Volume* s_vol  = (struct Solver_Volume*) dpg_s_vol;
	const struct DPG_Solver_Element* e = (struct DPG_Solver_Element*) vol->element;

	const int p = s_vol->p_ref,
	          curved = vol->curved;
	return get_Multiarray_Operator(e->cvcv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});
}

static void add_to_dL_ds__face_boundary_cmplx_step
	(struct Matrix_d*const dL_ds, const struct DPG_Solver_Volume* dpg_s_vol, const struct Simulation* sim)
{
/// \todo Move the simulation constructor out of the hot path!
	struct Simulation* sim_c = constructor_Simulation__no_mesh (sim->ctrl_name); // destructed
	convert_to_Test_Case_rc(sim_c,'c');

	// To avoid recomputing the derived \ref DPG_Solver_Element operators for each volume, the existing
	// \ref Simulation::elements list is used for the complex \ref Simulation as well.
	sim_c->elements = sim->elements;

	assert(sim_c->volumes == NULL && sim_c->faces == NULL);
	sim_c->volumes = constructor_Volumes_dpg_local(dpg_s_vol,sim_c);
	sim_c->faces   = constructor_Faces_dpg_local(dpg_s_vol,sim_c);

	constructor_derived_computational_elements_c(sim_c,IL_SOLVER);
	copy_members_computational_elements_dpg(dpg_s_vol,sim_c);

	constructor_derived_computational_elements_c(sim_c,IL_SOLVER_DPG);

// construct list of volumes as only the current volume, list of faces as only neighbouring boundary faces.
// derive and copy all the way down to DPG_Solver_*.



	const struct Solver_Volume* s_vol =  (struct Solver_Volume*) dpg_s_vol;
	struct Multiarray_d*const s_coef = s_vol->sol_coef;

	const ptrdiff_t size_s = compute_size(s_coef->order,s_coef->extents);

	struct Matrix_d* lhs = constructor_default_Matrix_d();
	lhs->layout    = dL_ds->layout;
	lhs->ext_0     = dL_ds->ext_0;
	lhs->ext_1     = size_s;
	lhs->owns_data = false;
	for (int i = 0; i < size_s; ++i) {
// complex step perturbation

		lhs->data = get_col_Matrix_d(i*size_s,dL_ds);
// possible transpose of lhs
//		add_to_rlhs__face_c(NULL,&
	}
	destructor_Matrix_d(lhs);

	destructor_derived_computational_elements_c(sim_c,IL_SOLVER);
	destructor_derived_computational_elements_c(sim_c,IL_BASE);
	convert_to_Test_Case_rc(sim_c,'r');
	sim_c->elements = NULL;
	destructor_Simulation(sim_c);
/// \todo Move the simulation constructor out of the hot path!
EXIT_UNSUPPORTED;
}

// Level 2 ********************************************************************************************************** //

static struct Intrusive_List* constructor_Volumes_dpg_local
	(const struct DPG_Solver_Volume*const dpg_s_vol, struct Simulation*const sim)
{
	struct Intrusive_List* volumes = constructor_empty_IL(IL_VOLUME,NULL); // returned

	push_back_IL(volumes,(struct Intrusive_Link*)constructor_copy_Volume((struct Volume*)dpg_s_vol,sim,false));
	sim->n_v = 1;

	return volumes;
}

static struct Intrusive_List* constructor_Faces_dpg_local
	(const struct DPG_Solver_Volume*const dpg_s_vol, struct Simulation*const sim)
{
	struct Intrusive_List* faces = constructor_empty_IL(IL_FACE,NULL); // returned

	const struct Volume*const vol_master = (struct Volume*) dpg_s_vol;

	struct Volume* vol = (struct Volume*) sim->volumes->first;
	assert(sim->volumes->first->next == NULL);


	ptrdiff_t n_f = 0;
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face_master = vol_master->faces[i][j];
		if (!face_master || !face_master->boundary)
			continue;

		push_back_IL(faces,(struct Intrusive_Link*)constructor_copy_Face(face_master,sim,false));
		++n_f;

		struct Face* face = (struct Face*) faces->last;

		const_cast_Face(&vol->faces[i][j],face);
		face->neigh_info[0].volume = vol;
	}}
	sim->n_f = n_f;

	return faces;
}

static void copy_members_computational_elements_dpg
	(const struct DPG_Solver_Volume*const dpg_s_vol, const struct Simulation*const sim)
{
	int category = -1;
	if (sim->volumes->name == IL_VOLUME_SOLVER) {
		assert(sim->faces->name == IL_FACE_SOLVER);
		category = IL_SOLVER;
	} else if (sim->volumes->name == IL_VOLUME_SOLVER_DPG) {
		assert(sim->faces->name == IL_FACE_SOLVER_DPG);
		category = IL_SOLVER_DPG;
	} else {
		EXIT_ERROR("Unsupported: %d\n",sim->volumes->name);
	}

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		if (category == IL_SOLVER)
			copy_members_r_to_c_Solver_Volume((struct Solver_Volume_c*)curr,(struct Solver_Volume*)dpg_s_vol);
		else if (category == IL_SOLVER_DPG)
			copy_members_r_to_c_DPG_Solver_Volume((struct DPG_Solver_Volume_c*)curr,dpg_s_vol);
	}
}
