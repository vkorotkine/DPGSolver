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

#include <assert.h>
#include <stdlib.h>
#include "petscmat.h"

#include "macros.h"
#include "definitions_dpg.h"
#include "definitions_elements.h"
#include "definitions_intrusive.h"


#include "def_templates_compute_all_rlhs_dpg.h"

#include "def_templates_face_solver_dpg.h"
#include "def_templates_face_solver.h"
#include "def_templates_volume_solver_dpg.h"
#include "def_templates_volume_solver.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_compute_face_rlhs.h"
#include "def_templates_compute_volume_rlhs.h"
#include "def_templates_flux.h"
#include "def_templates_operators.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

#define USE_EXACT_NORMAL_FLUX false ///< Flag for whether the exact normal flux should be used on the boundary faces.

struct S_Params_DPG;
struct Norm_DPG;

/** \brief Function pointer to functions constructing the \ref Norm_DPG whose members are used to evaluate the optimal
 *         test functions and their Jacobians.
 *
 *  \param dpg_s_vol The current volume.
 *  \param flux_r    \ref Flux_Ref_T.
 *  \param sim       \ref Simulation.
 */
typedef const struct Norm_DPG* (*constructor_norm_DPG_fptr)
	(const struct DPG_Solver_Volume_T* dpg_s_vol,
	 const struct Flux_Ref_T* flux_r,
	 const struct Simulation* sim
	);

/** \brief Pointer to functions computing rhs and lhs terms for the dpg solver.
 *
 *  \param s_params  \ref S_Params_DPG.
 *  \param flux_i    \ref Flux_Input_T.
 *  \param dpg_s_vol Pointer to the current volume.
 *  \param ssi       \ref Solver_Storage_Implicit.
 *  \param sim       \ref Simulation.
 *  \param sim_c     The complex \ref Simulation (may be NULL).
 */
typedef void (*compute_rlhs_dpg_fptr)
	(const struct S_Params_DPG* s_params,
	 struct Flux_Input_T* flux_i,
	 const struct DPG_Solver_Volume_T* dpg_s_vol,
	 struct Solver_Storage_Implicit* ssi,
	 const struct Simulation*const sim,
	 struct Simulation*const sim_c
	);

/// \brief Container for solver-related parameters.
struct S_Params_DPG {
	struct S_Params_Volume_Structor_T spvs; ///< \ref S_Params_Volume_Structor_T.

	constructor_norm_DPG_fptr constructor_norm_DPG; ///< Pointer to the appropriate function.
	compute_rlhs_dpg_fptr compute_rlhs;             ///< Pointer to the appropriate function.
};

/// \brief Container for DPG norm-related parameters.
struct Norm_DPG {
	const struct const_Matrix_T* N;     ///< The DPG norm operator.
	const struct const_Matrix_T* dN_ds; ///< The Jacobian of the DPG norm operator wrt the solution coefficients.
};

/** \brief Set the parameters of \ref S_Params_DPG.
 *  \return A statically allocated \ref S_Params_DPG container. */
static struct S_Params_DPG set_s_params_dpg
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Set the values of the global indices corresponding to the current unknown.
static void set_idxm
	(int* ind_idxm,                  ///< Pointer to the index in `idxm`.
	 struct Vector_i* idxm,          ///< Vector of global indices.
	 const int ind_dof,              ///< Index of the first dof corresponding to the `coef`.
	 const struct Multiarray_T* coef ///< Pointer to the coefficients under consideration.
	);

/// \brief Provides the internal face contribution to \ref add_to_rlhs__face_T.
static void add_to_rlhs__face_internal
	(const struct DPG_Solver_Volume_T* dpg_s_vol, ///< The current \ref DPG_Solver_Volume_T.
	 const struct DPG_Solver_Face_T* dpg_s_face,  ///< The current \ref DPG_Solver_Face_T.
	 struct Matrix_T* lhs,                        ///< The lhs matrix contribution for the current volume/faces.
	 struct Matrix_T* rhs,                        ///< The rhs matrix contribution for the current volume/faces.
	 int* ind_dof,                                ///< The index of the current dof under consideration.
	 const struct Simulation* sim                 ///< \ref Simulation.
	);

/// \brief Add the current face contribution to \ref Solver_Volume_T::flux_imbalance from both sides.
static void add_to_flux_imbalance
	(const struct Solver_Face_T*const s_face, ///< The current \ref Solver_Face_T.
	 const struct Simulation*const sim        ///< \ref Simulation.
	);

/** \brief Constructor for the \ref Numerical_Flux_T to be used for the dpg boundary condition imposition.
 *  \return See brief. */
static struct Numerical_Flux_T* constructor_Numerical_Flux_dpg
	(const struct DPG_Solver_Face_T*const dpg_s_face, ///< The current face.
	 const struct Simulation*const sim                ///< \ref Simulation.
	);

/// \brief Scale required \ref Numerical_Flux_T terms by the face Jacobian.
static void scale_by_Jacobian
	(const struct Numerical_Flux_T* num_flux, ///< \ref Numerical_Flux_T.
	 const struct Solver_Face_T* s_face       ///< The current \ref Solver_Face_T.
	);

/// \brief Increment the rhs terms with the contribution of the boundary face.
static void increment_rhs_boundary_face
	(struct Matrix_T* rhs,                    ///< Holds the rhs terms.
	 const struct Numerical_Flux_T* num_flux, ///< \ref Numerical_Flux_T.
	 const struct Solver_Face_T* s_face,      ///< The current \ref Solver_Face_T.
	 const struct Simulation* sim             ///< \ref Simulation.
	);

/// \brief Increment the lhs terms with the contribution of the boundary face.
static void increment_lhs_boundary_face
	(struct Matrix_T* lhs,                    ///< Holds the lhs terms.
	 const struct Numerical_Flux_T* num_flux, ///< \ref Numerical_Flux_T.
	 const struct Solver_Face_T* s_face,      ///< The current \ref Solver_Face_T.
	 const struct Simulation* sim             ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void compute_all_rlhs_dpg_T
	(const struct Simulation* sim, struct Solver_Storage_Implicit* ssi, struct Intrusive_List* volumes)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DPG);
	assert(sim->faces->name == IL_FACE_SOLVER_DPG);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DPG);

	// Unsteady cases not currently supported.
	assert(((struct Test_Case_T*)sim->test_case_rc->tc)->solver_method_curr == 'i');

	struct S_Params_DPG s_params = set_s_params_dpg(sim);
	struct Flux_Input_T* flux_i = constructor_Flux_Input_T(sim); // destructed

	/** A complex \ref Simulation is initialized here as it is used to compute Hessian terms relating to the
	 *  boundary conditions using the complex step method. While this is certainly less efficient than computing these
	 *  terms through the analytical Hessian, efficiency may not be significantly impacted as this procedure is only
	 *  used for faces on the boundary. \todo Profile and refer to test.
	 *
	 *  As the boundary conditions are currently imposed weakly through a numerical flux after computing the
	 *  appropriate boundary ghost state (exactly as is done in the DG solver), the analytical Hessian of both
	 *  boundary condition and numerical flux functions would be required if it is desired to remove the complex step
	 *  linearization. */
	struct Simulation* sim_c = NULL;
#if TYPE_RC == TYPE_REAL
	sim_c = constructor_Simulation__no_mesh(sim->ctrl_name); // destructed
	convert_to_Test_Case_rc(sim_c,'c');

	struct Test_Case_c* test_case = (struct Test_Case_c*) sim_c->test_case_rc->tc;
	/** Currently, the function pointer to the boundary condition/numerical flux computing functions are set based on
	 *  the value of solver_method_curr and functions do not necessarily allow for the computation of all required
	 *  terms. For example, setting solver_method_curr = 'e' here, \ref Test_Case_T::flux_comp_mem_e will indicate
	 *  that the Jacobian should be computed, but the function pointer will point to a function which does not support
	 *  the Jacobian computation and the values will all simply be 0.
	 *
	 * \todo Combine the separate (and redundant) numerical flux/boundary condition functions, similarly to what was
	 *       done for the standard flux functions. Once complete, change solver_method_curr to 'e' for the
	 *       Test_Case_c.
	 */
	test_case->solver_method_curr = 'i';

	/* To avoid recomputing the derived \ref DPG_Solver_Element operators for each volume, the existing
	 * \ref Simulation::elements list is used for the complex \ref Simulation as well. */
	sim_c->elements = sim->elements;
#endif

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct DPG_Solver_Volume_T* dpg_s_vol = (struct DPG_Solver_Volume_T*) curr;
		s_params.compute_rlhs(&s_params,flux_i,dpg_s_vol,ssi,sim,sim_c);
	}
	destructor_Flux_Input_T(flux_i);

#if TYPE_RC == TYPE_REAL
	convert_to_Test_Case_rc(sim_c,'r');
	sim_c->elements = NULL;
	destructor_Simulation(sim_c);
#endif
}

void compute_flux_imbalances_faces_dpg_T (struct Simulation*const sim)
{
	assert(list_is_derived_from("solver",'f',sim));
	assert(sim->elements->name == IL_ELEMENT_SOLVER);
	constructor_derived_Elements(sim,IL_ELEMENT_SOLVER_DPG); // destructed

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	test_case->solver_method_curr = 'e';

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) curr;
		add_to_flux_imbalance(s_face,sim);
	}

	destructor_derived_Elements(sim,IL_ELEMENT_SOLVER);
	test_case->solver_method_curr = 0;
}

const struct const_Matrix_T* constructor_lhs_l_internal_face_dpg_T
	(const struct DPG_Solver_Volume_T* dpg_s_vol, const struct DPG_Solver_Face_T* dpg_s_face)
{
	const struct Volume* vol           = (struct Volume*) dpg_s_vol;
	const struct Face* face            = (struct Face*) dpg_s_face;
	const struct Solver_Face_T* s_face = (struct Solver_Face_T*) dpg_s_face;

	const int side_index = compute_side_index_face(face,vol);
	const struct Operator* tw0_vt_fc_op = get_operator__tw0_vt_fc_T(side_index,s_face),
	                     * cv0_ff_fc_op = get_operator__cv0_ff_fc_T(s_face);

	struct Matrix_T* cv0_ff_fc = constructor_copy_Matrix_T_Matrix_R((struct Matrix_R*)cv0_ff_fc_op->op_std); // destructed

	const struct const_Multiarray_T* j_det = s_face->jacobian_det_fc;
	const struct const_Vector_T* j_det_V =
		constructor_copy_const_Vector_T_T(compute_size(j_det->order,j_det->extents),j_det->data); // destructed
	scale_Matrix_by_Vector_T('L',1.0,cv0_ff_fc,j_det_V,false);
	destructor_const_Vector_T(j_det_V);

	// Use "-ve" sign when looking from volume[0] as the face term was moved to the rhs of the equation. When looking
	// from volume[1], the "-ve" sign is cancelled by the "-ve" sign of the inverted normal vector.
	Real alpha = -1.0;
	if (side_index == 1) {
		permute_Matrix_T_fc(cv0_ff_fc,'R',side_index,s_face);
		alpha = 1.0;
	}

	const struct const_Matrix_T* lhs_l = constructor_mm_RT_const_Matrix_T
		('N','N',alpha,tw0_vt_fc_op->op_std,(struct const_Matrix_T*)cv0_ff_fc,'R'); // returned
	destructor_Matrix_T(cv0_ff_fc);

	return lhs_l;
}

ptrdiff_t compute_n_dof_nf_T (const struct Solver_Volume_T* s_vol)
{
	ptrdiff_t dof = 0;

	const struct Volume* vol = (struct Volume*) s_vol;
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face = vol->faces[i][j];
		if (!face)
			continue;

		const struct Solver_Face_T* s_face = (struct Solver_Face_T*) face;
		const ptrdiff_t size = s_face->nf_coef->extents[0];
		dof += size;

		assert((size > 0) || (face->boundary));
	}}
	return dof;
}

const struct const_Vector_i* constructor_petsc_idxm_dpg_T
	(const ptrdiff_t n_dof, const struct Solver_Volume_T* s_vol, const struct Simulation*const sim)
{
	struct Vector_i* idxm = constructor_empty_Vector_i(n_dof); // returned
	int ind_idxm = 0;

	// volume (sol_coef, grad_coef)
	set_idxm(&ind_idxm,idxm,(int)s_vol->ind_dof,s_vol->sol_coef);
	// grad_coef: To be done.

	// face (nf_coef, sol_coef)
	const struct Volume* vol = (struct Volume*) s_vol;
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face = vol->faces[i][j];
		if (!face)
			continue;

		const struct Solver_Face_T* s_face = (struct Solver_Face_T*) face;
		set_idxm(&ind_idxm,idxm,(int)s_face->ind_dof,s_face->nf_coef);
		// sol_coef: To be done.
	}}

	// volume (l_mult) - if applicable
	if (test_case_explicitly_enforces_conservation(sim))
		set_idxm(&ind_idxm,idxm,(int)s_vol->ind_dof_constraint,s_vol->l_mult);

	assert(ind_idxm == n_dof);

	return (struct const_Vector_i*) idxm;
}

void add_to_rlhs__face_T
	(struct Vector_T* rhs, struct Matrix_T** lhs_ptr, const struct DPG_Solver_Volume_T* dpg_s_vol,
	 const struct Simulation*const sim, const bool include_internal)
{
	struct Matrix_T* lhs = *lhs_ptr;

	const struct Volume* vol            = (struct Volume*) dpg_s_vol;
	const struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) dpg_s_vol;

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq,
	          n_vr = test_case->n_var;

	const ptrdiff_t n_dof_s  = (lhs->ext_1)/n_vr,
	                n_dof_nf = compute_n_dof_nf_T(s_vol);

	assert(include_internal == (rhs != NULL));
	struct Matrix_T* rhs_M = NULL;
	if (include_internal) {
		struct Matrix_T* lhs_add = constructor_zero_Matrix_T('R',lhs->ext_0,(n_dof_s+n_dof_nf)*n_vr); // moved
		set_block_Matrix_T(lhs_add,0,0,(struct const_Matrix_T*)lhs,0,0,lhs->ext_0,lhs->ext_1,'i');
		destructor_Matrix_T(lhs);
		lhs = lhs_add;

		rhs_M = constructor_move_Matrix_T_T('C',(rhs->ext_0)/n_eq,n_eq,false,rhs->data); // destructed

		// Internal faces
		int ind_dof = (int)(n_vr*n_dof_s);
		for (int i = 0; i < NFMAX;    ++i) {
		for (int j = 0; j < NSUBFMAX; ++j) {
			const struct Face* face = vol->faces[i][j];
			if (!face || face->boundary)
				continue;

			const struct DPG_Solver_Face_T* dpg_s_face = (struct DPG_Solver_Face_T*) face;
			add_to_rlhs__face_internal(dpg_s_vol,dpg_s_face,lhs,rhs_M,&ind_dof,sim);
		}}
	}

	// Boundary faces
	add_to_rlhs__face_boundary_T(dpg_s_vol,lhs,rhs_M,sim);

	destructor_conditional_Matrix_T(rhs_M);

	*lhs_ptr = lhs;
}

void add_to_rlhs__face_boundary_T
	(const struct DPG_Solver_Volume_T* dpg_s_vol, struct Matrix_T* lhs, struct Matrix_T* rhs,
	 const struct Simulation*const sim)
{

	const struct Volume*const vol = (struct Volume*) dpg_s_vol;
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face = vol->faces[i][j];
		if (!face || !face->boundary)
			continue;

		const struct DPG_Solver_Face_T* dpg_s_face = (struct DPG_Solver_Face_T*) face;
		struct Numerical_Flux_T* num_flux = constructor_Numerical_Flux_dpg(dpg_s_face,sim); // destructed

		const struct Solver_Face_T* s_face = (struct Solver_Face_T*) dpg_s_face;
		scale_by_Jacobian(num_flux,s_face);

		if (rhs != NULL)
			increment_rhs_boundary_face(rhs,num_flux,s_face,sim);
		increment_lhs_boundary_face(lhs,num_flux,s_face,sim);

		destructor_Numerical_Flux_T(num_flux);
	}}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_rlhs_dpg_fptr computing rhs and lhs terms for 1st order equations only.
static void compute_rlhs_1
	(const struct S_Params_DPG* s_params,         ///< See brief.
	 struct Flux_Input_T* flux_i,                 ///< See brief.
	 const struct DPG_Solver_Volume_T* dpg_s_vol, ///< See brief.
	 struct Solver_Storage_Implicit* ssi,         ///< See brief.
	 const struct Simulation*const sim,           ///< See brief.
	 struct Simulation*const sim_c                ///< See brief.
	);

/** \brief Version of \ref constructor_norm_DPG_fptr; see comments for \ref TEST_NORM_H0.
 *  \return See brief. */
static const struct Norm_DPG* constructor_norm_DPG__h0
	(const struct DPG_Solver_Volume_T* dpg_s_vol, ///< See brief.
	 const struct Flux_Ref_T* flux_r,             ///< See brief.
	 const struct Simulation* sim                 ///< See brief.
	);

/** \brief Version of \ref constructor_norm_DPG_fptr; see comments for \ref TEST_NORM_H1.
 *  \return See brief. */
static const struct Norm_DPG* constructor_norm_DPG__h1
	(const struct DPG_Solver_Volume_T* dpg_s_vol, ///< See brief.
	 const struct Flux_Ref_T* flux_r,             ///< See brief.
	 const struct Simulation* sim                 ///< See brief.
	);

/** \brief Version of \ref constructor_norm_DPG_fptr; see comments for \ref TEST_NORM_H1_UPWIND.
 *  \return See brief. */
static const struct Norm_DPG* constructor_norm_DPG__h1_upwind
	(const struct DPG_Solver_Volume_T* dpg_s_vol, ///< See brief.
	 const struct Flux_Ref_T* flux_r,             ///< See brief.
	 const struct Simulation* sim                 ///< See brief.
	);

/// \brief Set the numerical flux members to those corresponding to the use of the exact normal flux values.
static void set_exact_normal_flux
	(const struct Solver_Face_T*const s_face,       ///< \ref Solver_Face_T.
	 struct mutable_Numerical_Flux_T*const num_flux ///< \ref mutable_Numerical_Flux_T.
	);

static struct S_Params_DPG set_s_params_dpg (const struct Simulation* sim)
{
	struct S_Params_DPG s_params;

	set_S_Params_Volume_Structor_T(&s_params.spvs,sim);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	switch (test_case->solver_method_curr) {
	case 'i':
		if (test_case->has_1st_order && !test_case->has_2nd_order)
			s_params.compute_rlhs = compute_rlhs_1;
		else if (!test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_2;
		else if (test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_12;
		else
			EXIT_ERROR("Unsupported: %d %d\n",test_case->has_1st_order,test_case->has_2nd_order);
		break;
#if TYPE_RC == TYPE_REAL
	case 'e': // fallthrough
#endif
	default:
		EXIT_ERROR("Unsupported: %c\n",test_case->solver_method_curr);
		break;
	}

	switch (test_case->ind_test_norm) {
		case TEST_NORM_H0:        s_params.constructor_norm_DPG = constructor_norm_DPG__h0;        break;
		case TEST_NORM_H1:        s_params.constructor_norm_DPG = constructor_norm_DPG__h1;        break;
		case TEST_NORM_H1_UPWIND: s_params.constructor_norm_DPG = constructor_norm_DPG__h1_upwind; break;
		default:                  EXIT_ERROR("Unsupported: %d\n",test_case->ind_test_norm);        break;
	}

	return s_params;
}

static void set_idxm (int* ind_idxm, struct Vector_i* idxm, const int ind_dof, const struct Multiarray_T* coef)
{
	if (coef == NULL)
		return;

	const ptrdiff_t size = compute_size(coef->order,coef->extents);
	for (int i = 0; i < size; ++i) {
		idxm->data[*ind_idxm] = ind_dof+i;
		++*ind_idxm;
	}
}

static void add_to_rlhs__face_internal
	(const struct DPG_Solver_Volume_T* dpg_s_vol, const struct DPG_Solver_Face_T* dpg_s_face, struct Matrix_T* lhs,
	 struct Matrix_T* rhs, int* ind_dof, const struct Simulation* sim)
{
	/// As the rhs is **always** linear wrt the trace unknowns, the rhs and lhs are computed together.
	const struct const_Matrix_T* lhs_l = constructor_lhs_l_internal_face_dpg_T(dpg_s_vol,dpg_s_face); // destructed
	const struct Solver_Face_T* s_face  = (struct Solver_Face_T*) dpg_s_face;
	struct Matrix_T nf_coef = interpret_Multiarray_as_Matrix_T(s_face->nf_coef);
	mm_T('N','N',1.0,1.0,lhs_l,(struct const_Matrix_T*)&nf_coef,rhs);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq,
	          n_vr = test_case->n_var;

	const ptrdiff_t n_dof_test = (lhs->ext_0)/n_eq,
	                n_dof_nf   = nf_coef.ext_0;
	for (int vr = 0; vr < n_vr; ++vr) {
		set_block_Matrix_T(lhs,vr*n_dof_test,*ind_dof,lhs_l,0,0,lhs_l->ext_0,lhs_l->ext_1,'i');
		*ind_dof += (int)n_dof_nf;
	}
	destructor_const_Matrix_T(lhs_l);
}

static void add_to_flux_imbalance (const struct Solver_Face_T*const s_face, const struct Simulation*const sim)
{
	UNUSED(sim);
	const struct Face*const face                    = (struct Face*) s_face;
	const struct DPG_Solver_Face_T*const dpg_s_face = (struct DPG_Solver_Face_T*) s_face;

	const struct const_Vector_R* w_fc = get_operator__w_fc__s_e_T(s_face);
	const struct const_Vector_T jacobian_det_fc = interpret_const_Multiarray_as_Vector_T(s_face->jacobian_det_fc);

	if (!face->boundary) {
		const struct Operator* cv0_ff_fc = get_operator__cv0_ff_fc_T(s_face);
		const struct const_Matrix_T nf_coef_M =
			interpret_const_Multiarray_as_Matrix_T((struct const_Multiarray_T*)s_face->nf_coef);
		const struct const_Matrix_T* nf_M =
			constructor_mm_RT_const_Matrix_T('N','N',1.0,cv0_ff_fc->op_std,&nf_coef_M,'C'); // destructed
		scale_Matrix_by_Vector_T('L',1.0,(struct Matrix_T*)nf_M,&jacobian_det_fc,false);
		add_to_flux_imbalance_face_nf_w_T(nf_M,w_fc,s_face);
		destructor_const_Matrix_T(nf_M);
	} else {
		struct Numerical_Flux_T* num_flux = constructor_Numerical_Flux_dpg(dpg_s_face,sim); // destructed
		const struct const_Matrix_T nf_M = interpret_const_Multiarray_as_Matrix_T(num_flux->nnf);
		scale_Matrix_by_Vector_T('L',1.0,(struct Matrix_T*)&nf_M,&jacobian_det_fc,false);
		add_to_flux_imbalance_face_nf_w_T(&nf_M,w_fc,s_face);
		destructor_Numerical_Flux_T(num_flux);
	}
}

static struct Numerical_Flux_T* constructor_Numerical_Flux_dpg
	(const struct DPG_Solver_Face_T*const dpg_s_face, const struct Simulation*const sim)
{
	struct Numerical_Flux_Input_T* num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed

	const struct Solver_Face_T* s_face = (struct Solver_Face_T*) dpg_s_face;
	constructor_Numerical_Flux_Input_data_T(num_flux_i,s_face,sim); // destructed

	struct Numerical_Flux_T* num_flux = constructor_Numerical_Flux_T(num_flux_i); // returned
	destructor_Numerical_Flux_Input_data_T(num_flux_i);
	destructor_Numerical_Flux_Input_T(num_flux_i);

	if (USE_EXACT_NORMAL_FLUX) {
		assert(s_face->nf_fc != NULL);
		set_exact_normal_flux(s_face,(struct mutable_Numerical_Flux_T*)num_flux);
	}
	return num_flux;
}

static void scale_by_Jacobian (const struct Numerical_Flux_T* num_flux, const struct Solver_Face_T* s_face)
{
	const struct Face* face = (struct Face*) s_face;

	assert(face->boundary);
	assert(num_flux->neigh_info[0].dnnf_ds != NULL || num_flux->neigh_info[0].dnnf_dg != NULL);

	const struct const_Vector_T jacobian_det_fc = interpret_const_Multiarray_as_Vector_T(s_face->jacobian_det_fc);
	scale_Multiarray_by_Vector_T('L',1.0,(struct Multiarray_T*)num_flux->nnf,&jacobian_det_fc,false);

	if (num_flux->neigh_info[0].dnnf_ds)
		scale_Multiarray_by_Vector_T(
			'L',1.0,(struct Multiarray_T*)num_flux->neigh_info[0].dnnf_ds,&jacobian_det_fc,false);
	if (num_flux->neigh_info[0].dnnf_dg)
		EXIT_ADD_SUPPORT;
}

static void increment_rhs_boundary_face
	(struct Matrix_T* rhs, const struct Numerical_Flux_T* num_flux, const struct Solver_Face_T* s_face,
	 const struct Simulation* sim)
{
	assert(((struct Face*)s_face)->boundary);

	ptrdiff_t extents[2] = { rhs->ext_0, rhs->ext_1, };
	struct Multiarray_T rhs_Ma =
		{ .layout = 'C', .order = 2, .extents = extents, .owns_data = false, .data = rhs->data, };

	const struct Operator* tw0_vt_fc = get_operator__tw0_vt_fc_T(0,s_face);

	UNUSED(sim);
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';
	mm_NNC_Operator_Multiarray_T(-1.0,1.0,tw0_vt_fc,num_flux->nnf,&rhs_Ma,op_format,2,NULL,NULL);
}

static void increment_lhs_boundary_face
	(struct Matrix_T* lhs, const struct Numerical_Flux_T* num_flux, const struct Solver_Face_T* s_face,
	 const struct Simulation* sim)
{
	UNUSED(sim);
	assert(((struct Face*)s_face)->boundary);

	struct Matrix_T* lhs_ll = constructor_lhs_f_1_T((int[]){0,0},num_flux,s_face); // destructed

	set_block_Matrix_T(lhs,0,0,(struct const_Matrix_T*)lhs_ll,0,0,lhs_ll->ext_0,lhs_ll->ext_1,'a');
	destructor_Matrix_T(lhs_ll);
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for the rhs \ref Vector_T with volume contributions from 1st order equations included.
 *  \return See brief. */
static struct Vector_T* constructor_rhs_v_1
	(const struct Flux_Ref_T* flux_r,     ///< Defined for \ref compute_rlhs_dpg_fptr.
	 const struct Solver_Volume_T* s_vol, ///< Defined for \ref compute_rlhs_dpg_fptr.
	 const struct Simulation* sim         ///< Defined for \ref compute_rlhs_dpg_fptr.
	);

/// \brief Increment the rhs terms with the source contribution.
static void increment_rhs_source
	(struct Vector_T* rhs,                ///< Holds the values of the rhs.
	 const struct Solver_Volume_T* s_vol, ///< Defined for \ref compute_rlhs_dpg_fptr.
	 const struct Simulation* sim         ///< Defined for \ref compute_rlhs_dpg_fptr.
	);

/** \brief Constructor for the negated DPG rhs \ref Vector_T, including contributions from conservation enforcement if
 *         applicable.
 *  \return See brief. */
static const struct const_Vector_T* constructor_rhs_opt_neg
	(const struct const_Vector_T*const rhs_std,        ///< The DG-like component of the rhs vector.
	 const struct const_Matrix_T*const lhs_std,        ///< The DG-like component of the lhs matrix.
	 const struct const_Matrix_T*const opt_t_coef,     ///< The optimal test function coefficients.
	 const struct DPG_Solver_Volume_T*const dpg_s_vol, ///< The current volume.
	 const struct Simulation*const sim                 ///< \ref Simulation.
	);

/// \brief Destructor for \ref Norm_DPG.
static void destructor_Norm_DPG
	(const struct Norm_DPG* norm ///< Standard.
	);

#if TYPE_RC == TYPE_COMPLEX
/** \brief Version of \ref add_to_petsc_Mat_Vec_dpg setting a single column of Solver_Storage_Implicit::A to the values
 *         computed using the complex step rhs. */
static void add_to_petsc_Mat_dpg_c
	(const struct Solver_Volume_c* s_vol,  ///< See brief.
	 const struct const_Vector_c* rhs_neg, ///< See brief.
	 struct Solver_Storage_Implicit* ssi,  ///< See brief.
	 const struct Simulation*const sim     ///< See brief.
		);
#endif

static const struct Norm_DPG* constructor_norm_DPG__h0
	(const struct DPG_Solver_Volume_T* dpg_s_vol, const struct Flux_Ref_T* flux_r, const struct Simulation* sim)
{
	UNUSED(flux_r);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq;

	const struct const_Matrix_T* norm_op_H0 = dpg_s_vol->norm_op_H0;
	const ptrdiff_t ext_0 = norm_op_H0->ext_0;

	struct Matrix_T* N = constructor_zero_Matrix_T('R',n_eq*ext_0,n_eq*ext_0); // moved

	for (int eq = 0; eq < n_eq; ++eq)
		set_block_Matrix_T(N,eq*ext_0,eq*ext_0,norm_op_H0,0,0,norm_op_H0->ext_0,norm_op_H0->ext_1,'a');

	struct Norm_DPG* norm = malloc(sizeof* norm); // returned
	norm->N     = (struct const_Matrix_T*) N;
	norm->dN_ds = NULL;

	return norm;
}

static const struct Norm_DPG* constructor_norm_DPG__h1
	(const struct DPG_Solver_Volume_T* dpg_s_vol, const struct Flux_Ref_T* flux_r, const struct Simulation* sim)
{
	UNUSED(flux_r);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq;

	const struct const_Matrix_T* norm_op_H1 = dpg_s_vol->norm_op_H1;
	const ptrdiff_t ext_0 = norm_op_H1->ext_0;

	struct Matrix_T* N = constructor_zero_Matrix_T('R',n_eq*ext_0,n_eq*ext_0); // moved

	for (int eq = 0; eq < n_eq; ++eq)
		set_block_Matrix_T(N,eq*ext_0,eq*ext_0,norm_op_H1,0,0,norm_op_H1->ext_0,norm_op_H1->ext_1,'a');

	struct Norm_DPG* norm = malloc(sizeof* norm); // returned
	norm->N     = (struct const_Matrix_T*) N;
	norm->dN_ds = NULL;

	return norm;
}

static const struct Norm_DPG* constructor_norm_DPG__h1_upwind
	(const struct DPG_Solver_Volume_T* dpg_s_vol, const struct Flux_Ref_T* flux_r, const struct Simulation* sim)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq,
	          n_vr = test_case->n_var;

	struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) dpg_s_vol;
	const struct Multiarray_Operator cv1_vt_vc = get_operator__cv1_vt_vc_T(s_vol);

	const ptrdiff_t ext_0 = cv1_vt_vc.data[0]->op_std->ext_0,
	                ext_1 = cv1_vt_vc.data[0]->op_std->ext_1;

	struct Matrix_T* cv1r = constructor_empty_Matrix_T('R',n_vr*ext_0,n_eq*ext_1); // destructed

	struct Matrix_T* cv1r_l = constructor_empty_Matrix_T('R',ext_0,ext_1); // destructed
	const struct const_Multiarray_T* dfr_ds_Ma = flux_r->dfr_ds;
	struct Vector_T dfr_ds = { .ext_0 = dfr_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_T(cv1r_l,0.0);
		for (int dim = 0; dim < DIM; ++dim) {
			const ptrdiff_t ind =
				compute_index_sub_container(dfr_ds_Ma->order,1,dfr_ds_Ma->extents,(ptrdiff_t[]){eq,vr,dim});
			dfr_ds.data = (Type*)&dfr_ds_Ma->data[ind];
			mm_diag_T('L',1.0,1.0,cv1_vt_vc.data[dim]->op_std,(struct const_Vector_T*)&dfr_ds,cv1r_l,false);
		}
		set_block_Matrix_T(cv1r,eq*ext_0,vr*ext_1,
		                   (struct const_Matrix_T*)cv1r_l,0,0,cv1r_l->ext_0,cv1r_l->ext_1,'i');
	}}
	destructor_Matrix_T(cv1r_l);

	const struct const_Vector_R* w_vc = get_operator__w_vc__s_e_T(s_vol);
	const struct const_Vector_T J_vc  = interpret_const_Multiarray_as_Vector_T(s_vol->jacobian_det_vc);

	const struct const_Vector_T* J_inv_vc = constructor_inverse_const_Vector_T(&J_vc);                   // destructed
	const struct const_Vector_T* wJ_vc    = constructor_dot_mult_const_Vector_T_RT(1.0,w_vc,J_inv_vc,n_vr); // destructed
	destructor_const_Vector_T(J_inv_vc);

	const struct const_Matrix_T* n1_lt =
		constructor_mm_diag_const_Matrix_T(1.0,(struct const_Matrix_T*)cv1r,wJ_vc,'L',false); // destructed
	destructor_const_Vector_T(wJ_vc);

	// norm->N
	const struct const_Matrix_T* n1 =
		constructor_mm_const_Matrix_T('T','N',1.0,n1_lt,(struct const_Matrix_T*)cv1r,'R'); // destructed
	destructor_Matrix_T(cv1r);

	const struct const_Matrix_T* norm_op_H0 = dpg_s_vol->norm_op_H0;
	assert(norm_op_H0->ext_0 == ext_1);

	struct Matrix_T* N = constructor_empty_Matrix_T('R',n_eq*ext_1,n_eq*ext_1); // moved

	set_block_Matrix_T(N,0,0,n1,0,0,n1->ext_0,n1->ext_1,'i');
	for (int eq = 0; eq < n_eq; ++eq)
		set_block_Matrix_T(N,eq*ext_1,eq*ext_1,norm_op_H0,0,0,norm_op_H0->ext_0,norm_op_H0->ext_1,'a');
	destructor_const_Matrix_T(n1);

	// norm->dN_ds
	const struct const_Matrix_T* dN_ds = NULL;
#if TYPE_RC == TYPE_REAL
	dN_ds = constructor_norm_DPG_dN_ds__h1_upwind(dpg_s_vol,flux_r,n1_lt,sim); // moved
#endif
	destructor_const_Matrix_T(n1_lt);

	struct Norm_DPG* norm = malloc(sizeof* norm); // returned
	norm->N     = (struct const_Matrix_T*) N;
	norm->dN_ds = dN_ds;

	return norm;
}

static void compute_rlhs_1
	(const struct S_Params_DPG* s_params, struct Flux_Input_T* flux_i, const struct DPG_Solver_Volume_T* dpg_s_vol,
	 struct Solver_Storage_Implicit* ssi, const struct Simulation*const sim, struct Simulation*const sim_c)
{
	const struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) dpg_s_vol;

	struct Flux_Ref_T* flux_r = constructor_Flux_Ref_vol_T(&s_params->spvs,flux_i,s_vol); // destructed

	const struct Norm_DPG* norm = s_params->constructor_norm_DPG(dpg_s_vol,flux_r,sim); // destructed

	struct Vector_T* rhs_std = constructor_rhs_v_1(flux_r,s_vol,sim); // destructed
	struct Matrix_T* lhs_std = constructor_lhs_v_1_T(flux_r,s_vol);   // destructed

	add_to_rlhs__face_T(rhs_std,&lhs_std,dpg_s_vol,sim,true);
	increment_rhs_source(rhs_std,s_vol,sim);

	const struct const_Matrix_T* opt_t_coef =
		constructor_sysv_const_Matrix_T(norm->N,(struct const_Matrix_T*)lhs_std); // destructed

	const struct const_Vector_T*const rhs_opt_neg =
		constructor_rhs_opt_neg((struct const_Vector_T*)rhs_std,(struct const_Matrix_T*)lhs_std,
		                        opt_t_coef,dpg_s_vol,sim); // destructed

#if TYPE_RC == TYPE_REAL
	struct Matrix_T* lhs_opt =
		constructor_mm_Matrix_T('T','N',1.0,opt_t_coef,(struct const_Matrix_T*)lhs_std,'R'); // destructed

	add_to_lhs_opt__d_opt_t_ds(flux_r,dpg_s_vol,norm,opt_t_coef,(struct const_Vector_T*)rhs_std,
	                           (struct const_Matrix_T*)lhs_std,lhs_opt,sim,sim_c);
	add_to_lhs_opt__l_mult(&lhs_opt,(struct const_Matrix_T*)lhs_std,dpg_s_vol,sim,sim_c);

	add_to_petsc_Mat_Vec_dpg(s_vol,rhs_opt_neg,(struct const_Matrix_T*)lhs_opt,ssi,sim);
	destructor_Matrix_T(lhs_opt);
#elif TYPE_RC == TYPE_COMPLEX
	UNUSED(sim_c);

	add_to_petsc_Mat_dpg_c(s_vol,rhs_opt_neg,ssi,sim);
#endif
	destructor_Flux_Ref_T(flux_r);
	destructor_Norm_DPG(norm);

	destructor_Vector_T(rhs_std);
	destructor_Matrix_T(lhs_std);
	destructor_const_Vector_T(rhs_opt_neg);

	destructor_const_Matrix_T(opt_t_coef);
}

static void set_exact_normal_flux
	(const struct Solver_Face_T*const s_face, struct mutable_Numerical_Flux_T*const num_flux)
{
	assert(((struct Face*)s_face)->boundary == true);

	destructor_Multiarray_T(num_flux->nnf);
	num_flux->nnf = constructor_copy_Multiarray_T((struct Multiarray_T*)s_face->nf_fc);

	struct m_Neigh_Info_NF_T n_i = num_flux->neigh_info[0];
	if (n_i.dnnf_ds)
		set_to_value_Multiarray_T(n_i.dnnf_ds,0.0);
	if (n_i.dnnf_dg)
		set_to_value_Multiarray_T(n_i.dnnf_dg,0.0);
}

// Level 2 ********************************************************************************************************** //

/** \brief Get the appropriate sub-range of the \ref DPG_Solver_Element::ones_coef_vt operators.
 *  \return See brief. */
static const struct const_Vector_R* get_operator__ones_coef_vt
	(const struct DPG_Solver_Volume_T*const dpg_s_vol ///< The current volume.
	);

/** \brief Constructor for the transposed linearized contribution of the Lagrange multiplier constraint equations for
 *         conservation represented as a \ref const_Matrix_T\* with one column per Lagrange multiplier unknown.
 *  \return See brief. */
static const struct const_Matrix_T* constructor_l_mult_M
	(const struct const_Matrix_T*const lhs_std,        ///< Defined for \ref constructor_rhs_opt_neg.
	 const struct DPG_Solver_Volume_T*const dpg_s_vol, ///< Defined for \ref constructor_rhs_opt_neg.
	 const struct Simulation*const sim                 ///< Defined for \ref constructor_rhs_opt_neg.
	);

static struct Vector_T* constructor_rhs_v_1
	(const struct Flux_Ref_T* flux_r, const struct Solver_Volume_T* s_vol, const struct Simulation* sim)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq;

	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc_T(s_vol);

	const ptrdiff_t ext_0 = tw1_vt_vc.data[0]->op_std->ext_0;

	struct Vector_T* rhs = constructor_zero_Vector_T(ext_0*n_eq); // returned

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	ptrdiff_t extents[2] = { ext_0, n_eq, };
	struct Multiarray_T rhs_Ma =
		{ .layout = 'C', .order = 2, .extents = extents, .owns_data = false, .data = rhs->data, };
	for (ptrdiff_t dim = 0; dim < DIM; ++dim)
		mm_NNC_Operator_Multiarray_T(1.0,1.0,tw1_vt_vc.data[dim],flux_r->fr,&rhs_Ma,op_format,2,&dim,NULL);

	return rhs;
}

static void increment_rhs_source
	(struct Vector_T* rhs, const struct Solver_Volume_T* s_vol, const struct Simulation* sim)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq;

	ptrdiff_t extents[2] = { (rhs->ext_0)/n_eq, n_eq, };
	struct Multiarray_T rhs_Ma =
		{ .layout = 'C', .order = 2, .extents = extents, .owns_data = rhs->owns_data, .data = rhs->data };

	test_case->compute_source_rhs(sim,s_vol,&rhs_Ma);
}

static const struct const_Vector_T* constructor_rhs_opt_neg
	(const struct const_Vector_T*const rhs_std, const struct const_Matrix_T*const lhs_std,
	 const struct const_Matrix_T*const opt_t_coef, const struct DPG_Solver_Volume_T*const dpg_s_vol,
	 const struct Simulation*const sim)
{
	const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) dpg_s_vol;
	struct Vector_T*const rhs_opt_neg = constructor_mv_Vector_T('T',-1.0,opt_t_coef,rhs_std); // returned

	if (test_case_explicitly_enforces_conservation(sim)) {
		struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
		const int n_eq = test_case->n_eq;
		assert(test_case->ind_conservation == CONSERVATION_LAGRANGE_MULT);

		// Constraint effect on original "form".
		const struct const_Matrix_T*const lhs_l_mult_M = constructor_l_mult_M(lhs_std,dpg_s_vol,sim); // destructed

		const struct Multiarray_T*const l_mult = s_vol->l_mult;
		const struct const_Vector_T l_mult_V =
			{ .ext_0 = l_mult->extents[0], .owns_data = false, .data = l_mult->data };
		const struct const_Vector_T*const rhs_l_mult_V =
			constructor_mv_const_Vector_T('N',-1.0,lhs_l_mult_M,&l_mult_V); // destructed
		destructor_const_Matrix_T(lhs_l_mult_M);

		add_to_Vector_T(rhs_opt_neg,rhs_l_mult_V);
		destructor_const_Vector_T(rhs_l_mult_V);

		// Independent constraint equation.
		const struct const_Matrix_T rhs_std_M =
			{ .layout = 'C', .ext_0 = (rhs_std->ext_0)/n_eq, .ext_1 = n_eq, .owns_data = false,
			  .data = rhs_std->data, };

		const struct const_Vector_R*const ones_coef = get_operator__ones_coef_vt(dpg_s_vol);
		const struct const_Vector_T* ones_coef_T = constructor_copy_const_Vector_T_Vector_R(ones_coef); // dest.
		const struct const_Vector_T*const rhs_constraint =
			constructor_mv_const_Vector_T('T',-1.0,&rhs_std_M,ones_coef_T); // destructed
		destructor_const_Vector_T(ones_coef_T);

		push_back_Vector_Vector_T(rhs_opt_neg,rhs_constraint);
		destructor_const_Vector_T(rhs_constraint);
	}

	return (struct const_Vector_T*) rhs_opt_neg;
}

static void destructor_Norm_DPG (const struct Norm_DPG* norm)
{
	destructor_const_Matrix_T(norm->N);
/// \todo Make conditional destructor.
	if (norm->dN_ds)
		destructor_const_Matrix_T(norm->dN_ds);
	free((void*)norm);
}

#if TYPE_RC == TYPE_COMPLEX
static void add_to_petsc_Mat_dpg_c
	(const struct Solver_Volume_c* s_vol, const struct const_Vector_c* rhs_neg, struct Solver_Storage_Implicit* ssi,
	 const struct Simulation*const sim)
{
	const ptrdiff_t ext_0 = rhs_neg->ext_0;

	const struct const_Vector_i* idxm = constructor_petsc_idxm_dpg_c(ext_0,s_vol,sim); // destructed.

	PetscScalar rhs_c_data[ext_0];
	for (int i = 0; i < ext_0; ++i)
		rhs_c_data[i] = cimag((-rhs_neg->data[i])/CX_STEP);

	MatSetValues(ssi->A,(PetscInt)ext_0,idxm->data,1,&ssi->col,rhs_c_data,ADD_VALUES);

	destructor_const_Vector_i(idxm);
}
#endif

// Level 3 ********************************************************************************************************** //

static const struct const_Vector_R* get_operator__ones_coef_vt (const struct DPG_Solver_Volume_T*const dpg_s_vol)
{
	const struct Volume*const vol                 = (struct Volume*) dpg_s_vol;
	const struct Solver_Volume_T*const s_vol      = (struct Solver_Volume_T*) dpg_s_vol;
	const struct DPG_Solver_Element*const dpg_s_e = (struct DPG_Solver_Element*) vol->element;

	const int p = s_vol->p_ref;

	return get_const_Multiarray_Vector_d(dpg_s_e->ones_coef_vt,(ptrdiff_t[]){0,0,p,p});
}

static const struct const_Matrix_T* constructor_l_mult_M
	(const struct const_Matrix_T*const lhs_std, const struct DPG_Solver_Volume_T*const dpg_s_vol,
	 const struct Simulation*const sim)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq;
	assert(test_case->ind_conservation == CONSERVATION_LAGRANGE_MULT);

	const bool transpose = ( lhs_std->layout == 'C' ? false : true );
	if (transpose)
		transpose_Matrix_T((struct Matrix_T*)lhs_std,true);

	ptrdiff_t exts_lhs[2] = { (lhs_std->ext_0)/n_eq, (lhs_std->ext_1)*n_eq, };
	const struct const_Matrix_T lhs_std_M =
		{ .layout = 'C', .ext_0 = exts_lhs[0], .ext_1 = exts_lhs[1], .owns_data = false, .data = lhs_std->data, };

	const struct const_Vector_R*const ones_coef = get_operator__ones_coef_vt(dpg_s_vol);
	const struct const_Vector_T* ones_coef_T = constructor_copy_const_Vector_T_Vector_R(ones_coef); // dest.

	struct Vector_T*const lhs_l_mult = constructor_mv_Vector_T('T',1.0,&lhs_std_M,ones_coef_T); // destructed
	lhs_l_mult->owns_data = false;
	const struct const_Matrix_T*const lhs_l_mult_M =
		constructor_move_const_Matrix_T_T('R',(lhs_l_mult->ext_0)/n_eq,n_eq,true,lhs_l_mult->data); // returned
	destructor_Vector_T(lhs_l_mult);

	if (transpose)
		transpose_Matrix_T((struct Matrix_T*)lhs_std,true);

	return lhs_l_mult_M;
}

#include "undef_templates_compute_all_rlhs_dpg.h"

#include "undef_templates_face_solver_dpg.h"
#include "undef_templates_face_solver.h"
#include "undef_templates_volume_solver_dpg.h"
#include "undef_templates_volume_solver.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_compute_face_rlhs.h"
#include "undef_templates_compute_volume_rlhs.h"
#include "undef_templates_flux.h"
#include "undef_templates_operators.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_test_case.h"
