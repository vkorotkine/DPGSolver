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

#include "def_templates_compute_all_rlhs_dpg_d.h"
#include "def_templates_compute_face_rlhs.h"
#include "def_templates_compute_volume_rlhs.h"
#include "def_templates_flux.h"
#include "def_templates_multiarray_operator_d.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

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
typedef void (*compute_rlhs_fptr)
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
	compute_rlhs_fptr compute_rlhs;                 ///< Pointer to the appropriate function.
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

/** \brief Get the pointer to the appropriate \ref DPG_Solver_Element::cv0_ff_fc operator.
 *  \return See brief. */
static const struct Operator* get_operator__cv0_ff_fc
	(const int side_index,                    ///< The index of the side of the face under consideration.
	 const struct DPG_Solver_Face_T* dpg_s_face ///< The current face.
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

/// \brief Provides the boundary face contribution to \ref add_to_rlhs__face_T.
static void add_to_rlhs__face_boundary
	(const struct DPG_Solver_Volume_T* dpg_s_vol, ///< The current \ref DPG_Solver_Volume_T.
	 const struct DPG_Solver_Face_T* dpg_s_face,  ///< The current \ref DPG_Solver_Face_T.
	 struct Matrix_T* lhs,                        ///< The lhs matrix contribution for the current volume/faces.
	 struct Matrix_T* rhs,                        ///< The rhs matrix contribution for the current volume/faces.
	 const struct Simulation*const sim            ///< \ref Simulation.
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
//struct Volume*        vol   = (struct Volume*) curr;
//printf("v_ind: %d\n",vol->index);
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

struct Multiarray_Operator get_operator__cvt1_vt_vc__rlhs_T (const struct DPG_Solver_Volume_T* dpg_s_vol)
{
	struct Volume* vol            = (struct Volume*) dpg_s_vol;
	struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) vol;

	const struct DPG_Solver_Element* dpg_s_e = (struct DPG_Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;

	return set_MO_from_MO(dpg_s_e->cvt1_vt_vc[curved],1,(ptrdiff_t[]){0,0,p,p});
}

const struct const_Matrix_R* constructor_lhs_l_internal_face_dpg_T
	(const struct DPG_Solver_Volume_T* dpg_s_vol, const struct DPG_Solver_Face_T* dpg_s_face)
{
	const struct Volume* vol           = (struct Volume*) dpg_s_vol;
	const struct Face* face            = (struct Face*) dpg_s_face;
	const struct Solver_Face_T* s_face = (struct Solver_Face_T*) dpg_s_face;

	const int side_index = compute_side_index_face(face,vol);
//printf("%d %d %d\n",vol->index,face->index,side_index);
	const struct Operator* tw0_vt_fc_op = get_operator__tw0_vt_fc_T(side_index,s_face),
	                     * cv0_ff_fc_op = get_operator__cv0_ff_fc(side_index,dpg_s_face);

	struct Matrix_R* cv0_ff_fc = constructor_copy_Matrix_R((struct Matrix_R*)cv0_ff_fc_op->op_std); // destructed

	const struct const_Multiarray_R* j_det = s_face->jacobian_det_fc;
	const struct const_Vector_R* j_det_V =
		constructor_copy_const_Vector_R_R(compute_size(j_det->order,j_det->extents),j_det->data); // destructed
	scale_Matrix_R_by_Vector_R('L',1.0,cv0_ff_fc,j_det_V,false);
	destructor_const_Vector_R(j_det_V);

	// Use "-ve" sign when looking from volume[0] as the face term was moved to the rhs of the equation. When looking
	// from volume[1], the "-ve" sign is cancelled by the "-ve" sign of the inverted normal vector.
	Real alpha = -1.0;
	if (side_index == 1) {
		permute_Matrix_R_fc(cv0_ff_fc,'R',side_index,s_face);
		alpha = 1.0;
	}

	const struct const_Matrix_R* lhs_l = constructor_mm_const_Matrix_R
		('N','N',alpha,tw0_vt_fc_op->op_std,(struct const_Matrix_R*)cv0_ff_fc,'R'); // returned
	destructor_Matrix_R(cv0_ff_fc);

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
	(const ptrdiff_t n_dof, const struct Solver_Volume_T* s_vol)
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
		set_block_Matrix_T(lhs_add,(struct const_Matrix_T*)lhs,0,0,'i');
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
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face = vol->faces[i][j];
		if (!face || !face->boundary)
			continue;

		const struct DPG_Solver_Face_T* dpg_s_face = (struct DPG_Solver_Face_T*) face;
		add_to_rlhs__face_boundary(dpg_s_vol,dpg_s_face,lhs,rhs_M,sim);
	}}

	destructor_conditional_Matrix_T(rhs_M);

//print_Matrix_T(lhs);
	*lhs_ptr = lhs;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_rlhs_fptr computing rhs and lhs terms for 1st order equations only.
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

/** \brief Version of \ref constructor_norm_DPG_fptr; see comments for \ref TEST_NORM_H1_UPWIND.
 *  \return See brief. */
static const struct Norm_DPG* constructor_norm_DPG__h1_upwind
	(const struct DPG_Solver_Volume_T* dpg_s_vol, ///< See brief.
	 const struct Flux_Ref_T* flux_r,             ///< See brief.
	 const struct Simulation* sim                 ///< See brief.
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
		case TEST_NORM_H1_UPWIND: s_params.constructor_norm_DPG = constructor_norm_DPG__h1_upwind; break;
		default:                  EXIT_ERROR("Unsupported: %d\n",test_case->ind_test_norm);        break;
	}

	return s_params;
}

static const struct Operator* get_operator__cv0_ff_fc
	(const int side_index, const struct DPG_Solver_Face_T* dpg_s_face)
{
	const struct Face* face            = (struct Face*) dpg_s_face;
	const struct Solver_Face_T* s_face = (struct Solver_Face_T*) dpg_s_face;
	const struct Volume* vol           = face->neigh_info[side_index].volume;
	const struct DPG_Solver_Element* e = (struct DPG_Solver_Element*) vol->element;

	const int ind_e  = get_face_element_index(face),
	          p_f    = s_face->p_ref;
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );

	return get_Multiarray_Operator(e->cv0_ff_fc[curved],(ptrdiff_t[]){ind_e,ind_e,0,0,p_f,p_f});
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
	const struct const_Matrix_R* lhs_l = constructor_lhs_l_internal_face_dpg_T(dpg_s_vol,dpg_s_face); // destructed
//print_const_Matrix_d(tw0_vt_fc_op->op_std);
	const struct Solver_Face_T* s_face  = (struct Solver_Face_T*) dpg_s_face;
	struct Matrix_T nf_coef = interpret_Multiarray_as_Matrix_T(s_face->nf_coef);
	mm_RTT('N','N',1.0,1.0,lhs_l,(struct const_Matrix_T*)&nf_coef,rhs);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq,
	          n_vr = test_case->n_var;

	const ptrdiff_t n_dof_test = (lhs->ext_0)/n_eq,
	                n_dof_nf   = nf_coef.ext_0;
	for (int vr = 0; vr < n_vr; ++vr) {
		set_block_Matrix_T_R(lhs,lhs_l,vr*n_dof_test,*ind_dof,'i');
		*ind_dof += (int)n_dof_nf;
	}
	destructor_const_Matrix_d(lhs_l);
}

static void add_to_rlhs__face_boundary
	(const struct DPG_Solver_Volume_T* dpg_s_vol, const struct DPG_Solver_Face_T* dpg_s_face, struct Matrix_T* lhs,
	 struct Matrix_T* rhs, const struct Simulation*const sim)
{
	UNUSED(dpg_s_vol);
	struct Numerical_Flux_Input_T* num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed

	const struct Solver_Face_T* s_face = (struct Solver_Face_T*) dpg_s_face;
	constructor_Numerical_Flux_Input_data_T(num_flux_i,s_face,sim); // destructed

	struct Numerical_Flux_T* num_flux = constructor_Numerical_Flux_T(num_flux_i); // destructed
	destructor_Numerical_Flux_Input_data_T(num_flux_i);
	destructor_Numerical_Flux_Input_T(num_flux_i);

	scale_by_Jacobian(num_flux,s_face);

	if (rhs != NULL)
		increment_rhs_boundary_face(rhs,num_flux,s_face,sim);
	increment_lhs_boundary_face(lhs,num_flux,s_face,sim);

	destructor_Numerical_Flux_T(num_flux);
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for the rhs \ref Vector_T with volume contributions from 1st order equations included.
 *  \return See brief. */
static struct Vector_T* constructor_rhs_v_1
	(const struct Flux_Ref_T* flux_r,     ///< Defined for \ref compute_rlhs_fptr.
	 const struct Solver_Volume_T* s_vol, ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim         ///< Defined for \ref compute_rlhs_fptr.
	);

/// \brief Increment the rhs terms with the source contribution.
static void increment_rhs_source
	(struct Vector_T* rhs,                ///< Holds the values of the rhs.
	 const struct Solver_Volume_T* s_vol, ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim         ///< Defined for \ref compute_rlhs_fptr.
	);

/// \brief Destructor for \ref Norm_DPG.
static void destructor_Norm_DPG
	(const struct Norm_DPG* norm ///< Standard.
	);

static const struct Norm_DPG* constructor_norm_DPG__h0
	(const struct DPG_Solver_Volume_T* dpg_s_vol, const struct Flux_Ref_T* flux_r, const struct Simulation* sim)
{
	UNUSED(flux_r);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq;

	const struct const_Matrix_R* norm_op_H0 = dpg_s_vol->norm_op_H0;
	const ptrdiff_t ext_0 = norm_op_H0->ext_0;

	struct Matrix_T* N = constructor_zero_Matrix_T('R',n_eq*ext_0,n_eq*ext_0); // moved

	for (int eq = 0; eq < n_eq; ++eq)
		set_block_Matrix_T_R(N,norm_op_H0,eq*ext_0,eq*ext_0,'a');

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

	const struct Multiarray_Operator cvt1_vt_vc = get_operator__cvt1_vt_vc__rlhs_T(dpg_s_vol);

	const ptrdiff_t ext_0 = cvt1_vt_vc.data[0]->op_std->ext_0,
	                ext_1 = cvt1_vt_vc.data[0]->op_std->ext_1;

	struct Matrix_T* cvt1r = constructor_empty_Matrix_T('R',n_eq*ext_0,n_vr*ext_1); // destructed

	struct Matrix_T* cvt1r_l = constructor_empty_Matrix_T('R',ext_0,ext_1); // destructed
	const struct const_Multiarray_T* dfr_ds_Ma = flux_r->dfr_ds;
	struct Vector_T dfr_ds = { .ext_0 = dfr_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_T(cvt1r_l,0.0);
		for (int dim = 0; dim < DIM; ++dim) {
			const ptrdiff_t ind =
				compute_index_sub_container(dfr_ds_Ma->order,1,dfr_ds_Ma->extents,(ptrdiff_t[]){eq,vr,dim});
			dfr_ds.data = (Type*)&dfr_ds_Ma->data[ind];
			mm_diag_T('R',1.0,1.0,cvt1_vt_vc.data[dim]->op_std,(struct const_Vector_T*)&dfr_ds,cvt1r_l,false);
		}
		set_block_Matrix_T(cvt1r,(struct const_Matrix_T*)cvt1r_l,eq*ext_0,vr*ext_1,'i');
	}}
	destructor_Matrix_T(cvt1r_l);

	struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) dpg_s_vol;
	const struct const_Vector_R* w_vc = get_operator__w_vc__s_e_T(s_vol);
	const struct const_Vector_R J_vc  = interpret_const_Multiarray_as_Vector_R(s_vol->jacobian_det_vc);

	const struct const_Vector_R* J_inv_vc = constructor_inverse_const_Vector_R(&J_vc);               // destructed
	const struct const_Vector_R* wJ_vc    = constructor_dot_mult_const_Vector_R(w_vc,J_inv_vc,n_vr); // destructed
	destructor_const_Vector_R(J_inv_vc);

	const struct const_Matrix_T* n1_l =
		constructor_mm_diag_const_Matrix_T_R(1.0,(struct const_Matrix_T*)cvt1r,wJ_vc,'R',false); // destructed
	destructor_const_Vector_R(wJ_vc);

	const struct const_Matrix_T* n1 =
		constructor_mm_const_Matrix_T('N','T',1.0,n1_l,(struct const_Matrix_T*)cvt1r,'R'); // destructed
	destructor_const_Matrix_T(n1_l);
	destructor_Matrix_T(cvt1r);

	const struct const_Matrix_R* norm_op_H0 = dpg_s_vol->norm_op_H0;
	assert(norm_op_H0->ext_0 == ext_0);

	struct Matrix_T* norm_op = constructor_empty_Matrix_T('R',n_eq*ext_0,n_eq*ext_0); // moved

	set_block_Matrix_T(norm_op,n1,0,0,'i');
	for (int eq = 0; eq < n_eq; ++eq)
		set_block_Matrix_T_R(norm_op,norm_op_H0,eq*ext_0,eq*ext_0,'a');
	destructor_const_Matrix_T(n1);

	struct Norm_DPG* norm = malloc(sizeof* norm); // returned
	norm->N     = (struct const_Matrix_T*) norm_op;
	norm->dN_ds = constructor_default_const_Matrix_T();

	return norm;
}

static void compute_rlhs_1
	(const struct S_Params_DPG* s_params, struct Flux_Input_T* flux_i, const struct DPG_Solver_Volume_T* dpg_s_vol,
	 struct Solver_Storage_Implicit* ssi, const struct Simulation*const sim, struct Simulation*const sim_c)
{
	const struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) dpg_s_vol;

	struct Flux_Ref_T* flux_r = constructor_Flux_Ref_vol_T(&s_params->spvs,flux_i,s_vol,sim); // destructed

	const struct Norm_DPG* norm = s_params->constructor_norm_DPG(dpg_s_vol,flux_r,sim); // destructed

	struct Vector_T* rhs_std = constructor_rhs_v_1(flux_r,s_vol,sim); // destructed
	struct Matrix_T* lhs_std = constructor_lhs_v_1_T(flux_r,s_vol,sim); // destructed

	add_to_rlhs__face_T(rhs_std,&lhs_std,dpg_s_vol,sim,true);
	increment_rhs_source(rhs_std,s_vol,sim);

	const struct const_Matrix_T* opt_t =
		constructor_sysv_const_Matrix_T(norm->N,(struct const_Matrix_T*)lhs_std); // destructed

	const struct const_Vector_T* rhs_opt =
		constructor_mv_const_Vector_T('T',-1.0,opt_t,(struct const_Vector_T*)rhs_std); // destructed

#if TYPE_RC == TYPE_REAL
	struct Matrix_T* lhs_opt = constructor_mm_Matrix_T('T','N',1.0,opt_t,(struct const_Matrix_T*)lhs_std,'R'); // dest.

	add_to_lhs_opt__d_opt_t_ds(flux_r,dpg_s_vol,norm,opt_t,(struct const_Vector_T*)rhs_std,
	                           (struct const_Matrix_T*)lhs_std,lhs_opt,sim,sim_c);

	add_to_petsc_Mat_Vec_dpg(s_vol,rhs_opt,(struct const_Matrix_T*)lhs_opt,ssi,sim);
	destructor_Matrix_T(lhs_opt);
#elif TYPE_RC == TYPE_COMPLEX
	UNUSED(sim_c);
	add_to_petsc_Mat_dpg_c(s_vol,rhs_opt,ssi);
#endif
	destructor_Flux_Ref_T(flux_r);
	destructor_Norm_DPG(norm);

	destructor_Vector_T(rhs_std);
	destructor_Matrix_T(lhs_std);
	destructor_const_Vector_T(rhs_opt);

	destructor_const_Matrix_T(opt_t);
}

static void scale_by_Jacobian (const struct Numerical_Flux_T* num_flux, const struct Solver_Face_T* s_face)
{
	const struct Face* face = (struct Face*) s_face;

	assert(face->boundary);
	assert(num_flux->neigh_info[0].dnnf_ds != NULL || num_flux->neigh_info[0].dnnf_dg != NULL);

	const struct const_Vector_R jacobian_det_fc = interpret_const_Multiarray_as_Vector_R(s_face->jacobian_det_fc);
	scale_Multiarray_T_by_Vector_R('L',1.0,(struct Multiarray_T*)num_flux->nnf,&jacobian_det_fc,false);

	if (num_flux->neigh_info[0].dnnf_ds)
		scale_Multiarray_T_by_Vector_R(
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

	set_block_Matrix_T(lhs,(struct const_Matrix_T*)lhs_ll,0,0,'a');
#if 0
printf("lhs\n");
print_Matrix_T(lhs_ll);
print_Matrix_T(lhs);
#endif
	destructor_Matrix_T(lhs_ll);
}

// Level 2 ********************************************************************************************************** //

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

static void destructor_Norm_DPG (const struct Norm_DPG* norm)
{
	destructor_const_Matrix_T(norm->N);
/// \todo Make conditional destructor.
	if (norm->dN_ds)
		destructor_const_Matrix_T(norm->dN_ds);
	free((void*)norm);
}
