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

#include "solve_T.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "petscmat.h"
#include "petscvec.h"

#include "macros.h"
#include "definitions_test_case.h"
#include "definitions_intrusive.h"


#include "def_templates_solve.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

#include "def_templates_solve_dg.h"
#include "def_templates_solve_dpg.h"
#include "def_templates_solve_opg.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for a \ref Vector_T\* holding the 'n'umber of 'n'on-'z'ero entries in each row of the global
 *         system matrix.
 *  \return See brief. */
static struct Vector_i* constructor_nnz
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom in the volume computational elements.
 *  \return See brief. */
static ptrdiff_t compute_dof_volumes
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom in the face computational elements.
 *  \return See brief. */
static ptrdiff_t compute_dof_faces
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom in the 'L'agrange 'mult'iplier members of the volume
 *         computational elements.
 *  \return See brief. */
static ptrdiff_t compute_dof_volumes_l_mult
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom of the test space.
 *  \return See brief. */
static ptrdiff_t compute_dof_test_T
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Update \ref Solver_Volume_T::ind_dof_test for the method under consideration.
static void update_ind_dof_test_T
	(const struct Simulation*const sim ///< Standard.
	 );

// Interface functions ********************************************************************************************** //

struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit_T (const struct Simulation* sim)
{
	assert(sizeof(PetscInt) == sizeof(int)); // Ensure that all is working correctly if this is removed.

	update_ind_dof_T(sim);
	struct Vector_i* nnz = constructor_nnz(sim); // destructed
	const ptrdiff_t dof_solve = nnz->ext_0;

	struct Solver_Storage_Implicit* ssi = calloc(1,sizeof *ssi); // free
	switch (sim->method) {
		case METHOD_DG:  // fallthrough
		case METHOD_DPG: assert(dof_solve == compute_dof_T(sim));      break;
		case METHOD_OPG: // fallthrough
		case METHOD_OPGC0: assert(dof_solve == compute_dof_test_T(sim)); break;
		default:         EXIT_ERROR("Unsupported: %d.\n",sim->method); break;
	}

	MatCreateSeqAIJ(MPI_COMM_WORLD,(PetscInt)dof_solve,(PetscInt)dof_solve,0,nnz->data,&ssi->A); // destructed
	MatSetFromOptions(ssi->A);
	MatSetUp(ssi->A);

	VecCreateSeq(MPI_COMM_WORLD,(PetscInt)dof_solve,&ssi->b); // destructed
	VecSetFromOptions(ssi->b);
	VecSetUp(ssi->b);

	destructor_Vector_i(nnz);

	return ssi;
}

ptrdiff_t compute_dof_T (const struct Simulation* sim)
{
	assert((sim->method == METHOD_DG) || (sim->method == METHOD_DPG) || (sim->method == METHOD_OPG) ||
	       (sim->method == METHOD_OPGC0 || sim->method == METHOD_L2_PROJ));
	ptrdiff_t dof = 0;
	dof += compute_dof_volumes(sim);
	dof += compute_dof_faces(sim);
	dof += compute_dof_volumes_l_mult(sim);
	return dof;
}

void add_to_flux_imbalance_source_T
	(const struct const_Multiarray_T*const source_vc_w_J, const struct Solver_Volume_T*const s_vol,
	 const struct Simulation*const sim)
{
	UNUSED(sim);
	const struct const_Vector_R* w_vc = get_operator__w_vc__s_e_T(s_vol);

	const struct const_Matrix_T source_M = interpret_const_Multiarray_as_Matrix_T(source_vc_w_J);

	const struct const_Matrix_T* source_integral =
		constructor_mm_diag_const_Matrix_T_R(1.0,&source_M,w_vc,'L',false); // destructed

	const struct const_Vector_T* source_integral_sum =
		constructor_sum_const_Vector_T_const_Matrix_T('R',source_integral); // destructed
	destructor_const_Matrix_T(source_integral);

	const ptrdiff_t n_vr = s_vol->flux_imbalance->ext_0;
	for (int vr = 0; vr < n_vr; ++vr)
		s_vol->flux_imbalance->data[vr] += source_integral_sum->data[vr];
	destructor_const_Vector_T(source_integral_sum);
}

const struct Operator* get_operator__tw0_vt_vc_T (const struct Solver_Volume_T* s_vol)
{
	const struct Volume* vol       = (struct Volume*) s_vol;
	const struct Solver_Element* e = (struct Solver_Element*) vol->element;

	const int p = s_vol->p_ref,
	          curved = vol->curved;
	return get_Multiarray_Operator(e->tw0_vt_vc[curved],(ptrdiff_t[]){0,0,p,p});
}

void initialize_zero_memory_volumes_T (struct Intrusive_List* volumes)
{
	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) curr;

		struct Multiarray_T* ref_coef = NULL;
		switch (get_set_method(NULL)) {
		case METHOD_DG:  // fallthrough
		case METHOD_DPG:
			ref_coef = s_vol->sol_coef;
			break;
		case METHOD_OPG: // fallthrough
		case METHOD_OPGC0:
			ref_coef = s_vol->test_s_coef;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",get_set_method(NULL));
			break;
		}

		resize_Multiarray_T(s_vol->rhs,ref_coef->order,ref_coef->extents);
		set_to_value_Multiarray_T(s_vol->rhs,0.0);
	}
}

void update_ind_dof_T (const struct Simulation*const sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face_T* s_face = (struct Solver_Face_T*) curr;

		struct Multiarray_T* nf_coef = s_face->nf_coef;
		const ptrdiff_t size = compute_size(nf_coef->order,nf_coef->extents);
		if (size == 0) {
			const_cast_ptrdiff(&s_face->ind_dof,-1);
			continue;
		}

		const_cast_ptrdiff(&s_face->ind_dof,dof);
		dof += size;
	}

	if (test_case_explicitly_enforces_conservation(sim)) {
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;

			const_cast_ptrdiff(&s_vol->ind_dof_constraint,dof);

			struct Multiarray_T* l_mult = s_vol->l_mult;
			dof += compute_size(l_mult->order,l_mult->extents);
		}
	}

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;

		const_cast_ptrdiff(&s_vol->ind_dof,dof);

		struct Multiarray_T* sol_coef = s_vol->sol_coef;
		dof += compute_size(sol_coef->order,sol_coef->extents);
	}
	assert(dof == compute_dof(sim));

	update_ind_dof_test_T(sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom for the volume test functions.
 *  \return See brief. */
static ptrdiff_t compute_dof_test_volumes
	(const struct Simulation*const sim ///< \ref Simulation.
	);

static struct Vector_i* constructor_nnz (const struct Simulation*const sim)
{
	struct Vector_i* nnz = NULL;
	switch (sim->method) {
	case METHOD_DG:  nnz = constructor_nnz_dg_T(false,sim);  break;
	case METHOD_DPG: nnz = constructor_nnz_dpg_T(false,sim); break;
	case METHOD_OPG: // fallthrough
	case METHOD_OPGC0:
		nnz = constructor_nnz_opg_T(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",sim->method);
		break;
	}
	return nnz;
}

static ptrdiff_t compute_dof_volumes (const struct Simulation*const sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;
		dof += compute_size(s_vol->sol_coef->order,s_vol->sol_coef->extents);
	}
	return dof;
}

static ptrdiff_t compute_dof_faces (const struct Simulation*const sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face_T* s_face = (struct Solver_Face_T*) curr;
		dof += compute_size(s_face->nf_coef->order,s_face->nf_coef->extents);
	}
	return dof;
}

static ptrdiff_t compute_dof_volumes_l_mult (const struct Simulation*const sim)
{
	ptrdiff_t dof = 0;
	if (test_case_explicitly_enforces_conservation(sim)) {
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;
			dof += compute_size(s_vol->l_mult->order,s_vol->l_mult->extents);
		}
	}
	return dof;
}

static ptrdiff_t compute_dof_test_T (const struct Simulation*const sim)
{
	assert(sim->method == METHOD_OPG || sim->method == METHOD_OPGC0);
	ptrdiff_t dof = 0;
	dof += compute_dof_test_volumes(sim);
	return dof;
}

static void update_ind_dof_test_T (const struct Simulation*const sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;

		const_cast_ptrdiff(&s_vol->ind_dof_test,dof);

		struct Multiarray_T* test_s_coef = s_vol->test_s_coef;
		dof += compute_size(test_s_coef->order,test_s_coef->extents);
	}
}

// Level 1 ********************************************************************************************************** //

static ptrdiff_t compute_dof_test_volumes (const struct Simulation*const sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;
		dof += compute_size(s_vol->test_s_coef->order,s_vol->test_s_coef->extents);
	}
	return dof;
}

#include "undef_templates_solve.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_face_solver.h"
#include "undef_templates_volume_solver.h"

#include "undef_templates_solve_dg.h"
#include "undef_templates_solve_dpg.h"
#include "undef_templates_solve_opg.h"
#include "undef_templates_test_case.h"
