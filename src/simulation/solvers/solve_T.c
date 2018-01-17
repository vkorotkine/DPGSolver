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
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for a \ref Vector_T\* holding the 'n'umber of 'n'on-'z'ero entries in each row of the global
 *         system matrix.
 *  \return See brief. */
static struct Vector_i* constructor_nnz
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom in the volume computational elements.
 *  \return See brief. */
static ptrdiff_t compute_dof_volumes
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom in the face computational elements.
 *  \return See brief. */
static ptrdiff_t compute_dof_faces
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit_T (const struct Simulation* sim)
{
	assert(sizeof(PetscInt) == sizeof(int)); // Ensure that all is working correctly if this is removed.

	const ptrdiff_t dof = compute_dof_T(sim);
	update_ind_dof_T(sim);
	struct Vector_i* nnz = constructor_nnz(sim); // destructed

	struct Solver_Storage_Implicit* ssi = calloc(1,sizeof *ssi); // free

	MatCreateSeqAIJ(MPI_COMM_WORLD,(PetscInt)dof,(PetscInt)dof,0,nnz->data,&ssi->A); // destructed
	VecCreateSeq(MPI_COMM_WORLD,(PetscInt)dof,&ssi->b);                              // destructed

	destructor_Vector_i(nnz);

	return ssi;
}

ptrdiff_t compute_dof_T (const struct Simulation* sim)
{
	assert((sim->method == METHOD_DG) || (sim->method == METHOD_DPG)); // Ensure that all is working correctly if modified.
	ptrdiff_t dof = 0;
	dof += compute_dof_volumes(sim);
	dof += compute_dof_faces(sim);
	return dof;
}

void update_ind_dof_T (const struct Simulation* sim)
{
	switch (sim->method) {
	case METHOD_DG:  update_ind_dof_dg_T(sim);  break;
	case METHOD_DPG: update_ind_dof_dpg_T(sim); break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",sim->method);
		break;
	}
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
		constructor_sum_const_Vector_T_const_Matrix_T('C',source_integral); // destructed
	destructor_const_Matrix_T(source_integral);

	const ptrdiff_t n_vr = s_vol->flux_imbalance->ext_0;
	for (int vr = 0; vr < n_vr; ++vr)
		s_vol->flux_imbalance->data[vr] += source_integral_sum->data[vr];
	destructor_const_Vector_T(source_integral_sum);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Vector_i* constructor_nnz (const struct Simulation* sim)
{
	struct Vector_i* nnz = NULL;
	switch (sim->method) {
	case METHOD_DG:  nnz = constructor_nnz_dg_T(sim);  break;
	case METHOD_DPG: nnz = constructor_nnz_dpg_T(sim); break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",sim->method);
		break;
	}
	return nnz;
}

static ptrdiff_t compute_dof_volumes (const struct Simulation* sim)
{
	const bool enforcing_conservation = test_case_explicitly_enforces_conservation(sim);

	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;
		dof += compute_size(s_vol->sol_coef->order,s_vol->sol_coef->extents);
		if (enforcing_conservation)
			dof += compute_size(s_vol->l_mult->order,s_vol->l_mult->extents);
	}
	return dof;
}

static ptrdiff_t compute_dof_faces (const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face_T* s_face = (struct Solver_Face_T*) curr;
		dof += compute_size(s_face->nf_coef->order,s_face->nf_coef->extents);
	}
	return dof;
}
