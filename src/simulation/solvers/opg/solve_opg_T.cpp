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

#include "solve_opg_T.h"

#include <assert.h>

#include "macros.h"


#include "def_templates_solve_opg.h"

#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

#include "def_templates_multiarray.h"

#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom for the global test function computation.
 *  \return See brief. */
static ptrdiff_t compute_dof_test
	(const struct Simulation*const sim ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct Vector_i* constructor_nnz_opg_T (const struct Simulation* sim)
{
	assert(get_set_has_1st_2nd_order(NULL)[1] == false); // Add support.

	const ptrdiff_t dof = compute_dof_test(sim);
	struct Vector_i* nnz = constructor_zero_Vector_i(dof); // returned

	// Volume contribution (Diagonal)
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;

		struct Multiarray_T* test_s_coef = s_vol->test_s_coef;
		const ptrdiff_t size = compute_size(test_s_coef->order,test_s_coef->extents);
		increment_nnz(nnz,s_vol->ind_dof_test,size,size);
	}

	// Face contributions (Off-diagonal)
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face* face = (struct Face*) curr;
		if (face->boundary)
			continue;

		struct Solver_Volume_T* s_vol[2] = { (struct Solver_Volume_T*) face->neigh_info[0].volume,
		                                     (struct Solver_Volume_T*) face->neigh_info[1].volume, };

		struct Multiarray_T* test_s_coef[2] = { s_vol[0]->test_s_coef, s_vol[1]->test_s_coef, };
		const ptrdiff_t size[2] = { compute_size(test_s_coef[0]->order,test_s_coef[0]->extents),
		                            compute_size(test_s_coef[1]->order,test_s_coef[1]->extents), };

		increment_nnz(nnz,s_vol[0]->ind_dof_test,size[0],size[1]);
		increment_nnz(nnz,s_vol[1]->ind_dof_test,size[1],size[0]);
	}

	// Constraint - if applicable (Diagonal and Off-diagonal)
	if (test_case_explicitly_enforces_conservation(sim)) {
		EXIT_ADD_SUPPORT;
	}

	return nnz;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom in the volume computational elements for test functions.
 *  \return See brief. */
static ptrdiff_t compute_dof_volumes_test
	(const struct Simulation*const sim ///< \ref Simulation.
	);

static ptrdiff_t compute_dof_test (const struct Simulation*const sim)
{
	ptrdiff_t dof = 0;
	dof += compute_dof_volumes_test(sim);
	return dof;
}

// Level 1 ********************************************************************************************************** //

static ptrdiff_t compute_dof_volumes_test (const struct Simulation*const sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;
		dof += compute_size(s_vol->test_s_coef->order,s_vol->test_s_coef->extents);
	}
	return dof;
}

#include "undef_templates_solve_opg.h"

#include "undef_templates_face_solver.h"
#include "undef_templates_volume_solver.h"

#include "undef_templates_multiarray.h"

#include "undef_templates_test_case.h"
