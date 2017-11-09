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

#include "volume_solver_dg.h"

#include <string.h>

#include "macros.h"
#include "definitions_test_case.h"

#include "element_solver_dg.h"
#include "volume.h"
#include "volume_solver.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// Container holding flags for which members of \ref DG_Solver_Volume are needed.
struct Needed_Members {
	bool sol_coef_p, ///< Flag for \ref DG_Solver_Volume::sol_coef_p.
	     m_inv;      ///< Flag for \ref DG_Solver_Volume::m_inv.
};

/** \brief Return a statically allocated \ref Needed_Members container with values set.
 *  \return See brief. */
static struct Needed_Members set_needed_members
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for the inverse mass matrix of the input volume.
 *  \return See brief. */
static const struct const_Matrix_d* constructor_inverse_mass
	(const struct DG_Solver_Volume* volume ///< \ref DG_Solver_Volume.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_DG_Solver_Volume (struct Volume* volume_ptr, const struct Simulation* sim)
{
	struct Needed_Members needed_members = set_needed_members(sim);

	struct Solver_Volume* s_volume  = (struct Solver_Volume*) volume_ptr;
	struct DG_Solver_Volume* volume = (struct DG_Solver_Volume*) volume_ptr;

	const int order = s_volume->sol_coef->order;
	ptrdiff_t* extents = s_volume->sol_coef->extents;

	volume->rhs        = constructor_empty_Multiarray_d('C',order,extents); // destructed
	volume->sol_coef_p =
		( needed_members.sol_coef_p ? constructor_empty_Multiarray_d('C',order,extents) : NULL ); // destructed

	volume->m_inv = ( needed_members.m_inv ? constructor_inverse_mass(volume) : NULL ); // destructed
}

void destructor_derived_DG_Solver_Volume (struct Volume* volume_ptr)
{
	struct DG_Solver_Volume* volume = (struct DG_Solver_Volume*) volume_ptr;

	destructor_Multiarray_d(volume->rhs);
	if (volume->sol_coef_p)
		destructor_Multiarray_d(volume->sol_coef_p);
	if (volume->m_inv)
		destructor_const_Matrix_d(volume->m_inv);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Needed_Members set_needed_members (const struct Simulation* sim)
{
	const struct Test_Case* test_case = sim->test_case;
	struct Needed_Members needed_members =
		{ .sol_coef_p = false,
		  .m_inv      = false, };

	switch (test_case->solver_proc) {
	case SOLVER_E: // fallthrough
	case SOLVER_EI:
		if (!sim->collocated)
			needed_members.m_inv = true;
		switch (test_case->solver_type_e) {
		case SOLVER_E_SSP_RK_33: // fallthrough
		case SOLVER_E_LS_RK_54:
			needed_members.sol_coef_p = true;
			break;
		case SOLVER_E_EULER:
			// Do nothing
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",test_case->solver_type_e);
			break;
		}
		break;
	case SOLVER_I:
		// Do nothing
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->solver_proc);
		break;
	}

	return needed_members;
}

static const struct const_Matrix_d* constructor_inverse_mass (const struct DG_Solver_Volume* volume)
{
	struct Volume* vol          = (struct Volume*) volume;
	struct Solver_Volume* s_vol = (struct Solver_Volume*) volume;

	struct Solver_Element* s_e       = (struct Solver_Element*) vol->element;
	struct DG_Solver_Element* dg_s_e = (struct DG_Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;
	const struct Operator* cv0_vs_vc = get_Multiarray_Operator(s_e->cv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});
	const struct const_Vector_d* w_vc = get_const_Multiarray_Vector_d(dg_s_e->w_vc[curved],(ptrdiff_t[]){0,0,p,p});

	const struct const_Vector_d jacobian_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
	const struct const_Vector_d* wJ_vc = constructor_dot_mult_const_Vector_d(w_vc,&jacobian_det_vc); // destructed

	const struct const_Matrix_d* m_l = cv0_vs_vc->op_std;
	const struct const_Matrix_d* m_r = constructor_mm_diag_const_Matrix_d(1.0,m_l,wJ_vc,'L',false); // destructed
	destructor_const_Vector_d(wJ_vc);

	const struct const_Matrix_d* m = constructor_mm_const_Matrix_d('T','N',1.0,m_l,m_r,'R'); // destructed
	destructor_const_Matrix_d(m_r);

	const struct const_Matrix_d* m_inv = constructor_inverse_const_Matrix_d(m); // returned
	destructor_const_Matrix_d(m);

	return m_inv;
}
