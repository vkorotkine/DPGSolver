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
 *  \todo Attempt to template these functions.
 */

#include "test_complex_compute_face_rhs_dg.h"

#include "test_complex_boundary.h"
#include "test_complex_numerical_flux.h"
#include "test_complex_operators.h"
#include "test_complex_solve_dg.h"
#include "test_complex_test_case.h"

#include <assert.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "face_solver_dg_complex.h"
#include "volume_solver_dg_complex.h"

#include "complex_multiarray.h"
#include "multiarray.h"
#include "vector.h"

#include "compute_face_rlhs_dg.h"
#include "intrusive.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief `complex` version of \ref constructor_Numerical_Flux_Input_data.
void constructor_Numerical_Flux_Input_c_data
	(struct Numerical_Flux_Input_c* num_flux_i, ///< See brief.
	 const struct Face* face,                   ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/// \brief `complex` version of \ref destructor_Numerical_Flux_Input_data.
void destructor_Numerical_Flux_Input_c_data
	(struct Numerical_Flux_Input_c* num_flux_i ///< See brief.
	);

/// \brief `complex` version of \ref scale_by_Jacobian_e.
static void scale_by_Jacobian_e_c
	(const struct Numerical_Flux_c* num_flux, ///< See brief.
	 struct Face* face,                       ///< See brief.
	 const struct Simulation* sim             ///< See brief.
	);

/// \brief `complex` version of \ref compute_rhs_f_dg.
static void compute_rhs_f_dg_c
	(const struct Numerical_Flux_c* num_flux, ///< See brief.
	 struct Face* face,                       ///< See brief.
	 const struct Simulation* sim             ///< See brief.
	);

// Interface functions ********************************************************************************************** //

void compute_face_rhs_dg_c (const struct Simulation* sim, struct Intrusive_List* faces)
{
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);
	assert(sim->faces->name    == IL_FACE_SOLVER_DG_COMPLEX);
	assert(sim->volumes->name  == IL_VOLUME_SOLVER_DG_COMPLEX);

	constructor_derived_Complex_Test_Case((struct Simulation*)sim); // destructed
	struct Numerical_Flux_Input_c* num_flux_i = constructor_Numerical_Flux_Input_c(sim); // destructed

	for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next) {
		struct Face* face = (struct Face*) curr;

		constructor_Numerical_Flux_Input_c_data(num_flux_i,face,sim); // destructed

		struct Numerical_Flux_c* num_flux = constructor_Numerical_Flux_c(num_flux_i); // destructed
		destructor_Numerical_Flux_Input_c_data(num_flux_i);

		scale_by_Jacobian_e_c(num_flux,face,sim);

		compute_rhs_f_dg_c(num_flux,face,sim);
		destructor_Numerical_Flux_c(num_flux);
	}
	destructor_Numerical_Flux_Input_c(num_flux_i);
	destructor_derived_Complex_Test_Case((struct Simulation*)sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief `complex` version of \ref finalize_face_rhs_dg.
static void finalize_face_rhs_dg_c
	(const int side_index,                    ///< See brief.
	 const struct Numerical_Flux_c* num_flux, ///< See brief.
	 struct Face* face,                       ///< See brief.
	 const struct Simulation* sim             ///< See brief.
	);

void constructor_Numerical_Flux_Input_c_data
	(struct Numerical_Flux_Input_c* num_flux_i, const struct Face* face, const struct Simulation* sim)
{
	struct Complex_Test_Case* test_case = (struct Complex_Test_Case*) sim->test_case;
	struct Solver_Face* s_face                 = (struct Solver_Face*) face;
	struct Complex_DG_Solver_Face* c_dg_s_face = (struct Complex_DG_Solver_Face*) face;

	test_case->constructor_Boundary_Value_Input_c_face_fcl(&num_flux_i->bv_l,s_face,sim);          // destructed
	c_dg_s_face->constructor_Boundary_Value_c_fcl(&num_flux_i->bv_r,&num_flux_i->bv_l,s_face,sim); // destructed
}

void destructor_Numerical_Flux_Input_c_data (struct Numerical_Flux_Input_c* num_flux_i)
{
	destructor_Boundary_Value_Input_c(&num_flux_i->bv_l);
	destructor_Boundary_Value_c(&num_flux_i->bv_r);
}

static void scale_by_Jacobian_e_c
	(const struct Numerical_Flux_c* num_flux, struct Face* face, const struct Simulation* sim)
{
	UNUSED(sim);
	struct Solver_Face* s_face = (struct Solver_Face*)face;

	const struct const_Vector_d jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_c_by_Vector_d('L',1.0,(struct Multiarray_c*)num_flux->nnf,&jacobian_det_fc,false);
}

static void compute_rhs_f_dg_c
	(const struct Numerical_Flux_c* num_flux, struct Face* face, const struct Simulation* sim)
{
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);

	finalize_face_rhs_dg_c(0,num_flux,face,sim);
	if (!face->boundary) {
		permute_Multiarray_c_fc((struct Multiarray_c*)num_flux->nnf,'R',1,(struct Solver_Face*)face);

		scale_Multiarray_c((struct Multiarray_c*)num_flux->nnf,-1.0);
		finalize_face_rhs_dg_c(1,num_flux,face,sim);
	}
}

// Level 1 ********************************************************************************************************** //

static void finalize_face_rhs_dg_c
	(const int side_index, const struct Numerical_Flux_c* num_flux, struct Face* face, const struct Simulation* sim)
{
	UNUSED(sim);
	const struct Operator* tw0_vs_fc = get_operator__tw0_vs_fc__rlhs_dg(side_index,face);

	struct Complex_DG_Solver_Volume* c_dg_s_vol =
		(struct Complex_DG_Solver_Volume*) face->neigh_info[side_index].volume;

	mm_NNC_Operator_Multiarray_c(-1.0,1.0,tw0_vs_fc,num_flux->nnf,c_dg_s_vol->rhs,'d',2,NULL,NULL);
}
