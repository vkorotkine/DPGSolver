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

#include "compute_error.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "element_error.h"
#include "volume_solver.h"

#include "intrusive.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Destructor for a \ref Error_CE container.
static void destructor_Error_CE
	(struct Error_CE* error_ce ///< Standard.
	);

/** \brief Compute the volume of the input \ref Solver_Volume.
 *  \return See brief. */
static double compute_volume
	(const struct Solver_Volume* s_vol ///< \ref Solver_Volume.
	);

/** \brief Constructor for a \ref const_Vector_d\* holding the cubature weights multiplied by the Jacobian determinants.
 *  \return See brief. */
static const struct const_Vector_d* constructor_w_detJ
	(const struct Solver_Volume* s_vol ///< \ref Solver_Volume.
	);

// Interface functions ********************************************************************************************** //

void output_error (const struct Simulation* sim)
{
	assert(sim->volumes->name == IL_SOLVER_VOLUME);
	assert(sim->faces->name   == IL_SOLVER_FACE);

	constructor_derived_Elements((struct Simulation*)sim,IL_SOLUTION_ELEMENT);
	constructor_derived_Elements((struct Simulation*)sim,IL_ELEMENT_ERROR);

	struct Error_CE* error_ce = sim->test_case->constructor_Error_CE(sim);

// compute dof
// output dof, domain_volume, l2 errors with appropriate string header based on variables present.
// assert(mpi_rank == 0) and output specific error file for this case.
EXIT_UNSUPPORTED;

	destructor_derived_Elements((struct Simulation*)sim,IL_SOLUTION_ELEMENT);
	destructor_derived_Elements((struct Simulation*)sim,IL_ELEMENT);

	destructor_Error_CE(error_ce);
}

double compute_domain_volume (const struct Simulation* sim)
{
	double domain_volume = 0.0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		domain_volume += compute_volume((struct Solver_Volume*) curr);
	return domain_volume;
}

void increment_vol_errors_l2_2
	(struct Vector_d* errors_l2_2, const struct const_Multiarray_d* err_v, const struct Solver_Volume* s_vol)
{
	assert(errors_l2_2->ext_0 == err_v->extents[1]);

	const struct const_Vector_d* w_detJ = constructor_w_detJ(s_vol); // destructed
	const ptrdiff_t ext_0 = w_detJ->ext_0;

	const ptrdiff_t n_out = errors_l2_2->ext_0;
	for (int i = 0; i < n_out; ++i) {
		const double* err_data = get_col_const_Multiarray_d(i,err_v);
		for (int j = 0; j < ext_0; ++j)
			errors_l2_2->data[i] += w_detJ->data[j]*err_data[j]*err_data[j];
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void destructor_Error_CE (struct Error_CE* error_ce)
{
	destructor_const_Vector_d(error_ce->sol_L2);
}

static double compute_volume (const struct Solver_Volume* s_vol)
{
	const struct const_Vector_d* w_detJ = constructor_w_detJ(s_vol); // destructed
	const ptrdiff_t ext_0 = w_detJ->ext_0;

	double volume = 0.0;
	for (int i = 0; i < ext_0; ++i)
		volume += w_detJ->data[i];

	destructor_const_Vector_d(w_detJ);

	return volume;
}

static const struct const_Vector_d* constructor_w_detJ (const struct Solver_Volume* s_vol)
{
	struct Volume* b_vol = (struct Volume*)s_vol;
	struct Error_Element* e = (struct Error_Element*) b_vol->element;

	const int curved = b_vol->curved,
	          p      = s_vol->p_ref;
	const struct const_Vector_d* w_vc = get_const_Multiarray_Vector_d(e->w_vc[curved],(ptrdiff_t[]){0,0,p,p});
	const struct const_Vector_d jacobian_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
	assert(w_vc->ext_0 == jacobian_det_vc.ext_0);

	const ptrdiff_t ext_0 = w_vc->ext_0;

	struct Vector_d* w_detJ = constructor_empty_Vector_d(ext_0); // returned

	for (int i = 0; i < ext_0; ++i)
		w_detJ->data[i] = w_vc->data[i]*jacobian_det_vc.data[i];

	return (const struct const_Vector_d*) w_detJ;
}
