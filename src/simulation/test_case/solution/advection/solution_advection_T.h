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
 *  \brief Provides templated functions relating to the linear Advection solutions.
 */

#include "def_templates_solution_advection.h"
#include "def_templates_test_case.h"

#include <stdbool.h>

struct Test_Case_T;
struct Simulation;

/** \brief Function pointer to linear Advection velocity vector computing functions.
 *  \return An array holding the values of the vector at the given input nodes.
 *
 *  \param xyz The xyz coordinates.
 */
typedef const Real* (*compute_b_adv_fptr_T)
	(const Type*const xyz
	);

/// \brief Container for solution data relating to linear advection solutions.
struct Sol_Data__Advection_T {
	compute_b_adv_fptr_T compute_b_adv; ///< \ref compute_b_adv_fptr.

	bool use_constant_solution; ///< Flag for whether a constant solution should be used.

	Real u_scale; ///< Scaling constant for the solution.
};

/// \brief Set the solution function pointer members of an Advection \ref Test_Case_T.
void set_function_pointers_solution_advection_T
	(struct Test_Case_T* test_case,    ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
		);

/** \brief Return the statically allocated \ref Sol_Data__Advection container.
 *  \return See brief. */
struct Sol_Data__Advection_T get_sol_data_advection_T
	( );

/// \brief Read the required solution data into the \ref Sol_Data__Advection container.
void read_data_advection_T
	(struct Sol_Data__Advection_T*const sol_data ///< \ref Sol_Data__Advection.
	);

/** \brief Version of \ref compute_b_adv_fptr for constant advection velocity throughout the domain.
 *  \return See brief. */
const Real* compute_b_adv_constant_T
	(const Type*const xyz ///< See brief.
		);

/** \brief Version of \ref compute_b_adv_fptr for constant magnitude advection velocity with angle varying over a
 *         cylinder.
 *  \return See brief. */
const Real* compute_b_adv_vortex_T
	(const Type*const xyz ///< See brief.
		);

#include "undef_templates_solution_advection.h"
#include "undef_templates_test_case.h"
