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

#ifndef DPG__solution_advection_h__INCLUDED
#define DPG__solution_advection_h__INCLUDED
/** \file
 *  \brief Provides real functions relating to the linear Advection solutions.
 */

#include "def_templates_type_d.h"
#include "def_templates_solution_advection.h"
#include "def_templates_test_case.h"
#include "solution_advection_T.h"
#include "undef_templates_type.h"
#include "undef_templates_solution_advection.h"
#include "undef_templates_test_case.h"


#include "definitions_core.h"
struct Simulation;

/** \brief Function pointer to linear Advection velocity vector computing functions.
 *  \return An array holding the values of the vector at the given input nodes.
 *
 *  \param xyz The xyz coordinates.
 */
typedef const double* (*compute_b_adv_fptr)
	(const double*const xyz
	);

/// \brief Container for solution data relating to linear advection solutions.
struct Sol_Data__Advection {
	compute_b_adv_fptr compute_b_adv; ///< \ref compute_b_adv_fptr.

	double u_scale; ///< Scaling constant for the solution.
	double u_coef_polynomial4[5]; ///< Coefficients for the polynomial solution of maximal degree 4.

	/// Coefficients for the odd polynomial advection velocity vector.
	const struct const_Vector_d* b_coef_polynomial_odd;
};

/** \brief Return the statically allocated \ref Sol_Data__Advection container.
 *  \return See brief. */
struct Sol_Data__Advection get_sol_data_advection
	( );

/// \brief Read the required solution data into the \ref Sol_Data__Advection container.
void read_data_advection
	(struct Sol_Data__Advection*const sol_data ///< \ref Sol_Data__Advection.
	);

/** \brief Version of \ref compute_b_adv_fptr for constant advection velocity throughout the domain.
 *  \return See brief. */
const double* compute_b_adv_constant
	(const double*const xyz ///< See brief.
	);

/** \brief Version of \ref compute_b_adv_fptr for constant magnitude advection velocity with angle varying over a
 *         cylinder.
 *  \return See brief. */
const double* compute_b_adv_vortex
	(const double*const xyz ///< See brief.
	);

/** \brief Version of \ref compute_b_adv_fptr for polynomial advection velocity with behaviour similar to that of
 *         \ref compute_b_adv_vortex.
 *  \return See brief.
 *
 *  \note In the special case of only the first coefficient being non-zero, the advection velocity from
 */
const double* compute_b_adv_vortex_poly
	(const double*const xyz ///< See brief.
	);

#endif // DPG__solution_advection_h__INCLUDED
