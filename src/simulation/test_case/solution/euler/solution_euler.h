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

#ifndef DPG__solution_euler_h__INCLUDED
#define DPG__solution_euler_h__INCLUDED
/** \file
 *  \brief Provides real functions relating to the Euler solutions.
 */

#include "def_templates_type_d.h"
#include "solution_euler_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "solution_euler_T.h"
#include "undef_templates_type.h"

struct Multiarray_d;
struct const_Multiarray_d;

/// \brief Compute the entropy measure.
void compute_entropy
	(struct Multiarray_d* s,                ///< The container to hold the entropy data.
	 const struct const_Multiarray_d* vars, ///< The container of Euler variables.
	 const char var_type                    ///< The type of the variables.
	);

/// \brief Compute the mach number.
void compute_mach
	(struct Multiarray_d* mach,             ///< The container to hold the mach data.
	 const struct const_Multiarray_d* vars, ///< The container of Euler variables.
	 const char var_type                    ///< The type of the variables.
	);

/// \brief Compute the temperature.
void compute_temperature
	(struct Multiarray_d*const t,                ///< The container to hold the temperature data.
	 const struct const_Multiarray_d*const vars, ///< The container of Euler variables.
	 const char var_type                         ///< The type of the variables.
	);

/// \brief Compute the maximum wave speed.
void compute_max_wavespeed
	(struct Multiarray_d* V_p_c,            ///< The container to hold the maximum wave speed data.
	 const struct const_Multiarray_d* vars, ///< The container of Euler variables.
	 const char var_type                    ///< The type of the variables.
	);

/** \brief Compute the drag and lift coefficient **values** at the input nodes.
 *
 *  \note This function computes the values to be subsequently used to in the integral for the drag/lift coefficient
 *        computation.
 */
void compute_cd_cl_values
	(struct Multiarray_d* c_dl,                    ///< The container to hold the cd/cl data.
	 const struct const_Multiarray_d*const vars,   ///< The container of Euler variables.
	 const char var_type,                          ///< The type of the variables.
	 const struct const_Multiarray_d*const normals ///< The container of unit normal vectors.
	);

/** \brief Function to be used for constructing functionals of reference drag/lift coefficients with values specified in
 *         the test case solution file.
 *  \return See brief.
 *
 *  \note As the values returned from this function are generally subtracted from the computed __local__ drag/lift
 *        coefficient values, the value specified in the input file should correspond to the total divided by the
 *        surface area. If the surface area is not available, it is suggested to enter 0.0 for the reference values and
 *        compute the errors and convergence rates externally as the values returned by the code will be incorrect.
 */
const struct const_Multiarray_d* constructor_const_functionals_cd_cl_reference_constant
	(const struct const_Multiarray_d* xyz, ///< Defined for \ref constructor_sol_fptr_T.
	 const struct Simulation* sim          ///< Defined for \ref constructor_sol_fptr_T.
		);

#endif // DPG__solution_euler_h__INCLUDED
