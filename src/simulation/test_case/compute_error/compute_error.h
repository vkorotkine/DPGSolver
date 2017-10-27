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

#ifndef DPG__compute_error_h__INCLUDED
#define DPG__compute_error_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used for error computation and output.
 */

struct Error_CE;
struct Simulation;
struct Solver_Volume;
struct Vector_d;
struct const_Multiarray_d;

#include <stddef.h>
#include "definitions_alloc.h"

/** \brief Function pointer to error computing functions (which construct \ref Error_CE containers).
 *  \param sim \ref Simulation.
 */
typedef struct Error_CE* (*constructor_Error_CE_fptr)
	(const struct Simulation* sim
	);

/// Container holding information relating to the computational element errors.
struct Error_CE {
	const ptrdiff_t dof; ///< The number of degrees of freedom.
	const double domain_volume; ///< The volume of the domain.

	const struct const_Vector_d* sol_L2;         ///< The solution measured in L2.
	const struct const_Vector_i* expected_order; ///< The expected convergence order of each variable.

	const char* header_spec; ///< The header for the file specific to the variables for which the error was computed.
};

// Interface functions ********************************************************************************************** //

/** \brief Output the error of the solution.
 *
 *  The error, for a generic variable \f$v\f$, is computed as:
 *  \f[
 *  ||v||_{L^2(\Omega)} =
 *  \left(
 *  \frac{\int_{\Omega} (v-v_{\text{exact}})^2 d \Omega}{\int_{\Omega} d \Omega}
 *  \right)^{\frac{1}{2}}
 *  \f]
 */
void output_error
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Compute the volume of the domain.
 *  \return See brief. */
double compute_domain_volume
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Increment the global squared \f$L^2\f$ errors with the contribution from the current volume.
void increment_vol_errors_l2_2
	(struct Vector_d* errors_l2_2,           ///< Holds the global squared l2 errors.
	 const struct const_Multiarray_d* err_v, ///< Holds the error values for the current volume.
	 const struct Solver_Volume* s_vol       ///< Current \ref Solver_Volume.
	);

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom.
 *  \return See brief. */
ptrdiff_t compute_dof
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Return a statically allocated `char*` holding the name of the error file.
 *  \return See brief. */
const char* compute_error_file_name
	(const struct Simulation* sim ///< \ref Simulation.
	);

#endif // DPG__compute_error_h__INCLUDED
