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
 *  \brief Provides templated functions relating to the Navier-Stokes solutions.
 */

#include "def_templates_solution_navier_stokes.h"
#include "def_templates_multiarray.h"
#include "def_templates_test_case.h"

#include <stdbool.h>

struct Test_Case_T;
struct Simulation;
struct Multiarray_T;
struct const_Multiarray_T;

/** \brief Pointer to functions computing the value of the viscosity.
 *
 *  \param rho    The density.
 *  \param rhouvw The xyz momentum components.
 *  \param E      The total energy.
 */
typedef Type (*compute_mu_fptr_T)
	(const Type rho,
	 const Type*const rhouvw,
	 const Type E
	);

/// \brief Set the solution function pointer members of an Navier-Stokes \ref Test_Case_T.
void set_function_pointers_solution_navier_stokes_T
	(struct Test_Case_T*const test_case, ///< \ref Test_Case_T.
	 const struct Simulation*const sim   ///< \ref Simulation.
	);

/// \brief Convert gradients between supported variable types (see \ref convert_variables_T for options).
void convert_variables_gradients_T
	(struct Multiarray_T*const grad,            ///< The container holding the gradient data.
	 const struct const_Multiarray_T*const sol, ///< The container holding the solution data.
	 const char type_i,                         ///< The input variable type.
	 const char type_o                          ///< The output variable type.
	);

/** \brief Return the pointer to the appropriate \ref compute_mu_fptr_T specialization based on the viscosity type.
 *  \return See brief. */
compute_mu_fptr_T get_compute_mu_fptr_T
	( );

/// \brief Set the "viscosity_type" parameter based on the value in the input file.
void set_viscosity_type_T
	(int*const viscosity_type_ptr, ///< Pointer to the variable.
	 bool*const need_input         ///< Pointer to the flag for whether the input is still needed.
	);

/** \brief Version of \ref compute_mu_fptr_T for constant viscosity.
 *  \return See brief. */
Type compute_mu_constant_T
	(const Type rho,          ///< See brief.
	 const Type*const rhouvw, ///< See brief.
	 const Type E             ///< See brief.
	);

/** \brief Version of \ref compute_mu_fptr_T using the Sutherland formula (eq. (1.56), \cite Toro2009).
 *  \return See brief. */
Type compute_mu_sutherland_T
	(const Type rho,          ///< See brief.
	 const Type*const rhouvw, ///< See brief.
	 const Type E             ///< See brief.
	);

#include "undef_templates_solution_navier_stokes.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_test_case.h"
