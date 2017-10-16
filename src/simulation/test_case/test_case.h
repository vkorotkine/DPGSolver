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

#ifndef DPG__test_case_h__INCLUDED
#define DPG__test_case_h__INCLUDED
/** \file
 *  \brief Provides container(s) and functions relating to the test cases.
 */

#include "solution.h"

struct Simulation;

/** \brief Container for test case specific information.
 *
 *  This container is used to hold test case specific variables (such as those related to the pde under consideration)
 *  as well as function pointers to various functions such that the control flow is not broken during run-time.
 *
 *  As this container is a general template for many supported test cases, it is frequently the case that some of the
 *  function pointers are not necessary. In this case, they are set to point to functions which simply return
 *  immediately.
 */
struct Test_Case {
	const int pde_index; ///< Index corresponding to \ref Simulation::pde_name.

	const int n_var, ///< Number of variables in the PDE under consideration.
	          n_eq;  ///< Number of equations in the PDE under consideration.

	/// Function pointer to the function computing the initial \ref Solver_Volume::sol_coef.
	compute_sol_coef_v_fptr compute_init_sol_coef_v;

	/// Function pointer to the function computing the initial \ref Solver_Volume::grad_coef.
	compute_sol_coef_v_fptr compute_init_grad_coef_v;

	/// Function pointer to the function computing the initial \ref Solver_Face::sol_coef.
	compute_sol_coef_f_fptr compute_init_sol_coef_f;

	/// Function pointer to the function computing the initial \ref Solver_Face::grad_coef.
	compute_sol_coef_f_fptr compute_init_grad_coef_f;


	// Parameters for unsteady simulations.
	double time,       ///< Parameter storing the current time of the test case.
	       time_final; ///< Parameter storing the final time of the test case.
};

/** \brief Constructor for a \ref Test_Case.
 *  \return See brief. */
const struct Test_Case* constructor_Test_Case
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Test_Case.
void destructor_Test_Case
	(const struct Test_Case* test_case ///< \ref Test_Case.
	);

#endif // DPG__test_case_h__INCLUDED
