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

#include <stdbool.h>

/// \todo Potentially make a special header for function pointers to limit the dependency generation here.
#include "flux.h"
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

	const bool has_1st_order, ///< Flag for whether the pde under consideration has 1st order terms.
	           has_2nd_order; ///< Flag for whether the pde under consideration has 2nd order terms.

	const int n_var, ///< Number of variables in the PDE under consideration.
	          n_eq;  ///< Number of equations in the PDE under consideration.

	/// Function pointer to the function used to set \ref Solver_Volume::sol_coef.
	set_sol_coef_v_fptr set_sol_coef_v;

	/// Function pointer to the function used to set \ref Solver_Volume::grad_coef.
	set_sol_coef_v_fptr set_grad_coef_v;

	/// Function pointer to the function used to set \ref Solver_Face::sol_coef.
	set_sol_coef_f_fptr set_sol_coef_f;

	/// Function pointer to the function used to set \ref Solver_Face::grad_coef.
	set_sol_coef_f_fptr set_grad_coef_f;


	// Solver related parameters
	char solver_method_curr;     ///< The current solver method. Options: 'e'xplicit, 'i'mplicit.
	const bool display_progress; ///< Flag for whether the solver progress should be displayed (in stdout).

	/// The type of solver procedure to be used for the simulation. Options: See definitions_test_case.h.
	const int solver_proc;

	// Parameters for explicit simulations.
	const int solver_type_e; ///< The explicit solver type. Options: See definitions_test_case.h.

	double time;             ///< The current time.
	const double time_final; ///< The final time.
	const double dt;         ///< The time increment at each stage of the explicit solve.


	// Parameters for explicit to implicit simulations.
	const double exit_tol_e;   ///< The exit tolerance for the residual during the explicit solver stage.
	const double exit_ratio_e; ///< The exit ratio for the residual during the explicit solver stage.

	const bool flux_comp_mem_e[MAX_FLUX_OUT]; ///< \ref Flux_Input::compute_member for the explicit solver.
	const bool flux_comp_mem_i[MAX_FLUX_OUT]; ///< \ref Flux_Input::compute_member for the implicit solver.

	/// Function pointer to the function used to call the combination of 1st and 2nd order flux functions.
	compute_Flux_fptr compute_Flux;

/// \todo Merge the flux computation functions (without and with Jacobian) and remove one of these function pointers.
	/// Function pointers to the functions used to compute the 1st/2nd order fluxes for the explicit solver.
	compute_Flux_fptr compute_Flux_e[2];

	/** Function pointers to the functions used to compute the 1st/2nd order fluxes (and optionally flux Jacobians)
	 *  for the implicit solver. */
	compute_Flux_fptr compute_Flux_i[2];
};

/** \brief Constructor for a \ref Test_Case.
 *  \return See brief. */
struct Test_Case* constructor_Test_Case
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Test_Case.
void destructor_Test_Case
	(const struct Test_Case* test_case ///< \ref Test_Case.
	);

#endif // DPG__test_case_h__INCLUDED
