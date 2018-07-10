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

#ifndef DPG__optimizer_line_search_h__INCLUDED
#define DPG__optimizer_line_search_h__INCLUDED


struct Optimization_Case;


/** \brief Container for the line search optimization information
 */
struct Optimizer_Line_Search_Data {

	struct Matrix_d *Laplacian_Matrix; ///< Optional smoothing matrix

	// ===========================
	//     BFGS Data Structures
	// ===========================

	struct Matrix_d *B_k_inv, ///< BFGS matrix for Hessian approximation at k_th design iteration
					*s_k,  ///< BFGS matrix for step vector at k_th design iteration
					*grad_f_k; ///< BFGS matrix for the gradient of the objective at the k_th design iteration


	// ===========================
	//    Optimization Parameters
	// ===========================	

	double 	exit_L2_norm_gradient, ///< Exit condition based on the minimum L2 norm of the gradient of the objective
			exit_objective_value; ///< Exit condition based on the minimum value of the objective function value

	int exit_max_design_iterations; ///< Exit condition based on the maximum number of design iterations

	double 	step_size_init, ///< The initial step length for the search direction
			wolfe_condition_c, ///< The wolfe condition c value
			wolfe_condition_rho, ///< The wolfe condition rho value (how much to scale the search direction during the backtracking)
			wolfe_condition_alpha_min, ///< The minimum step length to use
			max_norm_search_vector; ///< The maximum norm of the search direction (if greater, will normalize by this amount)
};


/**	\brief Use a line search method (gradient descent or BFGS) to find the optimal shape
 */
void optimizer_line_search(
	struct Optimization_Case* optimization_case ///< Standard. The data structure holding all optimization data.
	);


#endif // DPG__optimizer_line_search_h__INCLUDED
