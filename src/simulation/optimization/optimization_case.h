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

#ifndef DPG__optimization_case_h__INCLUDED
#define DPG__optimization_case_h__INCLUDED

#include "objective_functions.h"

#include "petscvec.h"
#include "petscmat.h"

struct Multiarray_d;
struct Matrix_d;

/** \brief Container for test case specific information.
 *
 *  This container is used to hold optimization case specific variables as well as function 
 *  pointers to various functions such that the control flow is not broken during run-time.
 *  
 */

struct Optimization_Case {

	int num_design_pts_dofs; // The number of design control point degrees of freedom

	objective_function_fptr 	objective_function;  ///< Function pointer to the real objective function
	objective_function_fptr_c 	objective_function_c;  ///< Function pointer to the complex objective function


	// Pointers to the real (sim) and complex (sim_c) objects. 
	// A complex counterpart is needed for the complex step.
	struct Simulation *sim;
	struct Simulation *sim_c;


	// =================================
	//      Geometry Data Structures
	// =================================

	struct Geo_Data {

		// NURBS Parametric Domain orders in the xi direction (P) and 
		// eta direction (Q)
		int P, Q;

		// Knot vectors in the xi and eta direction
		struct Multiarray_d *knots_xi, *knots_eta; 

		// Control points and weights. A complex version is stored as well
		// in order to be used with the complex step
		struct Multiarray_d *control_points_and_weights;
		struct Multiarray_c *control_points_and_weights_c;

		// Control point connectivity information for the mesh
		struct Multiarray_i *control_point_connectivity;

		// The multiarray used to store which design points are used for the optimization.
		// Each row holds data of the form (control_pt_index, x dof bool, y dof bool). The
		// bools are 0 if the point cannot be moved in the given coordinate direction and 
		// 1 if it can.
		struct Multiarray_i *control_points_optimization;

	} geo_data;


	// =================================
	//     Adjoint Data Structures
	// =================================

	// 		RHS = [dI/dW] (will be set and used to initialize the petsc vector for the RHS)
	// 		Chi = Adjoint (will be set only after the adjoint is solved)
	// - Lifetime of structures:
	//		- All PETSc structures will be constructed in setup_adjoint and destructed in
	//			solve_adjoint
	//		- RHS will be constructed in setup_adjoint and destructed in solve_adjoint
	//		- Chi will be constructed in constructor_Optimization_Case and will be destructed
	//			in destructor_Optimization_Case (it is used in the whole optimization sequence)

	struct Multiarray_d *RHS, *Chi;

	Mat LHS_petsc; // [dR/dW]^T
	Vec RHS_petsc; // [dI/dW]
	Vec Chi_petsc; // Chi (Adjoint)

	// =================================
	//    Sensitivity Data Structures
	// =================================
	
	struct Multiarray_d  *dR_dXp,  // Sensitivity of the residual 
						 *dI_dXp;  // Sensitivity of the objective function


	// Gradient of the objective with respect to the design variables
	struct Multiarray_d *grad_I;


	// =================================
	// Optimization Minimizer Structures
	// =================================

	struct Matrix_d *Laplacian_Matrix;

	// BFGS Data structures
	struct Matrix_d *B_k_inv, *s_k, *grad_f_k;


};


/** \brief Container for test case specific information.
 *
 *  This container is used to hold optimization case specific variables as well as function 
 *  pointers to various functions such that the control flow is not broken during run-time.
 *  
 */

struct Optimization_Case* constructor_Optimization_Case (const struct Simulation* sim);

void destructor_Optimization_Case (struct Optimization_Case* optimization_case);

#endif // DPG__optimization_case_h__INCLUDED
