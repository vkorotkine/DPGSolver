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


#ifndef DPG__gradient_h__INCLUDED
#define DPG__gradient_h__INCLUDED

#include "adjoint.h"
#include "sensitivities.h"

struct Matrix_d;
struct Optimization_Case;

/** \brief Container for the gradient information. 
 */
struct Gradient_Data {
	struct Matrix_d *Gradient;  ///< The gradient of the specified functional with respect to the design variables
};


/** \brief Create the gradient_data data structure. Allocate the correct amount of memory for the 
 * Gradient_Data data structures (matrices), to be filled eventually.
 *
 * \return Gradient_Data data structure with the data structures ready to be filled.
 */
struct Gradient_Data* constructor_Gradient_Data (
	struct Adjoint_Data *adjoint_data, ///< The data structure holding the adjoint information
	struct Sensitivity_Data *sensitivity_data ///< The data structure holding the sensitivity information
	);


/** \brief Free the Gradient_Data data structure and clear allocated memory in the constructor
 */
void destructor_Gradient_Data (
	struct Gradient_Data *gradient_data ///< The sensitivity_data structure to be freed
	);


/** \brief Compute the gradient using the adjoint and the sensitivities. 
 *	That is, solve the equation
 *		
 *		Gradient = [dI/dXp]^T - Chi^T * [dR/dXp]
 *  
 *	where Chi^T is the transpose of the adjoint. Store the result in the Gradient_Data data structure (in Gradient).
 */
void compute_gradient(
	struct Gradient_Data *gradient_data, ///< The data structure to store the gradient in
	struct Adjoint_Data *adjoint_data, ///< The data structure with the adjoint information
	struct Sensitivity_Data *sensitivity_data ///< The data structure holding the sensitivity information
	);


void test_brute_force_gradient(struct Optimization_Case *optimization_case);

#endif // DPG__gradient_h__INCLUDED

