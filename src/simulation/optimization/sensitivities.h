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

#ifndef DPG__sensitivities_h__INCLUDED
#define DPG__sensitivities_h__INCLUDED

#include "functionals.h"

struct Optimization_Case;
struct Matrix_d;
struct Simulation;


/** \brief Container for the sensitivity information. 
 */
struct Sensitivity_Data {

	struct Matrix_d 	*dR_dXp,  ///< Sensitivity of the residual (R) with respect to the design points (Xp)
						*dI_dXp;  ///< Sensitivity of the functional (I) with respect to the design points (Xp)

};


/** \brief Create the sensitivity data structure. Allocate the correct amount of memory for the 
 * sensivitity data structures (matrices), to be filled eventually.
 *
 * \return Sensitivity_Data data structure with the data structures ready to be filled.
 */
struct Sensitivity_Data* constructor_Sensitivity_Data (
	struct Optimization_Case *optimization_case ///< Standard. The optimization data structure
	);


/** \brief Free the Sensitivity_Data data structure and clear allocated memory in the constructor
 */
void destructor_Sensitivity_Data (
	struct Sensitivity_Data* sensitivity_data ///< The sensitivity_data structure to be freed
	);


/** \brief 	Compute the sensitivities of the Residual [dR/dXp] and the objective
 *	function [dI/dXp] with respect to the design control points Xp. Store the
 *	results in the Matrix_d structures in the Sensitivity_Data data structure.
 *
 *  NOTE: All sensitivities computed using the complex step.
 */
void compute_sensitivities(
	struct Sensitivity_Data *sensitivity_data, ///< The sensitivity_data data structure to load data into
	struct Optimization_Case *optimization_case, ///< Standard. The optimization data structure
	functional_fptr functional, ///< The functional function pointer (real version)
	functional_fptr_c functional_c ///< The functional function pointer (complex version for complex step)
	);


#endif // DPG__sensitivities_h__INCLUDED
