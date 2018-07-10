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

#ifndef DPG__adjoint_h__INCLUDED
#define DPG__adjoint_h__INCLUDED

#include "functionals.h"

struct Optimization_Case;
struct Matrix_d;
struct Simulation;

/** \brief Container for the adjoint information. 
 * NOTE: The LHS is not stored as a Matrix_d since it will be directly loaded 
 * into the LHS_petsc structures in solve_adjoint.
 */
struct Adjoint_Data {

	struct Matrix_d *RHS, ///< RHS = [dI/dW] (will be set and used to initialize the petsc vector for the RHS)
					*Chi; ///< Chi = Adjoint (will be set only after the adjoint is solved)

};


/** \brief Create the adjoint data structure. Allocate the correct amount of memory for the 
 * adjoint data structures, to be filled eventually.
 *
 * \return Adjoint_Data data structure with the data structures ready to be filled for the 
 * adjoint solve
 */
struct Adjoint_Data* constructor_Adjoint_Data (
	struct Optimization_Case *optimization_case ///< Standard. The optimization data structure
	);


/** \brief Free the Adjoint_Data data structure and clear allocated memory in the constructor
 */
void destructor_Adjoint_Data (
	struct Adjoint_Data* adjoint_data ///< The adjoint_data structure to be freed
	);


/** \brief 	Setup the adjoint vectors in preparation of the adjoint sovle. 
 * 	This function will compute the RHS terms and the initial Chi (adjoint). 
 *  
 *	RHS: Computed using the complex step. Load data into Matrix_d RHS structure in Adjoint_Data
 * 	Chi: Set to be initially equal to the RHS. Chi will then be solved for iteratively during 
 * 		the adjoint solve. Fill the data into the Matrix_d in Adjoint_Data.
 */
void setup_adjoint(
	struct Adjoint_Data* adjoint_data, ///< The adjoint_data data structure to load data into
	struct Simulation *sim, ///< The real simulation object (holding the real volumes and faces)
	struct Simulation *sim_c, ///< The complex simulation object. Used for the complex step.
	functional_fptr functional, ///< The functional function pointer (real version)
	functional_fptr_c functional_c ///< The functional function pointer (complex version for complex step)
	);


/** \brief Solve the adjoint equation using the PETSc KSP solver. The final solution will
 * 	be loaded into the Chi Matrix_d data structure in adjoint_data. 
 * 	
 *	First all data is loaded into PETSc data structures. The LHS is obtained from the 
 * 	analytical linearization. The KSP solve is then called and all PETSc data structures
 *	are then destroyed
 */
void solve_adjoint(
	struct Adjoint_Data* adjoint_data, ///< Standard. The adjoint_data data structure
	struct Simulation *sim ///< The real simulation object (holding the real volumes and faces)
	);

#endif // DPG__adjoint_h__INCLUDED


