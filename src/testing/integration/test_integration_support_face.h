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

#ifndef DPG__test_integration_support_face_h__INCLUDED
#define DPG__test_integration_support_face_h__INCLUDED
/** \file
 *  \brief Provides support functions/structures for integration testing relating to face computational elements.
 *
 *  Functions included here are those which are, for example, used in multiple integration tests with face dependent
 *  data such as boundary and numerical flux integration tests.
 */

struct Numerical_Flux_Input;
struct Numerical_Flux_Input_c;
struct Solver_Face;
struct Solver_Face_c;
struct Simulation;

/** \brief Version of \ref constructor_Numerical_Flux_Input_data_T, forcing the gradients to be computed using
 *         interpolation if required.
 *
 *  Note that the gradient is not computed in this manner for the DG scheme as the partial gradient is used in the
 *  scheme.
 */
void constructor_Numerical_Flux_Input_data_with_gradients
	(struct Numerical_Flux_Input*const num_flux_i, ///< Standard.
	 const struct Solver_Face*const s_face,        ///< Standard.
	 const struct Simulation*const sim             ///< Standard.
		);

/// \brief Constructor for the data members of the complex \ref Numerical_Flux_Input_T container.
void constructor_Numerical_Flux_Input_c_data_members
	(struct Numerical_Flux_Input_c*const num_flux_c_i, ///< Complex \ref Numerical_Flux_Input_T.
	 struct Numerical_Flux_Input*const num_flux_i,     ///< Real    \ref Numerical_Flux_Input_T.
	 const char side                                   /**< Flag indicating on which side of face data members
	                                                    *   should be constructed. Options: 'l'eft, 'b'oth. */
		);

/// \brief Destructor for the 'l'eft data members of the complex \ref Numerical_Flux_Input_T container.
void destructor_Numerical_Flux_Input_c_data_members
	(struct Numerical_Flux_Input_c*const num_flux_c_i, ///< Complex \ref Numerical_Flux_Input_T.
	 const char side                                   /**< Flag indicating on which side of face data members
	                                                    *   should be constructed. Options: 'l'eft, 'b'oth. */
		);

/// \brief Version of \ref constructor_Numerical_Flux_Input_data_T constructing only the 'r'ight data members.
void constructor_Boundary_Value_c_data
	(struct Numerical_Flux_Input_c* num_flux_i, ///< See brief.
	 const struct Solver_Face_c* s_face,        ///< See brief.
	 const struct Simulation* sim               ///< See brief.
		);

/// \brief Version of \ref destructor_Numerical_Flux_Input_data_T destructing only the 'r'ight data members.
void destructor_Boundary_Value_c_data
	(struct Numerical_Flux_Input_c* num_flux_i ///< See brief.
		);

#endif // DPG__test_integration_support_face_h__INCLUDED
