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
 *  \brief Provides the interface for the templated \ref OPG_Solver_Face_T container and associated functions.
 *
 *  These faces are needed by the 'O'ptimal 'P'etrov 'G'alerkin solver functions.
 */

#include "def_templates_face_solver.h"
#include "def_templates_face_solver_opg.h"
#include "def_templates_matrix.h"

/// \brief Container for data relating to the OPG solver faces.
struct OPG_Solver_Face_T {
	struct Solver_Face_T face; ///< The base \ref Solver_Face_T.

	const struct const_Matrix_R* m_inv; ///< The inverse mass matrix.
};

/// \brief Constructor for a derived \ref OPG_Solver_Face_T.
void constructor_derived_OPG_Solver_Face_T
	(struct Face* face_ptr,       ///< Pointer to the face.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref OPG_Solver_Face_T.
void destructor_derived_OPG_Solver_Face_T
	(struct Face* face_ptr ///< Pointer to the face.
	);

/** \brief Get the pointer to the appropriate \ref OPG_Solver_Element::cv1_vt_fc operator.
 *  \return See brief. */
struct Multiarray_Operator get_operator__cv1_vt_fc_T
	(const int side_index,                           ///< The index of the side of the face under consideration.
	 const struct OPG_Solver_Face_T*const opg_s_face ///< Standard.
	);

#include "undef_templates_face_solver.h"
#include "undef_templates_face_solver_opg.h"
#include "undef_templates_matrix.h"
