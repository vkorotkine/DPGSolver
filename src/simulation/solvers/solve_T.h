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
 *  \brief Provides the interface to templated functions used to solve for the solution.
 */

#include <stddef.h>

#include "def_templates_solve.h"
#include "def_templates_multiarray.h"
#include "def_templates_volume_solver.h"

struct Simulation;
struct const_Multiarray_T;
struct Solver_Volume_T;
struct Intrusive_List;

/** \brief Constructor for a \ref Solver_Storage_Implicit container.
 *  \return See brief. */
struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit_T
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom.
 *  \return See brief. */
ptrdiff_t compute_dof_T
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Add the source contribution to \ref Solver_Volume_T::flux_imbalance.
void add_to_flux_imbalance_source_T
	(const struct const_Multiarray_T*const source_vc_w_J, /**< The source at the volume cubature nodes with Jacobian
	                                                       *   scaling. */
	 const struct Solver_Volume_T*const s_vol,            ///< \ref Solver_Volume_T.
	 const struct Simulation*const sim                    ///< \ref Simulation.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::tw0_vt_vc operator.
 *  \return See brief. */
const struct Operator* get_operator__tw0_vt_vc_T
	(const struct Solver_Volume_T* s_vol ///< The current volume.
	);

/// \brief Set the memory of the rhs and lhs (if applicable) terms to zero for the volumes.
void initialize_zero_memory_volumes_T
	(struct Intrusive_List* volumes ///< The list of volumes for which to set the memory.
	);

/** \brief Update the global 'd'egree 'o'f 'f'reedom indices in \ref Solver_Volume_T and \ref Solver_Face_T based on the
 *         size of the allocated variable containers. */
void update_ind_dof_T
	(const struct Simulation*const sim ///< Standard.
	 );

#include "undef_templates_solve.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_volume_solver.h"
