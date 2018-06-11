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
 *  \brief Provides the interface to functions used to set the default steady solution for diffusion test cases.
 */

#include "def_templates_solution.h"
#include "def_templates_solution_diffusion.h"
#include "def_templates_volume_solver.h"
#include "def_templates_multiarray.h"

struct Simulation;
struct Solution_Container_T;
struct Solver_Volume_T;
struct Multiarray_T;
struct const_Multiarray_R;

/// \brief Function to be used for \ref Test_Case_T::set_sol for the default steady diffusion solution.
void set_sol_diffusion_default_steady_T
	(const struct Simulation* sim,        ///< Defined for \ref set_sol_fptr_T.
	 struct Solution_Container_T sol_cont ///< Defined for \ref set_sol_fptr_T.
	);

/// \brief Function to be used for \ref Test_Case_T::set_grad for the default steady diffusion solution.
void set_grad_diffusion_default_steady_T
	(const struct Simulation* sim,        ///< Defined for \ref set_sol_fptr_T.
	 struct Solution_Container_T sol_cont ///< Defined for \ref set_sol_fptr_T.
	);

/** \brief \ref Test_Case_T::constructor_sol for the default steady diffusion solution.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_const_sol_diffusion_default_steady_T
	(const struct const_Multiarray_R* xyz, ///< Defined for \ref constructor_sol_fptr_T.
	 const struct Simulation* sim          ///< Defined for \ref constructor_sol_fptr_T.
	);

/** \brief \ref Test_Case_T::constructor_grad for the default steady diffusion solution gradient.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_const_grad_diffusion_default_steady_T
	(const struct const_Multiarray_R* xyz, ///< Defined for \ref constructor_sol_fptr_T.
	 const struct Simulation* sim          ///< Defined for \ref constructor_sol_fptr_T.
	);

/// \brief Version of \ref compute_source_rhs_fptr_T for the default steady diffusion solution.
void compute_source_rhs_diffusion_default_steady_T
	(const struct Simulation* sim,        ///< See brief.
	 const struct Solver_Volume_T* s_vol, ///< See brief.
	 struct Multiarray_T* rhs             ///< See brief.
	);

/** \brief Version of \ref compute_source_rhs_fptr_T adding flux imbalance contributions for the default steady diffusion
 *         solution. */
void add_to_flux_imbalance_source_diffusion_default_steady_T
	(const struct Simulation* sim,        ///< See brief.
	 const struct Solver_Volume_T* s_vol, ///< See brief.
	 struct Multiarray_T* rhs             ///< See brief.
	);

#include "undef_templates_solution.h"
#include "undef_templates_solution_diffusion.h"
#include "undef_templates_volume_solver.h"
#include "undef_templates_multiarray.h"
