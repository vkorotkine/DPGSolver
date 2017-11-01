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

#include "test_integration_linearization.h"

#include <assert.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "test_integration.h"

#include "face.h"
#include "element.h"
#include "volume.h"

#include "computational_elements.h"
#include "geometry.h"
#include "intrusive.h"
#include "mesh.h"
#include "simulation.h"
#include "solution.h"
#include "solve.h"
#include "solve_implicit.h"
#include "test_case.h"

#include "compute_grad_coef_dg.h"
#include "compute_volume_rlhs_dg.h"
#include "compute_face_rlhs_dg.h"

// Static function declarations ************************************************************************************* //

///\{ \name The available options and value for the linearization components to check.
#define CHECK_LIN_VOLUME 1 ///< Volume contributions only.
#define CHECK_LIN_FACE   2 ///< Face contributions only.
#define CHECK_LIN_ALL    3 ///< All contributions.

#define CHECK_LIN 1 ///< The selected option from those listed above.
///\}

/** \brief Function pointer to functions computing solution gradient coefficient related terms.
 *
 *  \param sim \ref Simulation.
 */
typedef void (*compute_grad_coef_fptr)
	(const struct Simulation* sim
	);

/** \brief Function pointer to functions computing solution coefficient related terms.
 *
 *  \param sim       \ref Simulation.
 *  \param s_store_i \ref Solver_Storage_Implicit.
 */
typedef void (*compute_sol_coef_fptr)
	(const struct Simulation* sim,
	 struct Solver_Storage_Implicit* s_store_i
	);

/// \brief Container for the function pointers to the functions computing the linearization components.
struct F_Ptrs {
	compute_grad_coef_fptr compute_grad_coef;  ///< Solution gradient coefficients and related terms.
	compute_sol_coef_fptr  compute_volume_lhs; ///< Solution coefficient terms contributed by the volume term.
	compute_sol_coef_fptr  compute_face_lhs;   ///< Solution coefficient terms contributed by the face   term.
};

/// \brief Construcor for the derived elements and computational elements for the specific solver method.
static void constructor_derived_solver
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destrucor for the derived elements and computational elements for the specific solver method.
static void destructor_derived_solver
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for a \ref F_Ptrs container, set based on \ref Simulation::method.
 *  \return See brief. */
static struct F_Ptrs* constructor_F_Ptrs
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref F_Ptrs container.
static void destructor_F_Ptrs
	(struct F_Ptrs* f_ptrs ///< Standard.
	);

/// \brief Compute \ref Solver_Storage_Implicit::A using the analytical linearization.
static void compute_lhs_analytical
	(const struct Simulation* sim,        ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi, ///< \ref Solver_Storage_Implicit.
	 const struct F_Ptrs* f_ptrs          ///< \ref F_Ptrs.
	);

/// \brief Compute \ref Solver_Storage_Implicit::A using the complex step method.
static void compute_lhs_cmplx_step
	(const struct Simulation* sim,       ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi, ///< \ref Solver_Storage_Implicit.
	 const struct F_Ptrs* f_ptrs          ///< \ref F_Ptrs.
	);

// Interface functions ********************************************************************************************** //

void test_integration_linearization (const char*const ctrl_name)
{
	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	const int p  = int_test_info->p_ref[1],
	          ml = int_test_info->ml[1];

	const int adapt_type = int_test_info->adapt_type;
	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,ctrl_name);
	struct Simulation* sim = NULL;
	structor_simulation(&sim,'c',adapt_type,p,ml,p,ml,ctrl_name_curr); // destructed

	constructor_derived_computational_elements(sim,IL_SOLVER); // destructed

	set_up_solver_geometry(sim);
	set_initial_solution(sim);

	sim->test_case->solver_method_curr = 'i';
	constructor_derived_solver(sim); // destructed

	struct Solver_Storage_Implicit* ssi[2] = { constructor_Solver_Storage_Implicit(sim),
	                                           constructor_Solver_Storage_Implicit(sim), }; // destructed

	struct F_Ptrs* f_ptrs = constructor_F_Ptrs(sim); // destructed

	compute_lhs_analytical(sim,ssi[0],f_ptrs);
	compute_lhs_cmplx_step(sim,ssi[0],f_ptrs);

	destructor_F_Ptrs(f_ptrs);

	for (int i = 0; i < 2; ++i)
		petsc_mat_vec_assemble(ssi[i]);

// test here.
EXIT_UNSUPPORTED;

	for (int i = 0; i < 2; ++i)
		destructor_Solver_Storage_Implicit(ssi[i]);

	destructor_derived_solver(sim);

	destructor_derived_computational_elements(sim,IL_BASE);

	structor_simulation(&sim,'d',adapt_type,p,ml,p,ml,NULL);

	destructor_Integration_Test_Info(int_test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_solver (struct Simulation* sim)
{
	int derived_element      = -1,
	    derived_comp_element = -1;
	switch (sim->method) {
	case METHOD_DG:
		derived_element      = IL_ELEMENT_SOLVER_DG;
		derived_comp_element = IL_SOLVER_DG;
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",sim->method);
		break;
	}
	constructor_derived_Elements(sim,derived_element);                    // destructed
	constructor_derived_computational_elements(sim,derived_comp_element); // destructed
}

static void destructor_derived_solver (struct Simulation* sim)
{
	destructor_derived_computational_elements(sim,IL_SOLVER);
	destructor_derived_Elements(sim,IL_ELEMENT);
}

static struct F_Ptrs* constructor_F_Ptrs (const struct Simulation* sim)
{
	struct F_Ptrs* f_ptrs = malloc(sizeof *f_ptrs); // destructed

	switch (sim->method) {
	case METHOD_DG:
		f_ptrs->compute_grad_coef  = compute_grad_coef_dg;
		f_ptrs->compute_volume_lhs = compute_volume_rlhs_dg;
		f_ptrs->compute_face_lhs   = compute_face_rlhs_dg;
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",sim->method);
		break;
	}

	return f_ptrs;
}

static void destructor_F_Ptrs (struct F_Ptrs* f_ptrs)
{
	free(f_ptrs);
}

static void compute_lhs_analytical
	(const struct Simulation* sim, struct Solver_Storage_Implicit* ssi, const struct F_Ptrs* f_ptrs)
{
	switch (CHECK_LIN) {
	case CHECK_LIN_VOLUME:
		f_ptrs->compute_volume_lhs(sim,ssi);
		break;
	case CHECK_LIN_FACE:
		f_ptrs->compute_face_lhs(sim,ssi);
		break;
	case CHECK_LIN_ALL:
		f_ptrs->compute_volume_lhs(sim,ssi);
		f_ptrs->compute_face_lhs(sim,ssi);
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",CHECK_LIN);
		break;
	}
}

static void compute_lhs_cmplx_step
	(const struct Simulation* sim, struct Solver_Storage_Implicit* ssi, const struct F_Ptrs* f_ptrs)
{
// To avoid polluting the rest of the code by linking the Test_Support library, defined the derived complex solver
// volumes/faces in the directory of the corresponding solver_method. Con/De'structor functions should also be defined
// there. This also applies for fluxes/boundary conditions/numerical fluxes.

// All other functions (such as compute_lhs_cmplx functions) can be defined in test_support_* files.

UNUSED(sim);
UNUSED(ssi);
UNUSED(f_ptrs);
EXIT_ADD_SUPPORT;
}
