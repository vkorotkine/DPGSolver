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
#include "definitions_test_integration.h"
#include "definitions_tol.h"

#include "test_base.h"
#include "test_integration.h"
#include "test_support_math_functions.h"
#include "test_complex_solve_dg.h"
#include "test_complex_solve_dpg.h"

#include "face.h"
#include "element.h"
#include "volume.h"
#include "volume_solver.h"

#include "computational_elements.h"
#include "const_cast.h"
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
#include "compute_all_rlhs_dpg.h"

// Static function declarations ************************************************************************************* //

///\{ \name Flag for whether PETSc Matrices should be output to files.
#define OUTPUT_PETSC_MATRICES false
///\}

/** \brief Function pointer to functions perturbing the solution to ensure that all terms in the linearization are
 *         activated.
 *
 *  \param sim \ref Simulation.
 */
typedef void (*perturb_solution_fptr)
	(const struct Simulation* sim
	);

/** \brief Function pointer to functions computing solution gradient coefficient related terms.
 *
 *  \param sim \ref Simulation.
 */
typedef void (*compute_grad_coef_fptr)
	(const struct Simulation* sim
	);

/** \brief Function pointer to function setting the initial solution for the complex solver computational elements.
 *
 *  \param sim \ref Simulation.
 */
typedef void (*set_initial_solution_complex_fptr)
	(const struct Simulation* sim
	);

/** \brief Function pointer to functions computing global degree of freedom coefficient related terms.
 *
 *  \param sim       \ref Simulation.
 *  \param s_store_i \ref Solver_Storage_Implicit.
 */
typedef void (*compute_dof_coef_fptr)
	(const struct Simulation* sim,
	 struct Solver_Storage_Implicit* s_store_i
	);

/// \brief Container for function pointers and data relating to the functions computing the linearization components.
struct F_Ptrs_and_Data {
	perturb_solution_fptr perturb_solution; ///< \ref perturb_solution_fptr.

	/// `derived_name` from \ref constructor_derived_Elements.
	const int derived_elem;

	/// `derived_category` from \ref constructor_derived_computational_elements for the analytical contributions.
	const int derived_comp_elem_analytical;

	/// `derived_category` from \ref constructor_derived_computational_elements for the complex step contributions.
	const int derived_comp_elem_cmplx_step;

	compute_grad_coef_fptr compute_grad_coef;  ///< Solution gradient coefficients and related terms.
	compute_dof_coef_fptr  compute_volume_lhs; ///< Solution coefficient terms contributed by the volume term.
	compute_dof_coef_fptr  compute_face_lhs;   ///< Solution coefficient terms contributed by the face   term.
	compute_dof_coef_fptr  compute_all_lhs;    ///< Solution coefficient terms contributed by the all    terms.

	set_initial_solution_complex_fptr set_initial_solution_complex; ///< \ref set_initial_solution_complex_fptr.
	compute_dof_coef_fptr compute_lhs_cmplx_step; ///< Compute the lhs terms using the complex step method.
};

/** \brief Constructor for a \ref F_Ptrs_and_Data container, set based on \ref Simulation::method.
 *  \return See brief. */
static struct F_Ptrs_and_Data* constructor_F_Ptrs_and_Data
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref F_Ptrs_and_Data container.
static void destructor_F_Ptrs_and_Data
	(struct F_Ptrs_and_Data* f_ptrs_data ///< Standard.
	);

/// \brief Compute \ref Solver_Storage_Implicit::A using the analytical linearization.
static void compute_lhs_analytical
	(struct Simulation* sim,                   ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi,      ///< \ref Solver_Storage_Implicit.
	 const struct F_Ptrs_and_Data* f_ptrs_data ///< \ref F_Ptrs_and_Data.
	);

/// \brief Compute \ref Solver_Storage_Implicit::A using the complex step method.
static void compute_lhs_cmplx_step
	(struct Simulation* sim,                   ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi,      ///< \ref Solver_Storage_Implicit.
	 const struct F_Ptrs_and_Data* f_ptrs_data ///< \ref F_Ptrs_and_Data.
	);

/// \brief Check the linearizations computed with the various methods.
static void check_linearizations
	(struct Test_Info*const test_info,                  ///< \ref Test_Info.
	 const struct Integration_Test_Info* int_test_info, ///< \ref Integration_Test_Info.
	 const struct Solver_Storage_Implicit* ssi[2],      ///< \ref Solver_Storage_Implicit for the various methods.
	 const struct Simulation* sim                       ///< \ref Simulation.
	);

/// \brief Output PETSc matrices to file for visualization.
static void output_petsc_matrices
	(const struct Solver_Storage_Implicit* ssi[2] ///< \ref Solver_Storage_Implicit for the various methods.
	);

// Interface functions ********************************************************************************************** //

void test_integration_linearization (struct Test_Info*const test_info, const char*const ctrl_name)
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

	struct F_Ptrs_and_Data* f_ptrs_data = constructor_F_Ptrs_and_Data(sim); // destructed
	f_ptrs_data->perturb_solution(sim);

	sim->test_case->solver_method_curr = 'i';
	constructor_derived_Elements(sim,IL_ELEMENT_SOLVER);         // destructed
	constructor_derived_Elements(sim,f_ptrs_data->derived_elem); // destructed

	struct Solver_Storage_Implicit* ssi[2] = { constructor_Solver_Storage_Implicit(sim),
	                                           constructor_Solver_Storage_Implicit(sim), }; // destructed

	compute_lhs_analytical(sim,ssi[0],f_ptrs_data);
	compute_lhs_cmplx_step(sim,ssi[1],f_ptrs_data);

	if (OUTPUT_PETSC_MATRICES)
		output_petsc_matrices((const struct Solver_Storage_Implicit**)ssi);

	destructor_F_Ptrs_and_Data(f_ptrs_data);

	for (int i = 0; i < 2; ++i)
		petsc_mat_vec_assemble(ssi[i]);

	check_linearizations(test_info,int_test_info,(const struct Solver_Storage_Implicit**)ssi,sim);

	for (int i = 0; i < 2; ++i)
		destructor_Solver_Storage_Implicit(ssi[i]);

	destructor_derived_Elements(sim,IL_ELEMENT_SOLVER);
	destructor_derived_Elements(sim,IL_ELEMENT);

	destructor_derived_computational_elements(sim,IL_BASE);

	structor_simulation(&sim,'d',adapt_type,p,ml,p,ml,NULL);

	destructor_Integration_Test_Info(int_test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct F_Ptrs_and_Data* constructor_F_Ptrs_and_Data (const struct Simulation* sim)
{
	struct F_Ptrs_and_Data* f_ptrs_data = malloc(sizeof *f_ptrs_data); // destructed

	switch (sim->method) {
	case METHOD_DG:
		f_ptrs_data->perturb_solution = perturb_solution_dg;

		const_cast_i(&f_ptrs_data->derived_elem,IL_ELEMENT_SOLVER_DG);
		const_cast_i(&f_ptrs_data->derived_comp_elem_analytical,IL_SOLVER_DG);
		const_cast_i(&f_ptrs_data->derived_comp_elem_cmplx_step,IL_SOLVER_DG_COMPLEX);

		f_ptrs_data->compute_grad_coef  = compute_grad_coef_dg;
		f_ptrs_data->compute_volume_lhs = compute_volume_rlhs_dg;
		f_ptrs_data->compute_face_lhs   = compute_face_rlhs_dg;

		f_ptrs_data->set_initial_solution_complex = set_initial_solution_complex_dg;
		f_ptrs_data->compute_lhs_cmplx_step       = compute_lhs_cmplx_step_dg;
		break;
	case METHOD_DPG:
		f_ptrs_data->perturb_solution = perturb_solution_dpg;

		const_cast_i(&f_ptrs_data->derived_elem,IL_ELEMENT_SOLVER_DPG);
		const_cast_i(&f_ptrs_data->derived_comp_elem_analytical,IL_SOLVER_DPG);
		const_cast_i(&f_ptrs_data->derived_comp_elem_cmplx_step,IL_SOLVER_DPG_COMPLEX);

		f_ptrs_data->compute_all_lhs = compute_all_rlhs_dpg;

		f_ptrs_data->set_initial_solution_complex = set_initial_solution_complex_dpg;
		f_ptrs_data->compute_lhs_cmplx_step       = compute_lhs_cmplx_step_dpg;
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",sim->method);
		break;
	}

	return f_ptrs_data;
}

static void destructor_F_Ptrs_and_Data (struct F_Ptrs_and_Data* f_ptrs_data)
{
	free(f_ptrs_data);
}

static void compute_lhs_analytical
	(struct Simulation* sim, struct Solver_Storage_Implicit* ssi, const struct F_Ptrs_and_Data* f_ptrs_data)
{
	constructor_derived_computational_elements(sim,f_ptrs_data->derived_comp_elem_analytical); // destructed
	switch (sim->method) {
	case METHOD_DG:
		switch (CHECK_LIN) {
		case CHECK_LIN_VOLUME:
			f_ptrs_data->compute_volume_lhs(sim,ssi);
			break;
		case CHECK_LIN_FACE:
			f_ptrs_data->compute_face_lhs(sim,ssi);
			break;
		case CHECK_LIN_ALL:
			f_ptrs_data->compute_volume_lhs(sim,ssi);
			f_ptrs_data->compute_face_lhs(sim,ssi);
			break;
		default:
			EXIT_ERROR("Unsupported: %d.\n",CHECK_LIN);
			break;
		}
		break;
	case METHOD_DPG:
		f_ptrs_data->compute_all_lhs(sim,ssi);
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",sim->method);
		break;
	}
	destructor_derived_computational_elements(sim,IL_SOLVER);

	petsc_mat_vec_assemble(ssi);
}

static void compute_lhs_cmplx_step
	(struct Simulation* sim, struct Solver_Storage_Implicit* ssi, const struct F_Ptrs_and_Data* f_ptrs_data)
{
	constructor_derived_computational_elements(sim,f_ptrs_data->derived_comp_elem_analytical); // destructed
	constructor_derived_computational_elements(sim,f_ptrs_data->derived_comp_elem_cmplx_step); // destructed

	f_ptrs_data->set_initial_solution_complex(sim);
	f_ptrs_data->compute_lhs_cmplx_step(sim,ssi);

	destructor_derived_computational_elements(sim,f_ptrs_data->derived_comp_elem_analytical);
	destructor_derived_computational_elements(sim,IL_SOLVER);
}

static void check_linearizations
	(struct Test_Info*const test_info, const struct Integration_Test_Info* int_test_info,
	 const struct Solver_Storage_Implicit* ssi[2], const struct Simulation* sim)
{
	bool pass        = true;
	const double tol = 2e1*EPS;

	const double diff = norm_diff_petsc_Mat(ssi[0]->A,ssi[1]->A);
	if (diff > tol) {
		pass = false;
		printf("Failed difference: % .3e (tol = % .3e).\n",diff,tol);
	}

	if (check_symmetric(sim)) {
		PetscBool symmetric = false;
		MatIsSymmetric(ssi[0]->A,tol,&symmetric);
		if (!symmetric) {
			pass = false;
			printf("Failed symmetric (tol = % .3e).\n",tol);
		}
	}

	char test_name[STRLEN_MAX];
	sprintf(test_name,"%s%s","Linearization - ",int_test_info->ctrl_name);
	test_increment_and_print_name(test_info,pass,test_name);
}

static void output_petsc_matrices (const struct Solver_Storage_Implicit* ssi[2])
{
	PetscViewer viewer;

	PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat_output0.m",&viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	MatView(ssi[0]->A,viewer);
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat_output1.m",&viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	MatView(ssi[1]->A,viewer);

	PetscViewerDestroy(&viewer);
}
