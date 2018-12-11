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

#include <assert.h>
#include "petscsys.h"

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_test_integration.h"
#include "definitions_test_case.h"
#include "definitions_tol.h"

#include "test_base.h"
#include "test_integration.h"
#include "test_support_math_functions.h"
#include "test_support_solve.h"
#include "test_support_solve_dg.h"
#include "test_support_solve_dpg.h"
#include "test_support_solve_opg.h"

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
#include "compute_volume_rlhs_opg.h"
#include "compute_face_rlhs_opg.h"

// Static function declarations ************************************************************************************* //

///\{ \name Flag for whether PETSc Matrices should be output to files.
#define OUTPUT_PETSC_MATRICES false
///\}

/** \brief Function pointer to functions computing solution gradient coefficient related terms.
 *
 *  \param sim     \ref Simulation.
 *  \param volumes The list of volumes.
 *  \param faces   The list of faces.
 */
typedef void (*compute_grad_coef_fptr)
	(const struct Simulation* sim,
	 struct Intrusive_List*const volumes,
	 struct Intrusive_List*const faces
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
 *  \param sim        \ref Simulation.
 *  \param ssi        \ref Solver_Storage_Implicit.
 *  \param comp_elems The list of volumes/faces.
 */
typedef void (*compute_dof_coef_fptr)
	(const struct Simulation* sim,
	 struct Solver_Storage_Implicit* ssi,
	 struct Intrusive_List* comp_elems
	);

/** \brief Function pointer to functions computing the complex step linearization.
 *
 *  \param sim \ref Simulation.
 *  \param ssi \ref Solver_Storage_Implicit.
 */
typedef void (*compute_cmplx_step_fptr)
	(const struct Simulation* sim,
	 struct Solver_Storage_Implicit* ssi
	);

/// \brief Container for function pointers and data relating to the functions computing the linearization components.
struct F_Ptrs_and_Data {

	/** The solver method to be used for the complex step rhs evaluation. Note that some methods require Jacobians
	 *  even when only the rhs terms are being evaluated. */
	const char solver_method_cmplx;

	/// `derived_name` from \ref constructor_derived_Elements for the method under consideration.
	const int derived_elem_method;

	/// `derived_category` from \ref constructor_derived_computational_elements_T for the method under consideration.
	const int derived_comp_elem_method;

	compute_grad_coef_fptr compute_grad_coef;  ///< Solution gradient coefficients and related terms.
	compute_dof_coef_fptr  compute_volume_lhs; ///< Solution coefficient terms contributed by the volume term.
	compute_dof_coef_fptr  compute_face_lhs;   ///< Solution coefficient terms contributed by the face   term.
	compute_dof_coef_fptr  compute_all_lhs;    ///< Solution coefficient terms contributed by the all    terms.

	compute_cmplx_step_fptr compute_lhs_cmplx_step; ///< Compute the lhs terms using the complex step method.
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
	(struct Test_Info*const test_info,             ///< \ref Test_Info.
	 const struct Solver_Storage_Implicit* ssi[2], ///< \ref Solver_Storage_Implicit for the various methods.
	 const struct Simulation* sim                  ///< \ref Simulation.
	);

/// \brief Output PETSc matrices to file for visualization.
static void output_petsc_matrices
	(const struct Solver_Storage_Implicit* ssi[2] ///< \ref Solver_Storage_Implicit for the various methods.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the solver linearization (\ref test_integration_linearization.cpp).
 *  \return 0 on success (when the LHS matrices computed using complex step and exact linearization match).
 *
 *  The linearization is verified by comparing with the output when using the complex step method. Details of the
 *  complex step method can be found in Squire et al. \cite Squire1998 and Martins et al. \cite Martins2003.
 *
 *  For easier debugging, specific contributions to the linearization can be checked based on the value of \ref
 *  CHECK_LIN.
 *
 *  For second order equations, note that the volume rhs contributes off-diagonal terms to the global system matrix due
 *  to the use of the fully corrected weak gradient.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	assert_condition_message(argc == 2,"Invalid number of input arguments");
	const char* ctrl_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	const int p  = int_test_info->p_ref[0],
	          ml = int_test_info->ml[0],
	          p_prev  = p-1,
	          ml_prev = ml-1;

	const int adapt_type = int_test_info->adapt_type;
	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);

	struct Solver_Storage_Implicit* ssi[2] = { NULL, NULL, };
	for (int i = 0; i < 2; ++i) {
		const char type_rc = ( i == 0 ? 'r' : 'c' );

		struct Simulation* sim = NULL;
		structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,type_rc,false); // destructed

		struct F_Ptrs_and_Data* f_ptrs_data = constructor_F_Ptrs_and_Data(sim); // destructed

		if (i == 0) {
			ssi[i] = constructor_Solver_Storage_Implicit(sim); // destructed

			struct Test_Case* test_case = (struct Test_Case*) sim->test_case_rc->tc;
			test_case->solver_method_curr = 'i';
			const_cast_b(&test_case->use_schur_complement,false); // Otherwise A is modified.

			perturb_solution(sim);
			compute_lhs_analytical(sim,ssi[i],f_ptrs_data);
		} else {
			ssi[i] = constructor_Solver_Storage_Implicit_c(sim); // destructed

			struct Test_Case_c* test_case = (struct Test_Case_c*) sim->test_case_rc->tc;
			test_case->solver_method_curr = f_ptrs_data->solver_method_cmplx;

			perturb_solution(sim);
			compute_lhs_cmplx_step(sim,ssi[i],f_ptrs_data);
		}
		petsc_mat_vec_assemble(ssi[i]);

		destructor_F_Ptrs_and_Data(f_ptrs_data);

		if (i != 0) {
			if (OUTPUT_PETSC_MATRICES)
				output_petsc_matrices((const struct Solver_Storage_Implicit**)ssi);

			check_linearizations(&test_info,(const struct Solver_Storage_Implicit**)ssi,sim);
			for (int j = 0; j < 2; ++j)
				destructor_Solver_Storage_Implicit(ssi[j]);
		}
		structor_simulation(&sim,'d',adapt_type,p,ml,p_prev,ml_prev,NULL,type_rc,false);
	}

	destructor_Integration_Test_Info(int_test_info);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief 'Con'/'De'structor for derived solver elements of the appropriate type.
static void structor_derived_Elements
	(struct Simulation* sim,                    ///< \ref Simulation.
	 const struct F_Ptrs_and_Data* f_ptrs_data, ///< \ref F_Ptrs_and_Data.
	 const char mode                            ///< Mode of operation. Options: 'c'onstructor, 'd'estructor.
	);

static struct F_Ptrs_and_Data* constructor_F_Ptrs_and_Data (const struct Simulation* sim)
{
	struct F_Ptrs_and_Data* f_ptrs_data = malloc(sizeof *f_ptrs_data); // destructed

	switch (sim->method) {
	case METHOD_DG:
		const_cast_c(&f_ptrs_data->solver_method_cmplx,'e');

		const_cast_i(&f_ptrs_data->derived_elem_method,IL_ELEMENT_SOLVER_DG);
		const_cast_i(&f_ptrs_data->derived_comp_elem_method,IL_SOLVER_DG);

		f_ptrs_data->compute_grad_coef  = compute_grad_coef_dg;
		f_ptrs_data->compute_volume_lhs = compute_volume_rlhs_dg;
		f_ptrs_data->compute_face_lhs   = compute_face_rlhs_dg;

		f_ptrs_data->compute_lhs_cmplx_step = compute_lhs_cmplx_step_dg;
		break;
	case METHOD_DPG:
		const_cast_c(&f_ptrs_data->solver_method_cmplx,'i'); // Jacobians are required.

		const_cast_i(&f_ptrs_data->derived_elem_method,IL_ELEMENT_SOLVER_DPG);
		const_cast_i(&f_ptrs_data->derived_comp_elem_method,IL_SOLVER_DPG);

		f_ptrs_data->compute_all_lhs = compute_all_rlhs_dpg;

		f_ptrs_data->compute_lhs_cmplx_step = compute_lhs_cmplx_step_dpg;
		break;
	case METHOD_OPG:
		const_cast_c(&f_ptrs_data->solver_method_cmplx,'e');

		const_cast_i(&f_ptrs_data->derived_elem_method,IL_ELEMENT_SOLVER_OPG);
		const_cast_i(&f_ptrs_data->derived_comp_elem_method,IL_SOLVER_OPG);

		f_ptrs_data->compute_volume_lhs = compute_volume_rlhs_opg;
		f_ptrs_data->compute_face_lhs   = compute_face_rlhs_opg;

		f_ptrs_data->compute_lhs_cmplx_step = compute_lhs_cmplx_step_opg;
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
	structor_derived_Elements(sim,f_ptrs_data,'c'); // destructed
	constructor_derived_computational_elements(sim,f_ptrs_data->derived_comp_elem_method); // destructed
	initialize_zero_memory_volumes(sim->volumes);
	switch (sim->method) {
	case METHOD_DG:
		f_ptrs_data->compute_grad_coef(sim,sim->volumes,sim->faces);
		switch (CHECK_LIN) {
		case CHECK_LIN_VOLUME:
			f_ptrs_data->compute_volume_lhs(sim,ssi,sim->volumes);
			break;
		case CHECK_LIN_FACE:
			f_ptrs_data->compute_face_lhs(sim,ssi,sim->faces);
			break;
		case CHECK_LIN_ALL:
			f_ptrs_data->compute_volume_lhs(sim,ssi,sim->volumes);
			f_ptrs_data->compute_face_lhs(sim,ssi,sim->faces);
			break;
		default:
			EXIT_ERROR("Unsupported: %d.\n",CHECK_LIN);
			break;
		}
		break;
	case METHOD_DPG:
		f_ptrs_data->compute_all_lhs(sim,ssi,sim->volumes);
		break;
	case METHOD_OPG:
		switch (CHECK_LIN) {
		case CHECK_LIN_VOLUME:
			f_ptrs_data->compute_volume_lhs(sim,ssi,sim->volumes);
			break;
		case CHECK_LIN_FACE:
			f_ptrs_data->compute_face_lhs(sim,ssi,sim->faces);
			break;
		case CHECK_LIN_ALL:
			f_ptrs_data->compute_volume_lhs(sim,ssi,sim->volumes);
			f_ptrs_data->compute_face_lhs(sim,ssi,sim->faces);
			break;
		default:
			EXIT_ERROR("Unsupported: %d.\n",CHECK_LIN);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",sim->method);
		break;
	}
	structor_derived_Elements(sim,f_ptrs_data,'d');
	destructor_derived_computational_elements(sim,IL_SOLVER);
}

static void compute_lhs_cmplx_step
	(struct Simulation* sim, struct Solver_Storage_Implicit* ssi, const struct F_Ptrs_and_Data* f_ptrs_data)
{
	structor_derived_Elements(sim,f_ptrs_data,'c'); // destructed
	constructor_derived_computational_elements_c(sim,f_ptrs_data->derived_comp_elem_method); // destructed
	f_ptrs_data->compute_lhs_cmplx_step(sim,ssi);
	structor_derived_Elements(sim,f_ptrs_data,'d');
	destructor_derived_computational_elements_c(sim,IL_SOLVER);
}

static void check_linearizations
	(struct Test_Info*const test_info, const struct Solver_Storage_Implicit* ssi[2], const struct Simulation* sim)
{
	UNUSED(test_info);
	bool pass        = true;
	const double tol = 4e-13;

	const double diff = norm_diff_petsc_Mat(ssi[0]->A,ssi[1]->A,false);
	if (diff > tol) {
		pass = false;
		printf("% .3e (tol = % .3e).\n",diff,tol);
		expect_condition(pass,"difference");
	}

	if (check_symmetric(sim) && CHECK_LIN == CHECK_LIN_ALL) {
		// Note: Symmetry is lost when terms are omitted from the full linearization.
		const double tol_symm = 1e-12;
		PetscBool symmetric[2] = { false, false, };
		for (int i = 0; i < 2; ++i)
			MatIsSymmetric(ssi[i]->A,tol_symm,&symmetric[i]);
		if (!(symmetric[0] && symmetric[1])) {
			pass = false;
			printf("symm: %d %d (tol = % .3e).\n",symmetric[0],symmetric[1],tol_symm);
			expect_condition(pass,"symmetric");
		}
	}

	if (pass && (CHECK_LIN != CHECK_LIN_ALL)) {
		pass = false;
		printf("Passing, but not checking full linearization (CHECK_LIN: %d).\n",CHECK_LIN);
		expect_condition(pass,"check full linearization");
	}

	assert_condition(pass);
}

static void output_petsc_matrices (const struct Solver_Storage_Implicit* ssi[2])
{
	PetscViewer viewer;

	PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat_output_analytical.m",&viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	MatView(ssi[0]->A,viewer);
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat_output_cmplx_step.m",&viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	MatView(ssi[1]->A,viewer);

	PetscViewerDestroy(&viewer);
}

// Level 1 ********************************************************************************************************** //

static void structor_derived_Elements
	(struct Simulation* sim, const struct F_Ptrs_and_Data* f_ptrs_data, const char mode)
{
	assert(mode == 'c' || mode == 'd');

	switch (mode) {
	case 'c':
		constructor_derived_Elements(sim,f_ptrs_data->derived_elem_method); // destructed
		break;
	case 'd':
		destructor_derived_Elements(sim,IL_ELEMENT_SOLVER);
		break;
	default:
		EXIT_ERROR("Unsupported: %c",mode);
		break;
	}
}
