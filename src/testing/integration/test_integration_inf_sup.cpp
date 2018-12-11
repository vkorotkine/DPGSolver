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
#include "slepceps.h"
#include "definitions_slepc.h"

#include "macros.h"
#include "definitions_adaptation.h"
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

#include "vector.h"

#include "computational_elements.h"
#include "const_cast.h"
#include "file_processing.h"
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

#define OUTPUT_PETSC_MATRICES_IS true  ///< Flag for whether PETSc Matrices should be output to files.
#define PRINT_EIG_VALS_ERRORS    false ///< Flag for whether eigenvalues and errors should be printed.

/** \brief Constructor for a \ref Gen_Eig_Data container holding the operators for the solver method under consideration.
 *  \return See brief.
 *
 *  Please see the comments in the main function of \ref test_integration_inf_sup.cpp for the details of the lhs and rhs
 *  operators which are being assembled here.
 */
static const struct Gen_Eig_Data* constructor_Gen_Eig_Data_inf_sup
	(struct Simulation*const sim ///< Standard.
	 );

/// \brief Output the data in the input \ref Gen_Eig_Data container to files.
static void output_Gen_Eig_Data
	(const struct Gen_Eig_Data*const ged, ///< Standard.
	 const int ml                         ///< The mesh level
	 );

/// \brief Compute and store the minimum eigenvalue for the generalized eigenvalue problem.
static void compute_and_store_min_eigenvalue
	(const struct Gen_Eig_Data*const ged, ///< Standard.
	 double*const eig_val_ptr             ///< Pointer to the memory location to store the eigenvalue.
	 );

/** \brief Set \ref Integration_Test_Info::conv_order_inf_sup_accept to the value specified for the test
 *         case (from the input file) or to 0.0 otherwise. */
static void set_convergence_order_inf_sup_accept
	(struct Integration_Test_Info*const int_test_info ///< \ref Integration_Test_Info.
	 );

/// \brief Check the inf-sup constants for the sequence of simulations.
static void check_inf_sup
	(bool*const pass,                                        ///< To be set based on the result of the test.
	 struct Test_Info*const test_info,                       ///< Standard.
	 const struct Integration_Test_Info*const int_test_info, ///< Standard.
	 const struct Vector_d*const eig_vals,                   ///< Vector of eigenvalues.
	 const struct Simulation*const sim                       ///< Standard.
	 );

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing relating to inf-sup constants for solver methods
 *        (\ref test_integration_inf_sup.cpp).
 *  \return 0 on success (when the inf-sup constant is bounded below with variation of the parameter under
 *          consideration).
 *
 *  This test computes the inf-sup constant for the solver method specified in the control file on a sequence of nested
 *  meshes with all other parameters held fixed or on the same mesh with variation of a parameter known to be
 *  problematic (e.g. tending towards zero viscosity). The inf-sup constants are computed following the method outlined
 *  by Bathe et al. (eq. (14), \cite Bathe2000) by solving the eigenvalue problem based on the bilinear form and the
 *  norms for the trial and test spaces:
 *  \f[
 *      (A' T^{-1} A) v = \lambda S v
 *  \f]
 *
 *  where \f$ \gamma = \sqrt{\lambda_{min}} \f$.
 *
 *  In the case of the method not being perfectly inf-sup stable, the test will pass with a warning if the rate at which
 *  the constant tends to zero is within \ref Integration_Test_Info::conv_order_inf_sup_accept. In general, that the
 *  inf-sup constant is bounded away from zero is a necessary condition for convergence (eq. (3) Bathe et al.
 *  \cite Bathe2000). However, if it tends to zero at a sufficiently slow rate or only tends to zero for relatively
 *  coarse meshes, the method may still be stable depending on the specified parameters.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	SlepcInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	assert_condition_message(argc == 3,"Invalid number of input arguments");
	const char* ctrl_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name); // destructed

	const int* p_ref  = int_test_info->p_ref;
	const int* ml_ref = int_test_info->ml;

	struct Simulation* sim = NULL;
	const char type_rc = 'r';

	int ml_prev = ml_ref[0]-1;
	int p_prev  = p_ref[0]-1;

	assert(int_test_info->adapt_type == ADAPT_0 || int_test_info->adapt_type == ADAPT_H);

	struct Vector_d*const eig_vals = constructor_empty_Vector_d(ml_ref[1]+1); // destructed
	set_to_value_Vector_d(eig_vals,-1.0);
	for (int ml_max = ml_ref[1], p = p_ref[0], ml = ml_ref[0]; ml <= ml_max; ++ml) {
		const int adapt_type = int_test_info->adapt_type;
		const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);
		structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,type_rc,false); // d.

		set_initial_v_test_sg_coef(sim); // Needed for dof indices for the test space.
		assert(!test_case_explicitly_enforces_conservation(sim)); // Ensure that all is working as expected.

		const struct Gen_Eig_Data*const ged = constructor_Gen_Eig_Data_inf_sup(sim); // destructed

		if (OUTPUT_PETSC_MATRICES_IS)
			output_Gen_Eig_Data(ged,ml);

		compute_and_store_min_eigenvalue(ged,&eig_vals->data[ml]);
		destructor_Gen_Eig_Data(ged);

		if ((ml == ml_max) && (p == p_ref[1])) {
			set_convergence_order_inf_sup_accept(int_test_info);
			bool pass = true;
			check_inf_sup(&pass,&test_info,int_test_info,eig_vals,sim);
			assert_condition(pass);

			structor_simulation(&sim,'d',ADAPT_0,p,ml,p_prev,ml_prev,NULL,type_rc,false);
		}

		ml_prev = ml;
	}
	destructor_Vector_d(eig_vals);
	destructor_Integration_Test_Info(int_test_info);

	SlepcFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const struct Gen_Eig_Data* constructor_Gen_Eig_Data_inf_sup (struct Simulation*const sim)
{
	const struct Gen_Eig_Data* ged = NULL;

	switch (get_set_method(NULL)) {
	case METHOD_DG:    ged = constructor_Gen_Eig_Data_inf_sup_dg(sim); break;
	case METHOD_DPG:   ged = constructor_Gen_Eig_Data_inf_sup_dpg(sim); break;
	/* case METHOD_OPGC0: ged = constructor_Gen_Eig_Data_inf_sup_opg(sim); break; */
	default: EXIT_ERROR("Unsupported: %d\n",get_set_method(NULL)); break;
	}

	return ged;
}

static void output_Gen_Eig_Data (const struct Gen_Eig_Data*const ged, const int ml)
{
	printf("Outputting Generalized Eigenvalue Data. Set OUTPUT_PETSC_MATRICES_IS to false for efficiency if not used.\n");

	char name_out[STRLEN_MIN];
	sprintf(name_out,"%s%d%s","A_ged_",ml,".m");
	output_petsc_mat(ged->A,name_out);

	sprintf(name_out,"%s%d%s","B_ged_",ml,".m");
	output_petsc_mat(ged->B,name_out);
}

static void compute_and_store_min_eigenvalue (const struct Gen_Eig_Data*const ged, double*const eig_val_ptr)
{
	SlepcEPS eps;
	EPSCreate(PETSC_COMM_WORLD,&eps);
	EPSSetOperators(eps,ged->A,ged->B);
	EPSSetProblemType(eps,EPS_GHEP);

	const PetscInt max_it = 10000;  // Was not converging with default value of 100.
	const PetscReal tol_eps = 1e-4; // Default value: 1e-8.
	EPSSetTolerances(eps,tol_eps,max_it);
	EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE);
	EPSSetFromOptions(eps);

	EPSSolve(eps);

	PetscInt nconv = 0;
	EPSGetConverged(eps,&nconv);
	if (nconv == 0) {
		PetscInt its = 0;
		EPSGetIterationNumber(eps,&its);
		EXIT_ERROR("Did not converge in %d iterations.\n",its);
	}

	PetscScalar kr = 0.0;
	PetscScalar ki = 0.0;
	Vec xr = NULL;
	Vec xi = NULL;

	const int ind_eig = 0; // Smallest is first.
	EPSGetEigenpair(eps,ind_eig,&kr,&ki,xr,xi);
	assert(ki == 0.0);
	*eig_val_ptr = kr;

	if (PRINT_EIG_VALS_ERRORS) {
		PetscReal error;
		EPSComputeError(eps,ind_eig,EPS_ERROR_RELATIVE,&error);
		PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",kr,error);
	}
	EPSDestroy(&eps);
}

static void set_convergence_order_inf_sup_accept (struct Integration_Test_Info*const int_test_info)
{
	char line[STRLEN_MAX];
	bool found_accept = false;

	FILE* input_file = fopen_input('t',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		if (strstr(line,"conv_order_inf_sup_accept")) {
			found_accept = true;
			read_skip_d(line,&int_test_info->conv_order_inf_sup_accept,1,false);
		}
	}
	fclose(input_file);

	if (!found_accept)
		int_test_info->conv_order_inf_sup_accept = 0.0;
}

static void check_inf_sup
	(bool*const pass, struct Test_Info*const test_info, const struct Integration_Test_Info*const int_test_info,
	 const struct Vector_d*const eig_vals, const struct Simulation*const sim)
{
	UNUSED(test_info); UNUSED(int_test_info); UNUSED(eig_vals); UNUSED(sim);
	*pass = false;

	print_Vector_d(eig_vals);
}
