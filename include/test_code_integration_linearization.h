// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_code_integration_linearization_h__INCLUDED
#define DPG__test_code_integration_linearization_h__INCLUDED
/// \file

#include <stdbool.h>

#include "petscmat.h"

/// \brief Container for information relating to linearization testing.
struct Test_Linearization {
	const int nargc; ///< Standard.
	char **argv_new, ///< New argv to be used for the test case under consideration.
	     *test_name; ///< Name of the test.

	bool         print_mat_to_file,        ///< Flag for whether the petsc `Mat` is printed to file.
	             print_timings,            /**< Flag for whether print timing comparison of complex step vs. analytical
	                                            linearization is done. */
	             check_full_linearization, ///< Flag for whether the full linearization is checked.
	             check_weak_gradients,     ///< Flag for whether the weak gradient linearization is checked.
	             omit_root;                ///< Flag for whether the root of the test name should be omitted.
	unsigned int n_ref,                    ///< Defined in \ref test_code_integration.h
	             modify_params,            ///< Defined in \ref test_code_integration.h
	             PGlobal,                  ///< Polynomial order to be used for the test.
	             ML,
	             PG_add,
	             IntOrder_mult,
	             IntOrder_add,
	             check_level;

	Mat A, A_cs, A_csc;
};

/** \test Performs integration testing for the solver linearization.
 *
 *	The linearization is verified by comparing with the output when using the complex step method.
 *
 *	For second order equations, the verification of the linearization of the weak gradients is also performed. Further,
 *	the volume rhs contributes off-diagonal terms to the global system matrix due to the use of the fully corrected
 *	gradient. A flag is provided to avoid the computation of these terms when checking only the diagonal volume/face
 *	contributions to the global system using compute_A_cs with \c assembly_type equal to \c 1 or \c 2.
 *
 *	Details of the complex step method can be found in Squire et al. \cite Squire1998 and Martins et al.
 *	\cite Martins2003.
 *
 *	By default, the complete linearization is checked here. However, linearizations of the individual contributions
 *	listed below may be checked separately:
 *		- diagonal volume contributions
 *		- diagonal face contributions
 *		- off-diagonal face (and volume for 2nd order equations) contributions
 *
 *	Further, the complete linearization can be checked either through the assembly of the individual contributions or
 *	more simply by using the assembled RHS directly.
 *
 *	\return The expected output is the correspondence of LHS matrices computed using complex step and exact
 *	        linearization indicated by a passing test.
 *
 *	\todo Update comments if necessary after HDG is implemented.
 */
void test_linearization
	(struct Test_Linearization *const data, ///< The data specific to the instance of the test.
	 char const *const TestName);           ///< The test name.

#endif // DPG__test_code_integration_linearization_h__INCLUDED
