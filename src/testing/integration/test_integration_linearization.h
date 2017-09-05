// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_integration_linearization_h__INCLUDED
#define DPG__test_integration_linearization_h__INCLUDED
/**	\file
 *	\brief Provides functionality for integration testing of the solver linearization.
 */

/** \brief Container for the information relating to the linearization testing.
 *	\todo Delete if unused.
 */
struct Test_Int_Linearization {
	const char *const ctrl_name; ///< Name (with relative path) of the control file for the test.
};

/**	\test Performs integration testing for the solver linearization.
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
void test_integration_linearization
	(const char*const ctrl_name ///< Defined in \ref Test_Int_Linearization.
	);


#endif // DPG__test_integration_linearization_h__INCLUDED
