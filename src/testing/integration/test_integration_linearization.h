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

#ifndef DPG__test_integration_linearization_h__INCLUDED
#define DPG__test_integration_linearization_h__INCLUDED
/** \file
 *  \brief Provides functionality for integration testing of the solver linearization.
 */

struct Test_Info;

/** \brief Container for the information relating to the linearization testing.
 *  \todo Delete if unused.
 */
struct Test_Int_Linearization {
	const char *const ctrl_name; ///< Name (with relative path) of the control file for the test.
};

/** \test Performs integration testing for the solver linearization.
 *  \return The expected output is the correspondence of LHS matrices computed using complex step and exact
 *          linearization indicated by a passing test.
 *
 *  The linearization is verified by comparing with the output when using the complex step method.
 *
 *  For second order equations, the verification of the linearization of the weak gradients is also performed. Further,
 *  the volume rhs contributes off-diagonal terms to the global system matrix due to the use of the fully corrected
 *  gradient. A flag is provided to avoid the computation of these terms when checking only the diagonal volume/face
 *  contributions to the global system using compute_A_cs with \c assembly_type equal to \c 1 or \c 2.
 *
 * Details of the complex step method can be found in Squire et al. \cite Squire1998 and Martins et al.
 * \cite Martins2003.
 *
 * By default, the complete linearization is checked here. However, linearizations of the individual contributions
 * listed below may be checked separately:
 *	- diagonal volume contributions
 *	- diagonal face contributions
 *	- off-diagonal face (and volume for 2nd order equations) contributions
 *
 * Further, the complete linearization can be checked either through the assembly of the individual contributions or
 * more simply by using the assembled RHS directly.
 *
 * \todo Update comments if necessary after HDG is implemented.
 */
void test_integration_linearization
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const char*const ctrl_name        ///< The name of the input control file.
	);


#endif // DPG__test_integration_linearization_h__INCLUDED
