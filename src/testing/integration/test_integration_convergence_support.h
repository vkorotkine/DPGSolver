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

#ifndef DPG__test_integration_convergence_support_h__INCLUDED
#define DPG__test_integration_convergence_support_h__INCLUDED
/** \file
 *  \brief Provides support functions/structures for convergence order integration testing.
 */

///\{ Options for the type of "solution" procedure to use for the convergence study.
#define CONV_STUDY_SOLVE   101 ///< Solve for the solution in the standard manner.
#define CONV_STUDY_RESTART 102 ///< Use the restarted solution to check convergence orders.
///\}

/// \brief Perform a convergence order study of the specified type based on the control file parameters.
void run_convergence_order_study
	(int argc,                 ///< Standard.
	 char** argv,              ///< Standard.
	 const int conv_study_type ///< The type of "solution" procedure to use for the study. Options: see above.
	);

#endif // DPG__test_integration_convergence_support_h__INCLUDED
