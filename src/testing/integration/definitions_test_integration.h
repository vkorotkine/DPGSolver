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

#ifndef DPG__definitions_test_integration_h__INCLUDED
#define DPG__definitions_test_integration_h__INCLUDED
/** \file
 *  \brief Provides the definitions related to the integration testing.
 */

///\{ \name The available options and current value for the linearization components to check.
#define CHECK_LIN_VOLUME 1 ///< Volume contributions only.
#define CHECK_LIN_FACE   2 ///< Face contributions only.
#define CHECK_LIN_ALL    3 ///< All contributions.

#define CHECK_LIN 1 ///< The currently selected value from the options above.
///\}

///\{ \name Magnitude of the complex step to take for linearization testing.
#define CX_STEP 1e-30
///\}

#endif // DPG__definitions_test_integration_h__INCLUDED
