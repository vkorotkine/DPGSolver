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

#ifndef DPG__definitions_test_case_h__INCLUDED
#define DPG__definitions_test_case_h__INCLUDED
/** \file
 *  \brief Provides the definitions relating to test cases.
 */

///\{ \name Definitions for the available PDEs.
#define PDE_ADVECTION     1
#define PDE_POISSON       2
#define PDE_EULER         3
#define PDE_NAVIER_STOKES 4
///\}

///\{ \name Definitions for variables related to the PDEs.
#define GAMMA  1.4
#define GM1    0.4
#define GM3   -1.6
///\}

#endif // DPG__definitions_test_case_h__INCLUDED
