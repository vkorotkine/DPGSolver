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
#define PDE_DIFFUSION     2
#define PDE_EULER         3
#define PDE_NAVIER_STOKES 4
///\}

///\{ \name Definitions for the available solver methods.
#define METHOD_DG   1
#define METHOD_HDG  2
#define METHOD_HDPG 3
#define METHOD_DPG  4
///\}

///\{ \name Definitions for the available solver procedures.
#define SOLVER_E  100 ///< Explicit.
#define SOLVER_I  200 ///< Implicit.
#define SOLVER_EI 300 ///< Explicit then implicit.
///\}

///\{ \name Definitions for the available solver types.
#define SOLVER_E_EULER     101 ///< Explicit forward Euler.
#define SOLVER_E_SSP_RK_33 102 ///< Explicit strong stability preserving Runge-Kutta (3-stage, 3rd order).
#define SOLVER_E_LS_RK_54  103 ///< Explicit low storage Runge-Kutta (5-stage, 4rd order).

#define SOLVER_I_DIRECT    201 ///< Implicit direct solver (LU; Cholesky if symmetric).
#define SOLVER_I_ITERATIVE 202 ///< Implicit iterative solver (specification provided in input Petsc options file).
///\}

#endif // DPG__definitions_test_case_h__INCLUDED
