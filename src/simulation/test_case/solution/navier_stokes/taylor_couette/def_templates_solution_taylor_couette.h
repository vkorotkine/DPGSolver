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
 *  \brief Provides the macro definitions used for c-style templating related to the solution functions for the
 *         navier-stokes equations (test case: taylor-couette).
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define set_sol_taylor_couette_T                set_sol_taylor_couette
#define set_grad_taylor_couette_T               set_grad_taylor_couette
#define constructor_const_sol_taylor_couette_T  constructor_const_sol_taylor_couette
#define constructor_const_grad_taylor_couette_T constructor_const_grad_taylor_couette
///\}

///\{ \name Static names
#define constructor_sol_taylor_couette constructor_sol_taylor_couette
#define constructor_grad_taylor_couette constructor_grad_taylor_couette
#define Sol_Data__tc Sol_Data__tc
#define get_sol_data get_sol_data
#define read_data_taylor_couette read_data_taylor_couette
#define set_data_taylor_couette set_data_taylor_couette
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define set_sol_taylor_couette_T                set_sol_taylor_couette_c
#define set_grad_taylor_couette_T               set_grad_taylor_couette_c
#define constructor_const_sol_taylor_couette_T  constructor_const_sol_taylor_couette_c
#define constructor_const_grad_taylor_couette_T constructor_const_grad_taylor_couette_c
///\}

///\{ \name Static names
#define constructor_sol_taylor_couette constructor_sol_taylor_couette_c
#define constructor_grad_taylor_couette constructor_grad_taylor_couette_c
#define Sol_Data__tc Sol_Data__tc_c
#define get_sol_data get_sol_data_c
#define read_data_taylor_couette read_data_taylor_couette_c
#define set_data_taylor_couette set_data_taylor_couette_c
///\}

#endif
