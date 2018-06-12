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
 *         euler equations (test case: periodic vortex).
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define set_sol_periodic_vortex_T set_sol_periodic_vortex
///\}

#define constructor_sol_periodic_vortex_T constructor_sol_periodic_vortex_T
#define Sol_Data__pv Sol_Data__pv
#define get_sol_data get_sol_data
#define set_xy_c set_xy_c
#define read_data_periodic_vortex read_data_periodic_vortex
#define set_data_periodic_vortex set_data_periodic_vortex

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define set_sol_periodic_vortex_T set_sol_periodic_vortex_c
///\}

#define constructor_sol_periodic_vortex_T constructor_sol_periodic_vortex_T_c
#define Sol_Data__pv Sol_Data__pv_c
#define get_sol_data get_sol_data_c
#define set_xy_c set_xy_c_c
#define read_data_periodic_vortex read_data_periodic_vortex_c
#define set_data_periodic_vortex set_data_periodic_vortex_c

#endif
