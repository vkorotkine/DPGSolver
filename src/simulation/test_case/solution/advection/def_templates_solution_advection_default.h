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
 *         linear advection equation (test case: default).
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define set_sol_advection_default_T               set_sol_advection_default
#define constructor_const_sol_advection_default_T constructor_const_sol_advection_default
#define compute_source_rhs_advection_default_T    compute_source_rhs_advection_default
#define add_to_flux_imbalance_source_advection_default_T add_to_flux_imbalance_source_advection_default
///\}

#define constructor_sol_advection_default_T constructor_sol_advection_default_T
#define constructor_source_advection_default_T constructor_source_advection_default_T
#define constructor_sol_advection_default_1d_T constructor_sol_advection_default_1d_T
#define constructor_sol_advection_default_2d_T constructor_sol_advection_default_2d_T
#define constructor_source_advection_default_1d_T constructor_source_advection_default_1d_T
#define constructor_source_advection_default_2d_T constructor_source_advection_default_2d_T

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define set_sol_advection_default_T               set_sol_advection_default_c
#define constructor_const_sol_advection_default_T constructor_const_sol_advection_default_c
#define compute_source_rhs_advection_default_T    compute_source_rhs_advection_default_c
#define add_to_flux_imbalance_source_advection_default_T add_to_flux_imbalance_source_advection_default_c
///\}

#define constructor_sol_advection_default_T constructor_sol_advection_default_T_c
#define constructor_source_advection_default_T constructor_source_advection_default_T_c
#define constructor_sol_advection_default_1d_T constructor_sol_advection_default_1d_T_c
#define constructor_sol_advection_default_2d_T constructor_sol_advection_default_2d_T_c
#define constructor_source_advection_default_1d_T constructor_source_advection_default_1d_T_c
#define constructor_source_advection_default_2d_T constructor_source_advection_default_2d_T_c

#endif
