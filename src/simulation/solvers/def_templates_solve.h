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
 *  \brief Provides the macro definitions used for c-style templating related to the general solve functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define compute_dof_T                         compute_dof
#define update_ind_dof_T                      update_ind_dof
#define constructor_Solver_Storage_Implicit_T constructor_Solver_Storage_Implicit
#define add_to_flux_imbalance_source_T        add_to_flux_imbalance_source
#define get_operator__tw0_vt_vc_T             get_operator__tw0_vt_vc
#define initialize_zero_memory_volumes_T      initialize_zero_memory_volumes
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define compute_dof_T                         compute_dof_c
#define update_ind_dof_T                      update_ind_dof_c
#define constructor_Solver_Storage_Implicit_T constructor_Solver_Storage_Implicit_c
#define add_to_flux_imbalance_source_T        add_to_flux_imbalance_source_c
#define get_operator__tw0_vt_vc_T             get_operator__tw0_vt_vc_c
#define initialize_zero_memory_volumes_T      initialize_zero_memory_volumes_c
///\}

#endif
