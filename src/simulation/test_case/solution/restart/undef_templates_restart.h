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
 *  \brief Undefine macro definitions for c-style templated containers/functions relating to restart.
 */

///\{ \name Data types
//#undef Solution_Container_T
///\}

///\{ \name Function names
#undef set_sol_restart_T
#undef constructor_const_sol_restart_T
///\}

#undef constructor_sol_restart
#undef Restart_Info
#undef get_Restart_Info
#undef constructor_indices_vol_background
#undef constructor_sol_restart_from_background
#undef constructor_volume_list
#undef initialize_volumes_restart
#undef initialize_ann_background
#undef update_ind_vol_b
#undef initialize_volumes_sol_coef
#undef get_xyz_node_m_ve
#undef read_sol_coef_bezier
#undef get_ind_ve
