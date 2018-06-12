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
 *  \brief Provides the macro definitions used for c-style templating related to the restart containers/functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
//#define Solution_Container_T Solution_Container
///\}

///\{ \name Function names
#define set_sol_restart_T               set_sol_restart_d
#define constructor_const_sol_restart_T constructor_const_sol_restart_d
///\}

#define constructor_sol_restart constructor_sol_restart
#define Restart_Info Restart_Info
#define get_Restart_Info get_Restart_Info
#define constructor_indices_vol_background constructor_indices_vol_background
#define constructor_sol_restart_from_background constructor_sol_restart_from_background
#define constructor_volume_list constructor_volume_list
#define initialize_volumes_restart initialize_volumes_restart
#define initialize_ann_background initialize_ann_background
#define update_ind_vol_b update_ind_vol_b
#define initialize_volumes_sol_coef initialize_volumes_sol_coef
#define get_xyz_node_m_ve get_xyz_node_m_ve
#define read_sol_coef_bezier read_sol_coef_bezier
#define get_ind_ve get_ind_ve

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
//#define Solution_Container_T Solution_Container_c
///\}

///\{ \name Function names
#define set_sol_restart_T               set_sol_restart_c
#define constructor_const_sol_restart_T constructor_const_sol_restart_c
///\}

#define constructor_sol_restart constructor_sol_restart_c
#define Restart_Info Restart_Info_c
#define get_Restart_Info get_Restart_Info_c
#define constructor_indices_vol_background constructor_indices_vol_background_c
#define constructor_sol_restart_from_background constructor_sol_restart_from_background_c
#define constructor_volume_list constructor_volume_list_c
#define initialize_volumes_restart initialize_volumes_restart_c
#define initialize_ann_background initialize_ann_background_c
#define update_ind_vol_b update_ind_vol_b_c
#define initialize_volumes_sol_coef initialize_volumes_sol_coef_c
#define get_xyz_node_m_ve get_xyz_node_m_ve_c
#define read_sol_coef_bezier read_sol_coef_bezier_c
#define get_ind_ve get_ind_ve_c

#endif
