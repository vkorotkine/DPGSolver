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
 *  \brief Provides the macro definitions used for c-style templating related to the dpg solver functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define update_ind_dof_dpg_T  update_ind_dof_dpg
#define constructor_nnz_dpg_T constructor_nnz_dpg
///\}

///\{ \name Static names
#define increment_nnz_off_diag increment_nnz_off_diag
#define increment_nnz_off_diag_constraint increment_nnz_off_diag_constraint
#define increment_nnz_off_diag_v increment_nnz_off_diag_v
#define increment_nnz_off_diag_f increment_nnz_off_diag_f
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define update_ind_dof_dpg_T  update_ind_dof_dpg_c
#define constructor_nnz_dpg_T constructor_nnz_dpg_c
///\}

///\{ \name Static names
#define increment_nnz_off_diag increment_nnz_off_diag_c
#define increment_nnz_off_diag_constraint increment_nnz_off_diag_constraint_c
#define increment_nnz_off_diag_v increment_nnz_off_diag_v_c
#define increment_nnz_off_diag_f increment_nnz_off_diag_f_c
///\}

#endif
