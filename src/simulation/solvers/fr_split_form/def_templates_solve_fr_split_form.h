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
 *  \brief Provides the macro definitions used for c-style templating related to the fr_split_form solver functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define update_ind_dof_fr_split_form_T     update_ind_dof_fr_split_form
#define constructor_nnz_frsf_T    constructor_nnz_fr_split_form
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define update_ind_dof_fr_split_form_T     update_ind_dof_fr_split_form_c
#define constructor_nnz_frsf_T    constructor_nnz_fr_split_form_c
///\}

#endif
