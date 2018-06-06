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
 *  \brief Undefine macro definitions for c-style templated containers/functions relating to solution.
 */

///\{ \name Data types
#undef Solution_Container_T
///\}

///\{ \name Function pointers
#undef constructor_sol_fptr_T
#undef mutable_constructor_sol_fptr_T
#undef set_sol_fptr_T
#undef compute_source_rhs_fptr_T
///\}

///\{ \name Function names
#undef constructor_const_sol_invalid_T
#undef set_initial_solution_T
#undef set_sg_do_nothing_T
#undef set_sg_zero_T
#undef constructor_xyz_sol_T
#undef compute_coef_from_val_vs_T
#undef constructor_sol_v_T
#undef compute_source_rhs_do_nothing_T
#undef add_to_flux_imbalance_source_do_nothing_T
#undef update_Solution_Container_sol_T
#undef update_Solution_Container_grad_T
#undef constructor_xyz_vc_interp_T
#undef constructor_Solver_Face__nf_coef_T
#undef set_grad_zero_T
#undef constructor_const_grad_zero_T
#undef set_to_zero_residual_T
#undef set_initial_v_test_sg_coef_T
///\}

///\{ \name Function names (pde specific)
#undef convert_variables_T
///\}
