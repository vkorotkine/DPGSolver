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
 *  \brief Provides the macro definitions used for c-style templating related to the solution functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Solution_Container_T Solution_Container
///\}

///\{ \name Function pointers
#define constructor_sol_fptr_T         constructor_sol_fptr
#define mutable_constructor_sol_fptr_T mutable_constructor_sol_fptr
#define set_sol_fptr_T                 set_sol_fptr
#define compute_source_rhs_fptr_T      compute_source_rhs_fptr
///\}

///\{ \name Function names
#define constructor_const_sol_invalid_T constructor_const_sol_invalid
#define set_initial_solution_T          set_initial_solution
#define set_sg_do_nothing_T             set_sg_do_nothing
#define set_sg_zero_T                   set_sg_zero
#define constructor_xyz_sol_T           constructor_xyz_sol
#define compute_coef_from_val_vs_T      compute_coef_from_val_vs
#define constructor_sol_v_T             constructor_sol_v
#define compute_source_rhs_do_nothing_T compute_source_rhs_do_nothing
#define add_to_flux_imbalance_source_do_nothing_T add_to_flux_imbalance_source_do_nothing
#define update_Solution_Container_sol_T update_Solution_Container_sol
#define update_Solution_Container_grad_T update_Solution_Container_grad
#define constructor_xyz_vc_interp_T     constructor_xyz_vc_interp
#define constructor_Solver_Face__nf_coef_T constructor_Solver_Face__nf_coef
///\}

///\{ \name Function names (pde specific)
#define convert_variables_T convert_variables
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Solution_Container_T Solution_Container_c
///\}

///\{ \name Function pointers
#define constructor_sol_fptr_T         constructor_sol_fptr_c
#define mutable_constructor_sol_fptr_T mutable_constructor_sol_fptr_c
#define set_sol_fptr_T                 set_sol_fptr_c
#define compute_source_rhs_fptr_T      compute_source_rhs_fptr_c
///\}

///\{ \name Function names
#define constructor_const_sol_invalid_T constructor_const_sol_invalid_c
#define set_initial_solution_T          set_initial_solution_c
#define set_sg_do_nothing_T             set_sg_do_nothing_c
#define set_sg_zero_T                   set_sg_zero_c
#define constructor_xyz_sol_T           constructor_xyz_sol_c
#define compute_coef_from_val_vs_T      compute_coef_from_val_vs_c
#define constructor_sol_v_T             constructor_sol_v_c
#define compute_source_rhs_do_nothing_T compute_source_rhs_do_nothing_c
#define add_to_flux_imbalance_source_do_nothing_T add_to_flux_imbalance_source_do_nothing_c
#define update_Solution_Container_sol_T update_Solution_Container_sol_c
#define update_Solution_Container_grad_T update_Solution_Container_grad_c
#define constructor_xyz_vc_interp_T     constructor_xyz_vc_interp_c
#define constructor_Solver_Face__nf_coef_T constructor_Solver_Face__nf_coef_c
///\}

///\{ \name Function names (pde specific)
#define convert_variables_T convert_variables_c
///\}

#endif
