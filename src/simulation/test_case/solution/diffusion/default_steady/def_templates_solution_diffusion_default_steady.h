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
 *         diffusion equation (test case: default_steady).
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define set_sol_diffusion_default_steady_T                      set_sol_diffusion_default_steady
#define set_grad_diffusion_default_steady_T                     set_grad_diffusion_default_steady
#define constructor_const_sol_diffusion_default_steady_T        constructor_const_sol_diffusion_default_steady
#define constructor_const_grad_diffusion_default_steady_T       constructor_const_grad_diffusion_default_steady
#define compute_source_rhs_diffusion_default_steady_T           compute_source_rhs_diffusion_default_steady
#define add_to_flux_imbalance_source_diffusion_default_steady_T add_to_flux_imbalance_source_diffusion_default_steady
///\}

///\{ \name Static names
#define constructor_sol_diffusion_default_steady constructor_sol_diffusion_default_steady
#define constructor_grad_diffusion_default_steady constructor_grad_diffusion_default_steady
#define constructor_source_diffusion_default_steady constructor_source_diffusion_default_steady
#define constructor_sol_diffusion_default_steady_1d constructor_sol_diffusion_default_steady_1d
#define constructor_grad_diffusion_default_steady_1d constructor_grad_diffusion_default_steady_1d
#define constructor_source_diffusion_default_steady_1d constructor_source_diffusion_default_steady_1d
#define constructor_sol_diffusion_default_steady_2d constructor_sol_diffusion_default_steady_2d
#define constructor_grad_diffusion_default_steady_2d constructor_grad_diffusion_default_steady_2d
#define constructor_source_diffusion_default_steady_2d constructor_source_diffusion_default_steady_2d
#define Sol_Data__dd Sol_Data__dd
#define get_sol_data get_sol_data
#define read_data_default_diffusion read_data_default_diffusion
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define set_sol_diffusion_default_steady_T                      set_sol_diffusion_default_steady_c
#define set_grad_diffusion_default_steady_T                     set_grad_diffusion_default_steady_c
#define constructor_const_sol_diffusion_default_steady_T        constructor_const_sol_diffusion_default_steady_c
#define constructor_const_grad_diffusion_default_steady_T       constructor_const_grad_diffusion_default_steady_c
#define compute_source_rhs_diffusion_default_steady_T           compute_source_rhs_diffusion_default_steady_c
#define add_to_flux_imbalance_source_diffusion_default_steady_T add_to_flux_imbalance_source_diffusion_default_steady_c
///\}

///\{ \name Static names
#define constructor_sol_diffusion_default_steady constructor_sol_diffusion_default_steady_c
#define constructor_grad_diffusion_default_steady constructor_grad_diffusion_default_steady_c
#define constructor_source_diffusion_default_steady constructor_source_diffusion_default_steady_c
#define constructor_sol_diffusion_default_steady_1d constructor_sol_diffusion_default_steady_1d_c
#define constructor_grad_diffusion_default_steady_1d constructor_grad_diffusion_default_steady_1d_c
#define constructor_source_diffusion_default_steady_1d constructor_source_diffusion_default_steady_1d_c
#define constructor_sol_diffusion_default_steady_2d constructor_sol_diffusion_default_steady_2d_c
#define constructor_grad_diffusion_default_steady_2d constructor_grad_diffusion_default_steady_2d_c
#define constructor_source_diffusion_default_steady_2d constructor_source_diffusion_default_steady_2d_c
#define Sol_Data__dd Sol_Data__dd_c
#define get_sol_data get_sol_data_c
#define read_data_default_diffusion read_data_default_diffusion_c
///\}

#endif
