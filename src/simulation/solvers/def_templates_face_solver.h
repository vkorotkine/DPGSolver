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
 *  \brief Provides the macro definitions used for c-style templating related to the \ref Solver_Face_T
 *         containers/functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Solver_Face_T Solver_Face
///\}

///\{ \name Function names
#define constructor_derived_Solver_Face_T     constructor_derived_Solver_Face
#define destructor_derived_Solver_Face_T      destructor_derived_Solver_Face
#define set_function_pointers_face_num_flux_T set_function_pointers_face_num_flux
#define get_operator__w_fc__s_e_T             get_operator__w_fc__s_e
#define constructor_mass_face_T constructor_mass_face_d
///\}

///\{ \name Static names
#define set_function_pointers_num_flux_bc set_function_pointers_num_flux_bc
#define set_function_pointers_num_flux_bc_advection set_function_pointers_num_flux_bc_advection
#define set_function_pointers_num_flux_bc_diffusion set_function_pointers_num_flux_bc_diffusion
#define set_function_pointers_num_flux_bc_euler set_function_pointers_num_flux_bc_euler
#define set_function_pointers_num_flux_bc_navier_stokes set_function_pointers_num_flux_bc_navier_stokes
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Solver_Face_T Solver_Face_c
///\}

///\{ \name Function names
#define constructor_derived_Solver_Face_T     constructor_derived_Solver_Face_c
#define destructor_derived_Solver_Face_T      destructor_derived_Solver_Face_c
#define set_function_pointers_face_num_flux_T set_function_pointers_face_num_flux_c
#define get_operator__w_fc__s_e_T             get_operator__w_fc__s_e_c
#define constructor_mass_face_T constructor_mass_face_c
///\}

///\{ \name Static names
#define set_function_pointers_num_flux_bc set_function_pointers_num_flux_bc_c
#define set_function_pointers_num_flux_bc_advection set_function_pointers_num_flux_bc_advection_c
#define set_function_pointers_num_flux_bc_diffusion set_function_pointers_num_flux_bc_diffusion_c
#define set_function_pointers_num_flux_bc_euler set_function_pointers_num_flux_bc_euler_c
#define set_function_pointers_num_flux_bc_navier_stokes set_function_pointers_num_flux_bc_navier_stokes_c
///\}

#endif
