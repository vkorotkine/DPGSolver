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
 *         euler equations (test case: free stream).
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define set_sol_free_stream_T               set_sol_free_stream
#define constructor_const_sol_free_stream_T constructor_const_sol_free_stream
///\}

///\{ \name Static names
#define constructor_sol_free_stream constructor_sol_free_stream
#define mutable_constructor_sol_fs_fptr_T mutable_constructor_sol_fs_fptr_T
#define Sol_Data__fs Sol_Data__fs
#define get_sol_data get_sol_data
#define read_data_free_stream read_data_free_stream
#define constructor_sol_free_stream_const constructor_sol_free_stream_const
#define constructor_sol_free_stream_trig_x constructor_sol_free_stream_trig_x
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define set_sol_free_stream_T               set_sol_free_stream_c
#define constructor_const_sol_free_stream_T constructor_const_sol_free_stream_c
///\}

///\{ \name Static names
#define constructor_sol_free_stream constructor_sol_free_stream_c
#define mutable_constructor_sol_fs_fptr_T mutable_constructor_sol_fs_fptr_T_c
#define Sol_Data__fs Sol_Data__fs_c
#define get_sol_data get_sol_data_c
#define read_data_free_stream read_data_free_stream_c
#define constructor_sol_free_stream_const constructor_sol_free_stream_const_c
#define constructor_sol_free_stream_trig_x constructor_sol_free_stream_trig_x_c
///\}

#endif
