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
 *  \brief Provides the macro definitions used for c-style templating related to the `double complex` multiarray constructor
 *         functions.
 */

///\{ \name Function names
#define constructor_empty_Multiarray_T             constructor_empty_Multiarray_c             
#define constructor_empty_Multiarray_T_dyn_extents constructor_empty_Multiarray_c_dyn_extents 

#define constructor_copy_Multiarray_T       constructor_copy_Multiarray_c       
#define constructor_copy_const_Multiarray_T constructor_copy_const_Multiarray_c 

#define constructor_move_Multiarray_T_T           constructor_move_Multiarray_c_c           
#define constructor_move_const_Multiarray_T_T     constructor_move_const_Multiarray_c_c     
#define constructor_move_Multiarray_T_dyn_extents constructor_move_Multiarray_c_dyn_extents 

#define destructor_Multiarray_T       destructor_Multiarray_c       
#define destructor_const_Multiarray_T destructor_const_Multiarray_c 
///\}
