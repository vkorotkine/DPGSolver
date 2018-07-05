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
 *  \brief Provides the macro definitions used for c-style templating related to the multiarray constructor functions.
 */

#if defined TYPE_RC
#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define constructor_default_Multiarray_T              constructor_default_Multiarray_d
#define constructor_default_const_Multiarray_T        constructor_default_const_Multiarray_d
#define constructor_default_Multiarray_Matrix_T       constructor_default_Multiarray_Matrix_d
#define constructor_default_const_Multiarray_Matrix_T constructor_default_const_Multiarray_Matrix_d

#define constructor_empty_Multiarray_T                constructor_empty_Multiarray_d
#define constructor_empty_const_Multiarray_T          constructor_empty_const_Multiarray_d
#define constructor_empty_Multiarray_T_dyn_extents    constructor_empty_Multiarray_d_dyn_extents
#define constructor_empty_Multiarray_Vector_T         constructor_empty_Multiarray_Vector_d
#define constructor_empty_const_Multiarray_Vector_T   constructor_empty_const_Multiarray_Vector_d
#define constructor_empty_const_Multiarray_Vector_T_V constructor_empty_const_Multiarray_Vector_d_V
#define constructor_empty_Multiarray_Matrix_T         constructor_empty_Multiarray_Matrix_d
#define constructor_empty_const_Multiarray_Matrix_T   constructor_empty_const_Multiarray_Matrix_d
#define constructor_empty_Multiarray_Matrix_T_V       constructor_empty_Multiarray_Matrix_d_V
#define constructor_empty_const_Multiarray_Matrix_T_V constructor_empty_const_Multiarray_Matrix_d_V

#define constructor_zero_Multiarray_T             constructor_zero_Multiarray_d
#define constructor_zero_Multiarray_T_dyn_extents constructor_zero_Multiarray_d_dyn_extents

#define constructor_copy_Multiarray_T          constructor_copy_Multiarray_d
#define constructor_copy_const_Multiarray_T    constructor_copy_const_Multiarray_d
#define constructor_copy_Multiarray_Vector_T_T constructor_copy_Multiarray_Vector_d_d
#define constructor_copy_Multiarray_Vector_T   constructor_copy_Multiarray_Vector_d
#define const_constructor_copy_Multiarray_T    const_constructor_copy_Multiarray_d
#define constructor_copy_Multiarray_T_Multiarray_R       constructor_copy_Multiarray_d_Multiarray_d
#define constructor_copy_const_Multiarray_T_Multiarray_R constructor_copy_const_Multiarray_d_Multiarray_d

#define constructor_move_Multiarray_T_T                  constructor_move_Multiarray_d_d
#define constructor_move_const_Multiarray_T_T            constructor_move_const_Multiarray_d_d
#define constructor_move_Multiarray_T_dyn_extents        constructor_move_Multiarray_d_dyn_extents
#define constructor_move_const_Multiarray_T_dyn_extents  constructor_move_const_Multiarray_d_dyn_extents
#define constructor_move_Multiarray_Vector_T_dyn_extents constructor_move_Multiarray_Vector_d_dyn_extents
#define constructor_move_Multiarray_Matrix_T_dyn_extents constructor_move_Multiarray_Matrix_d_dyn_extents
#define constructor_move_Multiarray_T_Matrix_T           constructor_move_Multiarray_d_Matrix_d
#define constructor_move_const_Multiarray_T_Matrix_T     constructor_move_const_Multiarray_d_Matrix_d
#define const_constructor_move_Multiarray_T              const_constructor_move_Multiarray_d
#define const_constructor_move_const_Multiarray_T        const_constructor_move_const_Multiarray_d
#define const_constructor_move_Multiarray_Vector_T       const_constructor_move_Multiarray_Vector_d
#define const_constructor_move_Multiarray_Matrix_T       const_constructor_move_Multiarray_Matrix_d

#define constructor_sum_Multiarrays_Multiarray_T             constructor_sum_Multiarrays_Multiarray_d
#define constructor_sum_Multiarrays_const_Multiarray_T       constructor_sum_Multiarrays_const_Multiarray_d
#define constructor_MaM1_V_const_Multiarray_T                constructor_MaM1_V_const_Multiarray_d
#define set_Multiarray_Matrix_from_Multiarray_Matrix_T       set_Multiarray_Matrix_from_Multiarray_Matrix_d
#define set_const_Multiarray_Matrix_from_Multiarray_Matrix_T set_const_Multiarray_Matrix_from_Multiarray_Matrix_d
#define constructor_mm_NN1C_Multiarray_T                     constructor_mm_NN1C_Multiarray_d
#define constructor_mm_NN1C_const_Multiarray_T               constructor_mm_NN1C_const_Multiarray_d
#define constructor_mm_NN1C_Multiarray_TT constructor_mm_NN1C_Multiarray_dd
#define constructor_mm_tp_NN1C_const_Multiarray_T            constructor_mm_tp_NN1C_const_Multiarray_d

#define destructor_Multiarray_T                   destructor_Multiarray_d
#define destructor_const_Multiarray_T             destructor_const_Multiarray_d
#define destructor_conditional_Multiarray_T       destructor_conditional_Multiarray_d
#define destructor_conditional_const_Multiarray_T destructor_conditional_const_Multiarray_d
#define destructor_Multiarray_Vector_T            destructor_Multiarray_Vector_d
#define destructor_const_Multiarray_Vector_T      destructor_const_Multiarray_Vector_d
#define destructor_conditional_Multiarray_Vector_T       destructor_conditional_Multiarray_Vector_d
#define destructor_conditional_const_Multiarray_Vector_T destructor_conditional_const_Multiarray_Vector_d
#define destructor_Multiarray_Matrix_T            destructor_Multiarray_Matrix_d
#define destructor_const_Multiarray_Matrix_T      destructor_const_Multiarray_Matrix_d
#define destructor_const_Multiarray2_Matrix_T     destructor_const_Multiarray2_Matrix_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define constructor_default_Multiarray_T              constructor_default_Multiarray_c
#define constructor_default_const_Multiarray_T        constructor_default_const_Multiarray_c
#define constructor_default_Multiarray_Matrix_T       constructor_default_Multiarray_Matrix_c
#define constructor_default_const_Multiarray_Matrix_T constructor_default_const_Multiarray_Matrix_c

#define constructor_empty_Multiarray_T                constructor_empty_Multiarray_c
#define constructor_empty_const_Multiarray_T          constructor_empty_const_Multiarray_c
#define constructor_empty_Multiarray_T_dyn_extents    constructor_empty_Multiarray_c_dyn_extents
#define constructor_empty_Multiarray_Vector_T         constructor_empty_Multiarray_Vector_c
#define constructor_empty_const_Multiarray_Vector_T   constructor_empty_const_Multiarray_Vector_c
#define constructor_empty_const_Multiarray_Vector_T_V constructor_empty_const_Multiarray_Vector_c_V
#define constructor_empty_Multiarray_Matrix_T         constructor_empty_Multiarray_Matrix_c
#define constructor_empty_const_Multiarray_Matrix_T   constructor_empty_const_Multiarray_Matrix_c
#define constructor_empty_Multiarray_Matrix_T_V       constructor_empty_Multiarray_Matrix_c_V
#define constructor_empty_const_Multiarray_Matrix_T_V constructor_empty_const_Multiarray_Matrix_c_V

#define constructor_zero_Multiarray_T             constructor_zero_Multiarray_c
#define constructor_zero_Multiarray_T_dyn_extents constructor_zero_Multiarray_c_dyn_extents

#define constructor_copy_Multiarray_T                    constructor_copy_Multiarray_c
#define constructor_copy_const_Multiarray_T              constructor_copy_const_Multiarray_c
#define constructor_copy_Multiarray_Vector_T_T           constructor_copy_Multiarray_Vector_c_c
#define constructor_copy_Multiarray_Vector_T   constructor_copy_Multiarray_Vector_c
#define const_constructor_copy_Multiarray_T              const_constructor_copy_Multiarray_c
#define constructor_copy_Multiarray_T_Multiarray_R       constructor_copy_Multiarray_c_Multiarray_d
#define constructor_copy_const_Multiarray_T_Multiarray_R constructor_copy_const_Multiarray_c_Multiarray_d

#define constructor_move_Multiarray_T_T                  constructor_move_Multiarray_c_c
#define constructor_move_const_Multiarray_T_T            constructor_move_const_Multiarray_c_c
#define constructor_move_Multiarray_T_dyn_extents        constructor_move_Multiarray_c_dyn_extents
#define constructor_move_const_Multiarray_T_dyn_extents  constructor_move_const_Multiarray_c_dyn_extents
#define constructor_move_Multiarray_Vector_T_dyn_extents constructor_move_Multiarray_Vector_c_dyn_extents
#define constructor_move_Multiarray_Matrix_T_dyn_extents constructor_move_Multiarray_Matrix_c_dyn_extents
#define constructor_move_Multiarray_T_Matrix_T           constructor_move_Multiarray_c_Matrix_c
#define constructor_move_const_Multiarray_T_Matrix_T     constructor_move_const_Multiarray_c_Matrix_c
#define const_constructor_move_Multiarray_T              const_constructor_move_Multiarray_c
#define const_constructor_move_const_Multiarray_T        const_constructor_move_const_Multiarray_c
#define const_constructor_move_Multiarray_Vector_T       const_constructor_move_Multiarray_Vector_c
#define const_constructor_move_Multiarray_Matrix_T       const_constructor_move_Multiarray_Matrix_c

#define constructor_sum_Multiarrays_Multiarray_T             constructor_sum_Multiarrays_Multiarray_c
#define constructor_sum_Multiarrays_const_Multiarray_T       constructor_sum_Multiarrays_const_Multiarray_c
#define constructor_MaM1_V_const_Multiarray_T                constructor_MaM1_V_const_Multiarray_c
#define set_Multiarray_Matrix_from_Multiarray_Matrix_T       set_Multiarray_Matrix_from_Multiarray_Matrix_c
#define set_const_Multiarray_Matrix_from_Multiarray_Matrix_T set_const_Multiarray_Matrix_from_Multiarray_Matrix_c
#define constructor_mm_NN1C_Multiarray_T                     constructor_mm_NN1C_Multiarray_c
#define constructor_mm_NN1C_const_Multiarray_T               constructor_mm_NN1C_const_Multiarray_c
#define constructor_mm_NN1C_Multiarray_TT constructor_mm_NN1C_Multiarray_cc
#define constructor_mm_tp_NN1C_const_Multiarray_T            constructor_mm_tp_NN1C_const_Multiarray_c

#define destructor_Multiarray_T                   destructor_Multiarray_c
#define destructor_const_Multiarray_T             destructor_const_Multiarray_c
#define destructor_conditional_Multiarray_T       destructor_conditional_Multiarray_c
#define destructor_conditional_const_Multiarray_T destructor_conditional_const_Multiarray_c
#define destructor_Multiarray_Vector_T            destructor_Multiarray_Vector_c
#define destructor_const_Multiarray_Vector_T      destructor_const_Multiarray_Vector_c
#define destructor_conditional_Multiarray_Vector_T       destructor_conditional_Multiarray_Vector_c
#define destructor_conditional_const_Multiarray_Vector_T destructor_conditional_const_Multiarray_Vector_c
#define destructor_Multiarray_Matrix_T            destructor_Multiarray_Matrix_c
#define destructor_const_Multiarray_Matrix_T      destructor_const_Multiarray_Matrix_c
#define destructor_const_Multiarray2_Matrix_T     destructor_const_Multiarray2_Matrix_c
///\}

#endif

///\{ \name Function names
#define constructor_default_Multiarray_R              constructor_default_Multiarray_d

#define constructor_empty_Multiarray_R       constructor_empty_Multiarray_d
#define constructor_empty_const_Multiarray_R constructor_empty_const_Multiarray_d

#define constructor_copy_Multiarray_R          constructor_copy_Multiarray_d
#define constructor_copy_const_Multiarray_R    constructor_copy_const_Multiarray_d
#define const_constructor_copy_Multiarray_R    const_constructor_copy_Multiarray_d

#define const_constructor_move_Multiarray_R              const_constructor_move_Multiarray_d
#define const_constructor_move_const_Multiarray_R        const_constructor_move_const_Multiarray_d

#define destructor_Multiarray_R                   destructor_Multiarray_d
#define destructor_const_Multiarray_R             destructor_const_Multiarray_d
#define destructor_conditional_Multiarray_R       destructor_conditional_Multiarray_d
#define destructor_conditional_const_Multiarray_R destructor_conditional_const_Multiarray_d
///\}


#elif defined TYPE_I

#if TYPE_I == TYPE_II

///\{ \name Function names
#define constructor_default_Multiarray_T              constructor_default_Multiarray_i
#define constructor_default_const_Multiarray_T        constructor_default_const_Multiarray_i
#define constructor_default_Multiarray_Matrix_T       constructor_default_Multiarray_Matrix_i
#define constructor_default_const_Multiarray_Matrix_T constructor_default_const_Multiarray_Matrix_i

#define constructor_empty_Multiarray_T                constructor_empty_Multiarray_i
#define constructor_empty_Multiarray_T_dyn_extents    constructor_empty_Multiarray_i_dyn_extents
#define constructor_empty_Multiarray_Vector_T         constructor_empty_Multiarray_Vector_i
#define constructor_empty_const_Multiarray_Vector_T   constructor_empty_const_Multiarray_Vector_i
#define constructor_empty_const_Multiarray_Vector_T_V constructor_empty_const_Multiarray_Vector_i_V
#define constructor_empty_Multiarray_Vector_T         constructor_empty_Multiarray_Vector_i
#define constructor_empty_const_Multiarray_Vector_T   constructor_empty_const_Multiarray_Vector_i
#define constructor_empty_const_Multiarray_Vector_T_V constructor_empty_const_Multiarray_Vector_i_V
#define constructor_empty_Multiarray_Matrix_T         constructor_empty_Multiarray_Matrix_i
#define constructor_empty_const_Multiarray_Matrix_T   constructor_empty_const_Multiarray_Matrix_i
#define constructor_empty_Multiarray_Matrix_T_V       constructor_empty_Multiarray_Matrix_i_V
#define constructor_empty_const_Multiarray_Matrix_T_V constructor_empty_const_Multiarray_Matrix_i_V

#define constructor_zero_Multiarray_T             constructor_zero_Multiarray_i
#define constructor_zero_Multiarray_T_dyn_extents constructor_zero_Multiarray_i_dyn_extents

#define constructor_copy_Multiarray_T          constructor_copy_Multiarray_i
#define constructor_copy_const_Multiarray_T    constructor_copy_const_Multiarray_i
#define constructor_copy_Multiarray_Vector_T_T constructor_copy_Multiarray_Vector_i_i
#define constructor_copy_Multiarray_Vector_T   constructor_copy_Multiarray_Vector_i
#define const_constructor_copy_Multiarray_T    const_constructor_copy_Multiarray_i

#define constructor_move_Multiarray_T_T                  constructor_move_Multiarray_i_i
#define constructor_move_const_Multiarray_T_T            constructor_move_const_Multiarray_i_i
#define constructor_move_Multiarray_T_dyn_extents        constructor_move_Multiarray_i_dyn_extents
#define constructor_move_const_Multiarray_T_dyn_extents  constructor_move_const_Multiarray_i_dyn_extents
#define constructor_move_Multiarray_Vector_T_dyn_extents constructor_move_Multiarray_Vector_i_dyn_extents
#define constructor_move_Multiarray_Matrix_T_dyn_extents constructor_move_Multiarray_Matrix_i_dyn_extents
#define constructor_move_Multiarray_T_Matrix_T           constructor_move_Multiarray_i_Matrix_i
#define constructor_move_const_Multiarray_T_Matrix_T     constructor_move_const_Multiarray_i_Matrix_i
#define const_constructor_move_Multiarray_T              const_constructor_move_Multiarray_i
#define const_constructor_move_const_Multiarray_T        const_constructor_move_const_Multiarray_i
#define const_constructor_move_Multiarray_Vector_T       const_constructor_move_Multiarray_Vector_i
#define const_constructor_move_Multiarray_Matrix_T       const_constructor_move_Multiarray_Matrix_i

#define destructor_Multiarray_T               destructor_Multiarray_i
#define destructor_const_Multiarray_T         destructor_const_Multiarray_i
#define destructor_Multiarray_Vector_T        destructor_Multiarray_Vector_i
#define destructor_const_Multiarray_Vector_T  destructor_const_Multiarray_Vector_i
#define destructor_conditional_Multiarray_Vector_T       destructor_conditional_Multiarray_Vector_i
#define destructor_conditional_const_Multiarray_Vector_T destructor_conditional_const_Multiarray_Vector_i
#define destructor_Multiarray_Matrix_T        destructor_Multiarray_Matrix_i
#define destructor_const_Multiarray_Matrix_T  destructor_const_Multiarray_Matrix_i
#define destructor_const_Multiarray2_Matrix_T destructor_const_Multiarray2_Matrix_i
///\}

#endif

#endif
