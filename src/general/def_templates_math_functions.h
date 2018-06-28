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
 *  \brief Provides the macro definitions used for c-style templating related to the math functions.
 */

#if defined TYPE_RC
#if TYPE_RC == TYPE_REAL

///\{ \name Aliases for std library functions
#define abs_T  fabs
#define real_T
#define sqrt_T sqrt
#define pow_T  pow
#define log_T  log
///\}

///\{ \name Function names
#define equal_T     equal_d
#define norm_T      norm_d
#define norm_R_from_T norm_d_from_d
#define norm_diff_T norm_diff_d
#define norm_diff_RT norm_diff_dd
#define norm_diff_inf_no_rel_T norm_diff_inf_no_rel_d
#define max_abs_T   max_abs_d
#define z_yxpz_T    z_yxpz
#define z_yxpz_RTT  z_yxpz
#define average_T   average_d
#define minimum_T   minimum_d
#define maximum_abs_T maximum_abs_d
#define maximum_RT    maximum_dd
#define add_to_T    add_to_d
#define dot_T       dot_d
#define dot_R_from_RT dot_d_from_dd
#define min_abs_real_T min_abs_real_T
#define max_abs_real_T max_abs_real_T
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Aliases for std library functions
#define abs_T  cabs
#define real_T creal
#define sqrt_T csqrt
#define pow_T  cpow
#define log_T  clog
///\}

///\{ \name Function names
#define equal_T     equal_c
#define norm_T      norm_c
#define norm_R_from_T norm_d_from_c
#define norm_diff_T norm_diff_c
#define norm_diff_RT norm_diff_dc
#define norm_diff_inf_no_rel_T norm_diff_inf_no_rel_c
#define max_abs_T   max_abs_c
#define z_yxpz_T    z_yxpz_c
#define z_yxpz_RTT  z_yxpz_dcc
#define average_T   average_c
#define minimum_T   minimum_c
#define maximum_abs_T maximum_abs_c
#define maximum_RT    maximum_dc
#define add_to_T    add_to_c
#define dot_T       dot_c
#define dot_R_from_RT dot_d_from_dc
#define min_abs_real_T min_abs_real_T_c
#define max_abs_real_T max_abs_real_T_c
///\}

#endif

///\{ \name Real Data types/Function names
#define abs_R  fabs
#define sqrt_R sqrt
#define pow_R  pow

#define equal_R     equal_d
#define norm_R      norm_d
#define norm_diff_R norm_diff_d
#define dot_R       dot_d
#define max_abs_R   max_abs_d
///\}



#elif defined TYPE_I
#if TYPE_I == TYPE_II

///\{ \name Aliases for std library functions
#define abs_T  abs
///\}

#endif
#endif
