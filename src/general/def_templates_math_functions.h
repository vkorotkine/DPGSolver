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

#if TYPE_RC == TYPE_REAL

///\{ \name Aliases for std library functions
#define abs_T  fabs
#define real_T
#define sqrt_T sqrt
#define pow_T  pow
///\}

///\{ \name Function names
#define equal_T     equal_d
#define norm_T      norm_d
#define norm_diff_T norm_diff_d
#define max_abs_T   max_abs_d
#define z_yxpz_T    z_yxpz
#define z_yxpz_RTT  z_yxpz
#define average_T   average_d
#define minimum_T   minimum_d
#define maximum_abs_T maximum_abs_d
#define add_to_T    add_to_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Aliases for std library functions
#define abs_T  cabs
#define real_T creal
#define sqrt_T csqrt
#define pow_T  cpow
///\}

///\{ \name Function names
#define equal_T     equal_c
#define norm_T      norm_c
#define norm_diff_T norm_diff_c
#define max_abs_T   max_abs_c
#define z_yxpz_T    z_yxpz_c
#define z_yxpz_RTT  z_yxpz_dcc
#define average_T   average_c
#define minimum_T   minimum_c
#define maximum_abs_T maximum_abs_c
#define add_to_T    add_to_c
///\}

#endif
