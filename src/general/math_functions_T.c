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
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include "macros.h"
#include "definitions_tol.h"

#include "def_templates_math_functions.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

bool equal_T (const Type x0, const Type x1, const Real tol)
{
#if TYPE_RC == TYPE_REAL
	if ((fabs(x0) < tol && fabs(x0-x1) < tol) || (fabs((x0-x1)/x0) < tol))
#elif TYPE_RC == TYPE_COMPLEX
	if (equal_R(creal(x0),creal(x1),tol) && equal_R(cimag(x0),cimag(x1),tol))
#endif
		return true;
	return false;
}

Type norm_T (const ptrdiff_t n_entries, const Type*const data, const char*const norm_type)
{
	Type norm = 0.0;
	if (strcmp(norm_type,"L2") == 0) {
		for (ptrdiff_t i = 0; i < n_entries; ++i)
			norm += data[i]*data[i];
		return sqrt_T(norm);
	} else if (strcmp(norm_type,"Inf") == 0) {
		for (ptrdiff_t i = 0; i < n_entries; ++i) {
			if (abs_T(data[i]) > abs_T(norm))
				norm = abs_T(data[i]);
		}
		return norm;
	}
	EXIT_UNSUPPORTED;
}

Real norm_R_from_T (const ptrdiff_t n_entries, const Type*const data, const char*const norm_type)
{
	Real norm = 0.0;
	if (strcmp(norm_type,"L2") == 0) {
		for (ptrdiff_t i = 0; i < n_entries; ++i)
			norm += real_T(data[i])*real_T(data[i]);
		return sqrt_R(norm);
	} else if (strcmp(norm_type,"Inf") == 0) {
		for (ptrdiff_t i = 0; i < n_entries; ++i) {
			if (abs_R(real_T(data[i])) > abs_R(norm))
				norm = abs_R(real_T(data[i]));
		}
		return norm;
	}
	EXIT_UNSUPPORTED;
}

Type norm_diff_T
	(const ptrdiff_t n_entries, const Type*const data_0, const Type*const data_1, const char*const norm_type)
{
	Type norm_num = 0.0,
	     norm_den = 0.0;

	if (strstr(norm_type,"Inf")) {
		for (ptrdiff_t i = 0; i < n_entries; ++i) {
			const Type diff = abs_T(data_0[i]-data_1[i]);
			if (abs_T(diff) > abs_T(norm_num))
				norm_num = diff;

			const Type max = max_abs_T(data_0[i],data_1[i]);
			if (abs_T(max) > abs_T(norm_den))
				norm_den = max;
		}
	} else {
		EXIT_UNSUPPORTED;
	}

	return ( abs_T(norm_den) > 1e2*EPS ? norm_num/norm_den : norm_num );
}

Real norm_diff_RT
	(const ptrdiff_t n_entries, const Real*const data_0, const Type*const data_1, const char*const norm_type)
{
	Real norm_num = 0.0,
	     norm_den = 0.0;

	if (strstr(norm_type,"Inf")) {
		for (ptrdiff_t i = 0; i < n_entries; ++i) {
			const Real diff = abs_R(data_0[i]-real_T(data_1[i]));
			if (abs_R(diff) > abs_R(norm_num))
				norm_num = diff;

			const Real max = max_abs_R(data_0[i],real_T(data_1[i]));
			if (abs_R(max) > abs_R(norm_den))
				norm_den = max;
		}
	} else {
		EXIT_UNSUPPORTED;
	}

	return ( abs_R(norm_den) > 1e2*EPS ? norm_num/norm_den : norm_num );
}

Type max_abs_T (const Type a, const Type b)
{
	const Type a_abs = abs_T(a),
	           b_abs = abs_T(b);
	return ( abs_T(a_abs) > abs_T(b_abs) ? a_abs : b_abs );
}

void z_yxpz_T (const int n, const Type* x, const Type* y, Type* z)
{
	for (int i = 0; i < n; ++i)
		z[i] += y[i]*x[i];
}
#if TYPE_RC == TYPE_COMPLEX
void z_yxpz_RTT (const int n, const Real* x, const Type* y, Type* z)
{
	for (int i = 0; i < n; ++i)
		z[i] += y[i]*x[i];
}
#endif

Type average_T (const Type*const data, const ptrdiff_t n_entries)
{
	Type sum = 0.0;
	for (int i = 0; i < n_entries; ++i)
		sum += data[i];
	return sum/(double)n_entries;
}

Type minimum_T (const Type*const data, const ptrdiff_t n_entries)
{
	Type min = DBL_MAX;
	for (int i = 0; i < n_entries; ++i) {
		if (real_T(data[i]) < real_T(min))
			min = data[i];
	}
	return min;
}

Type maximum_abs_T (const Type*const data, const ptrdiff_t n_entries)
{
	Type max = DBL_MIN;
	for (int i = 0; i < n_entries; ++i) {
		if (abs_T(data[i]) > abs_T(max))
			max = abs_T(data[i]);
	}
	return max;
}

Real maximum_RT (const Type*const data, const ptrdiff_t n_entries)
{
	Real max = DBL_MIN;
	for (int i = 0; i < n_entries; ++i) {
		if (real_T(data[i]) > max)
			max = real_T(data[i]);
	}
	return max;
}

void add_to_T (Type*const data, const Type c_add, const ptrdiff_t n_entries)
{
	for (int i = 0; i < n_entries; ++i)
		data[i] += c_add;
}

Type dot_T (const ptrdiff_t n, const Type*const a, const Type*const b)
{
	Type res = 0.0;
	for (int i = 0; i < n; ++i)
		res += a[i]*b[i];
	return res;
}

Real dot_R_from_RT (const ptrdiff_t n, const Real*const a, const Type*const b)
{
	Real res = 0.0;
	for (int i = 0; i < n; ++i)
		res += a[i]*real_T(b[i]);
	return res;
}

Type max_abs_real_T (const Type a, const Type b)
{
	Type c = ( abs_T(a) > abs_T(b) ? a : b );
	if (real_T(c) < 0.0)
		c *= -1.0;
	return c;
}

Type min_abs_real_T (const Type a, const Type b)
{
	Type c = ( abs_T(a) < abs_T(b) ? a : b );
	if (real_T(c) < 0.0)
		c *= -1.0;
	return c;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_math_functions.h"
