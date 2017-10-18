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
/// \file

#include "const_cast.h"

// Standard data types ********************************************************************************************** //

void const_cast_i (const int* dest, const int src)
{
	*(int*) dest = src;
}

void const_cast_i1 (const int* dest, const int* src, const int n_src)
{
	for (int i = 0; i < n_src; ++i)
		const_cast_i(&dest[i],src[i]);
}

void const_cast_ptrdiff (const ptrdiff_t* dest, const ptrdiff_t src)
{
	*(ptrdiff_t*) dest = src;
}

void const_cast_bool (const bool* dest, const bool src)
{
	const_cast_b(dest,src);
}

void const_cast_b (const bool* dest, const bool src)
{
	*(bool*) dest = src;
}

void const_cast_b1 (const bool* dest, const bool* src, const int n_src)
{
	for (int i = 0; i < n_src; ++i)
		const_cast_bool(&dest[i],src[i]);
}

void const_cast_c (const char* dest, const char src)
{
	*(char*) dest = src;
}

void const_cast_c1 (const char*const* dest, const char*const src)
{
	*(char**) dest = (char*) src;
}

void const_cast_d (const double* dest, const double src)
{
	*(double*) dest = src;
}

void const_cast_d1 (const double*const* dest, const double*const src)
{
	*(double**) dest = (double*) src;
}

void const_cast_void1 (const void*const* dest, const void*const src)
{
	*(void**) dest = (void*) src;
}

// Custom data types ************************************************************************************************ //
