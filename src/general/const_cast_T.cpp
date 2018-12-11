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

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void const_cast_T (const Type* dest, const Type src)
{
	*(Type*) dest = src;
}

void const_cast_T_n (const Type* dest, const Type* src, const int n)
{
	for (int i = 0; i < n; ++i)
		const_cast_T(&dest[i],src[i]);
}

void const_cast_T1 (const Type*const* dest, const Type*const src)
{
	*(Type**) dest = (Type*) src;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
