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

#ifndef DPG__const_cast_h__INCLUDED
#define DPG__const_cast_h__INCLUDED
/** \file
 *  \brief Provides functions for casting to `const` lvalues for standard datatypes.
 *
 *  \note These functions should only be used for dynamically allocated `const` variables.
 */

#include <stddef.h>
#include <stdbool.h>

/// \brief `char` version of \ref const_cast_T.
void const_cast_c
	(const char* dest, ///< See brief.
	 const char src    ///< See brief.
	);

/// \brief `char` version of \ref const_cast_T1.
void const_cast_c1
	(const char*const* dest, ///< See brief.
	 const char*const  src   ///< See brief.
	);

/// \brief `int` version of \ref const_cast_T.
void const_cast_i
	(const int* dest, ///< See brief.
	 const int src    ///< See brief.
	);

/// \brief `int` version of \ref const_cast_T_n.
void const_cast_i_n
	(const int* dest, ///< See brief.
	 const int* src,  ///< See brief.
	 const int n_src  ///< See brief.
	);

/// \brief `double` version of \ref const_cast_T.
void const_cast_d
	(const double* dest, ///< See brief.
	 const double src    ///< See brief.
	);

/// \brief `double` version of \ref const_cast_T1.
void const_cast_d1
	(const double*const* dest, ///< See brief.
	 const double*const  src   ///< See brief.
	);

/// \brief `ptrdiff_t` version of \ref const_cast_T.
void const_cast_ptrdiff
	(const ptrdiff_t* dest, ///< See brief.
	 const ptrdiff_t src    ///< See brief.
	);

/// \brief `bool` version of \ref const_cast_T.
void const_cast_b
	(const bool* dest, ///< See brief.
	 const bool src    ///< See brief.
	);

#endif // DPG__const_cast_h__INCLUDED
