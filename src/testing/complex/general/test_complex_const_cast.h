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

#ifndef DPG__test_complex_const_cast_h__INCLUDED
#define DPG__test_complex_const_cast_h__INCLUDED
/** \file
 *  \brief Provides `complex` version of functions defined in \ref const_cast.h.
 */

#include <complex.h>

/// \brief `double complex` version of \ref const_cast_T.
void const_cast_dc
	(const double complex* dest, ///< See brief.
	 const double complex src    ///< See brief.
	);

/// \brief `double complex` version of \ref const_cast_T1.
void const_cast_dc1
	(const double complex*const* dest, ///< See brief.
	 const double complex*const  src   ///< See brief.
	);

#endif // DPG__test_complex_const_cast_h__INCLUDED
