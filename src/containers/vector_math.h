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

#ifndef DPG__vector_math_h__INCLUDED
#define DPG__vector_math_h__INCLUDED
/** \file
 *  \brief Provides Vector_\* math functions.
 */

struct Vector_d;

/// \brief Invert each of the entries of the input \ref Vector_d\*.
void invert_Vector_d
	(struct Vector_d* a ///< Input vector.
	);

#endif // DPG__vector_math_h__INCLUDED
