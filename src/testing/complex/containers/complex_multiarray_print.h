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

#ifndef DPG__complex_multiarray_print_h__INCLUDED
#define DPG__complex_multiarray_print_h__INCLUDED
/** \file
 *  \brief Provides \ref Multiarray_c printing functions.
 */

struct Multiarray_c;
struct const_Multiarray_c;

/// \brief `complex` version of \ref print_Multiarray_c.
void print_Multiarray_c
	(const struct Multiarray_c*const a ///< See brief.
	);

/// \brief `const` version of \ref print_Multiarray_c.
void print_const_Multiarray_c
	(const struct const_Multiarray_c*const a ///< See brief.
	);

#endif // DPG__complex_multiarray_print_h__INCLUDED
