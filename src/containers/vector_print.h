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

#ifndef DPG__vector_print_h__INCLUDED
#define DPG__vector_print_h__INCLUDED
/** \file
 *  \brief Provides Vector_\* printing functions.
 */

#include <stdio.h>

struct Vector_i;
struct Vector_d;
struct const_Vector_i;
struct const_Vector_d;

/// \brief Print a \ref Vector_i\* to the terminal.
void print_Vector_i
	(const struct Vector_i*const a ///< Standard.
	);

/// \brief Print a \ref const_Vector_i\* to the terminal.
void print_const_Vector_i
	(const struct const_Vector_i*const a ///< Standard.
	);

/// \brief Print a \ref Vector_d\* to the terminal displaying entries below the tolerance as 0.0.
void print_Vector_d_tol
	(const struct Vector_d*const a, ///< Standard.
	 const double tol               ///< The tolerance.
	);

/// \brief `const` version of \ref print_Vector_d_tol.
void print_const_Vector_d_tol
	(const struct const_Vector_d*const a, ///< Defined for \ref print_Vector_d_tol.
	 const double tol                     ///< Defined for \ref print_Vector_d_tol.
	);

/// \brief Print a \ref Vector_d\* calling \ref print_Vector_d_tol with a default tolerance.
void print_Vector_d
	(const struct Vector_d*const a ///< Defined for \ref print_Vector_d_tol.
	);

/// \brief `const` version of \ref print_Vector_d.
void print_const_Vector_d
	(const struct const_Vector_d*const a ///< Defined for \ref print_Vector_d.
	);

/// \brief Print a \ref const_Vector_i to a file with the input number of tabs before each row.
void fprint_const_Vector_i
	(FILE* file,                    ///< The file.
	 const int n_tab,               ///< The number of tabs.
	 const struct const_Vector_i* a ///< Standard.
	);

/// \brief `mutable` version of \ref fprint_const_Vector_i.
void fprint_Vector_i
	(FILE* file,        ///< Defined for \ref fprint_const_Vector_i.
	 const int n_tab,   ///< Defined for \ref fprint_const_Vector_i.
	 struct Vector_i* a ///< Defined for \ref fprint_const_Vector_i.
	);

#endif // DPG__vector_print_h__INCLUDED
