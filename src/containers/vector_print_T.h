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
 *  \brief Provides Vector_\* printing functions.
 */

#include <stdio.h>

struct Vector_T;
struct const_Vector_T;

/// \brief Print a \ref Vector_T\* to the terminal displaying entries below the tolerance as 0.0.
void print_Vector_T_tol
	(const struct Vector_T*const a, ///< Standard.
	 const Real tol                 ///< The tolerance.
	);

/// \brief `const` version of \ref print_Vector_T_tol.
void print_const_Vector_T_tol
	(const struct const_Vector_T*const a, ///< Defined for \ref print_Vector_T_tol.
	 const Real tol                       ///< Defined for \ref print_Vector_T_tol.
	);

/// \brief Print a \ref Vector_T\* calling \ref print_Vector_T_tol with a default tolerance.
void print_Vector_T
	(const struct Vector_T*const a ///< Defined for \ref print_Vector_T_tol.
	);

/// \brief `const` version of \ref print_Vector_T.
void print_const_Vector_T
	(const struct const_Vector_T*const a ///< Defined for \ref print_Vector_T.
	);

#ifndef TYPE_RC
/// \brief Print a \ref const_Vector_T to a file with the input number of tabs before each row.
void fprint_const_Vector_T
	(FILE* file,                    ///< The file.
	 const int n_tab,               ///< The number of tabs.
	 const struct const_Vector_T* a ///< Standard.
	);

/// \brief `mutable` version of \ref fprint_const_Vector_T.
void fprint_Vector_T
	(FILE* file,        ///< Defined for \ref fprint_const_Vector_T.
	 const int n_tab,   ///< Defined for \ref fprint_const_Vector_T.
	 struct Vector_T* a ///< Defined for \ref fprint_const_Vector_T.
	);
#endif
