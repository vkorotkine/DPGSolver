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
 *  \brief Provides Matrix_\* printing functions.
 */

#include <stdio.h>
#include <complex.h>

#include "def_templates_matrix.h"

struct Matrix_T;
struct const_Matrix_T;

/// \brief Print a \ref Matrix_T\* to the terminal displaying entries below the tolerance as 0.0.
void print_Matrix_T_tol
	(const struct Matrix_T*const a, ///< Standard.
	 const Real tol                 ///< The tolerance.
	);

/// \brief `const` version of \ref print_Matrix_T_tol.
void print_const_Matrix_T_tol
	(const struct const_Matrix_T*const a, ///< See brief.
	 const Real tol                       ///< See brief.
	);

/// \brief Print a \ref Matrix_T\* calling \ref print_Matrix_T_tol with a default tolerance.
void print_Matrix_T
	(const struct Matrix_T*const a ///< Defined for \ref print_Matrix_T_tol.
	);

/// \brief `const` version of \ref print_Matrix_T.
void print_const_Matrix_T
	(const struct const_Matrix_T*const a ///< Defined for \ref print_Matrix_T.
	);

/// \brief Print a \ref const_Matrix_T\* to a file.
void print_to_file_const_Matrix_T
	(FILE*const file,                    ///< The file.
	 const struct const_Matrix_T*const a ///< The matrix to print.
	);

/// \brief `mutable` version of \ref print_to_file_const_Matrix_T.
void print_to_file_Matrix_T
	(FILE*const file,              ///< The file.
	 const struct Matrix_T*const a ///< The matrix to print.
	);

#if TYPE_RC == TYPE_COMPLEX
/// \brief Print a real value to the terminal with default format.
void print_real
	(const double complex val ///< The complex value.
	);

/// \brief Print an imaginary value to the terminal with default format.
void print_imag
	(const double complex val ///< The complex value.
	);
#endif

#include "undef_templates_matrix.h"
