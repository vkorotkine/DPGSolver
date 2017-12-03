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

#ifndef DPG__matrix_print_h__INCLUDED
#define DPG__matrix_print_h__INCLUDED
/** \file
 *  \brief Provides Matrix_\* printing functions.
 */

struct Matrix_d;
struct Matrix_i;
struct const_Matrix_d;
struct const_Matrix_i;

/// \brief Print a \ref Matrix_T\* to the terminal displaying entries below the tolerance as 0.0.
void print_Matrix_d_tol
	(const struct Matrix_d*const a, ///< Standard.
	 const double tol               ///< The tolerance.
	);

/// \brief `const` version of \ref print_Matrix_d_tol.
void print_const_Matrix_d_tol
	(const struct const_Matrix_d*const a, ///< See brief.
	 const double tol                     ///< See brief.
	);

/// \brief Print a \ref Matrix_T\* calling \ref print_Matrix_d_tol with a default tolerance.
void print_Matrix_d
	(const struct Matrix_d*const a ///< Defined for \ref print_Matrix_d_tol.
	);

/// \brief `const` version of \ref print_Matrix_d.
void print_const_Matrix_d
	(const struct const_Matrix_d*const a ///< Defined for \ref print_Matrix_d.
	);

/// \brief Print a \ref Matrix_T\* to the terminal.
void print_Matrix_i
	(const struct Matrix_i*const a ///< Standard.
	);

/// \brief Print a \ref const_Matrix_T\* to the terminal.
void print_const_Matrix_i
	(const struct const_Matrix_i*const a ///< Standard.
	);

#endif // DPG__matrix_print_h__INCLUDED
