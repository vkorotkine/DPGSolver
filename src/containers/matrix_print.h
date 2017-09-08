// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__matrix_print_h__INCLUDED
#define DPG__matrix_print_h__INCLUDED
/**	\file
 *	\brief Provides Matrix_\* printing functions.
 */

struct Matrix_d;
struct Matrix_i;
struct const_Matrix_d;
struct const_Matrix_i;

/// \brief Print a \ref Matrix_d\* to the terminal displaying entries below the tolerance as 0.0.
void print_Matrix_d
	(const struct Matrix_d*const a, ///< Standard.
	 const double tol               ///< The tolerance.
	);

/// \brief Print a \ref const_Matrix_d\* as in \ref print_Matrix_d.
void print_const_Matrix_d
	(const struct const_Matrix_d*const a, ///< Standard.
	 const double tol                     ///< The tolerance.
	);

/// \brief Print a \ref Matrix_i\* to the terminal.
void print_Matrix_i
	(const struct Matrix_i*const a ///< Standard.
	);

/// \brief Print a \ref const_Matrix_i\* to the terminal.
void print_const_Matrix_i
	(const struct const_Matrix_i*const a ///< Standard.
	);

#endif // DPG__matrix_print_h__INCLUDED
