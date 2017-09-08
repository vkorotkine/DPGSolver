// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__vector_print_h__INCLUDED
#define DPG__vector_print_h__INCLUDED
/**	\file
 *	\brief Provides Vector_\* printing functions.
 */

struct Vector_i;
struct Vector_d;
struct const_Vector_i;

/// \brief Print a \ref Vector_i\* to the terminal.
void print_Vector_i
	(const struct Vector_i*const a ///< Standard.
	);

/// \brief Print a \ref const_Vector_i\* to the terminal.
void print_const_Vector_i
	(const struct const_Vector_i*const a ///< Standard.
	);

/// \brief Print a \ref Vector_d\* to the terminal displaying entries below the tolerance as 0.0.
void print_Vector_d
	(const struct Vector_d*const a, ///< Standard.
	 const double tol               ///< The tolerance.
	);

#endif // DPG__vector_print_h__INCLUDED
