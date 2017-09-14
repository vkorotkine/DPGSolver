// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__multiarray_print_h__INCLUDED
#define DPG__multiarray_print_h__INCLUDED
/** \file
 *  \brief Provides Multiarray_\* printing functions.
 */

struct Multiarray_d;
struct Multiarray_Vector_i;
struct Multiarray_Matrix_d;
struct const_Multiarray_d;
struct const_Multiarray_Vector_i;
struct const_Multiarray_Matrix_d;

/// \brief Print a \ref Multiarray_Vector_i\* to the terminal.
void print_Multiarray_Vector_i
	(const struct Multiarray_Vector_i*const a ///< Standard.
	);

/// \brief `const` version of \ref print_Multiarray_Vector_i.
void print_const_Multiarray_Vector_i
	(const struct const_Multiarray_Vector_i*const a ///< Defined for \ref print_Multiarray_Vector_i.
	);

/// \brief Print a \ref Multiarray_d\* to the terminal displaying entries below the tolerance as 0.0.
void print_Multiarray_d
	(const struct Multiarray_d*const a, ///< Standard.
	 const double tol                   ///< The tolerance.
	);

/// \brief `const` version of \ref print_Multiarray_d.
void print_const_Multiarray_d
	(const struct const_Multiarray_d*const a, ///< Defined for \ref print_Multiarray_d.
	 const double tol                         ///< Defined for \ref print_Multiarray_d.
	);

/// \brief Print a \ref Multiarray_Matrix_d\* to the terminal.
void print_Multiarray_Matrix_d
	(const struct Multiarray_Matrix_d*const a, ///< Standard.
	 const double tol                          ///< The tolerance.
	);

/// \brief `const` version of \ref print_Multiarray_Matrix_d.
void print_const_Multiarray_Matrix_d
	(const struct const_Multiarray_Matrix_d*const a, ///< Defined for \ref print_Multiarray_Matrix_d.
	 const double tol                                ///< Defined for \ref print_Multiarray_Matrix_d.
	);

#endif // DPG__multiarray_print_h__INCLUDED
