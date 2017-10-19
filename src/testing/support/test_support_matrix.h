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

#ifndef DPG__test_support_matrix_h__INCLUDED
#define DPG__test_support_matrix_h__INCLUDED
/** \file
 *  \brief Provides support functions for testing relating to the containers defined in \ref matrix.h.
 */

#include <stdio.h>
#include <stdbool.h>

struct const_Vector_d;
struct const_Matrix_d;

// Constructor functions ******************************************************************************************** //

/** \brief Constructor for a \ref Matrix_i\* from data in the input file of the given name.
 *  \return Standard. */
struct Matrix_i* constructor_file_name_Matrix_i
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief Constructor for a \ref Matrix_d\* from data in the input file of the given name.
 *  \return Standard. */
struct Matrix_d* constructor_file_name_Matrix_d
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief Constructor for a \ref const_Matrix_d\* from data in the input file of the given name.
 *  \return Standard. */
const struct const_Matrix_d* constructor_file_name_const_Matrix_d
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief Constructor for a \ref Matrix_i\* from the current line in the input file.
 *  \return Standard. */
struct Matrix_i* constructor_file_Matrix_i
	(FILE* data_file,           ///< The pointer to the file from which to read the data.
	 const bool check_container ///< Flag for whether the container type should be checked.
	);

/** \brief Constructor for a \ref Matrix_d\* from the current line in the input file.
 *  \return Standard. */
struct Matrix_d* constructor_file_Matrix_d
	(FILE* data_file,           ///< The pointer to the file from which to read the data.
	 const bool check_container ///< Flag for whether the container type should be checked.
	);

/** \brief `const` version of \ref constructor_copy_transpose_Matrix_d.
 *  \return Standard. */
const struct const_Matrix_d* constructor_copy_transpose_const_Matrix_d
	(const struct const_Matrix_d* a, ///< The input matrix.
	 const bool mem_only             ///< Defined for \ref transpose_Matrix_d.
	);

// Math functions *************************************************************************************************** //

/// \brief `const` version of \ref transpose_Matrix_d.
void transpose_const_Matrix_d
	(const struct const_Matrix_d* a, ///< Defined for \ref transpose_Matrix_d.
	 const bool mem_only             ///< Defined for \ref transpose_Matrix_d.
	);

/// \brief `const` version of \ref scale_Matrix_by_Vector_d.
void scale_const_Matrix_by_Vector_d
	(const char side,                     ///< Defined for \ref scale_Matrix_by_Vector_d.
	 const double alpha,                  ///< Defined for \ref scale_Matrix_by_Vector_d.
	 const struct const_Matrix_d*const a, ///< Defined for \ref scale_Matrix_by_Vector_d.
	 const struct const_Vector_d*const b, ///< Defined for \ref scale_Matrix_by_Vector_d.
	 const bool invert_diag               ///< Defined for \ref scale_Matrix_by_Vector_d.
	);

// Math functions *************************************************************************************************** //

// Difference functions ********************************************************************************************* //

/** \brief Check the difference between entries in the input \ref Matrix_i\*s.
 *  \return The `true` if inputs differ; `false` otherwise. */
bool diff_Matrix_i
	(const struct Matrix_i*const a, ///< Input 0.
	 const struct Matrix_i*const b  ///< Input 1.
	);

/** \brief Check the relative difference between entries in the input \ref Matrix_d\*s up to the input tolerance.
 *  \return The `true` if inputs differ; `false` otherwise. */
bool diff_Matrix_d
	(const struct Matrix_d*const a, ///< Input 0.
	 const struct Matrix_d*const b, ///< Input 1.
	 const double tol               ///< The tolerance.
	);

/** \brief `const` version of \ref diff_Matrix_d.
 *  \return See brief. */
bool diff_const_Matrix_d
	(const struct const_Matrix_d*const a, ///< Input 0.
	 const struct const_Matrix_d*const b, ///< Input 1.
	 const double tol                     ///< The tolerance.
	);

// Printing functions *********************************************************************************************** //

/// \brief Print the difference of the input \ref Matrix_i\*s.
void print_diff_Matrix_i
	(const struct Matrix_i*const a, ///< Input 0.
	 const struct Matrix_i*const b  ///< Input 1.
	);

/// \brief Print the relative difference of the input \ref Matrix_d\*s, outputting 0 if less than the tolerance.
void print_diff_Matrix_d
	(const struct Matrix_d*const a, ///< Input 0.
	 const struct Matrix_d*const b, ///< Input 1.
	 const double tol               ///< The tolerance.
	);

/// \brief `const` version of \ref print_diff_Matrix_d.
void print_diff_const_Matrix_d
	(const struct const_Matrix_d*const a, ///< Input 0.
	 const struct const_Matrix_d*const b, ///< Input 1.
	 const double tol                     ///< The tolerance.
	);

#endif // DPG__test_support_matrix_h__INCLUDED
