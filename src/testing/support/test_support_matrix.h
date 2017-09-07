// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_support_matrix_h__INCLUDED
#define DPG__test_support_matrix_h__INCLUDED
/**	\file
 *	\brief Provides support functions for testing relating to the containers defined in \ref matrix.h.
 */

#include <stdio.h>
#include <stdbool.h>

struct const_Matrix_d;

/** \brief Constructor for a \ref Matrix_i\* from data in the input file of the given name.
 *	\return Standard. */
struct Matrix_i* constructor_file_name_Matrix_i
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief Constructor for a \ref Matrix_d\* from data in the input file of the given name.
 *	\return Standard. */
struct Matrix_d* constructor_file_name_Matrix_d
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief Constructor for a \ref Matrix_i\* from the current line in the input file.
 *	\return Standard. */
struct Matrix_i* constructor_file_Matrix_i
	(FILE* data_file,           ///< The pointer to the file from which to read the data.
	 const bool check_container ///< Flag for whether the container type should be checked.
	);

/** \brief Constructor for a \ref Matrix_d\* from the current line in the input file.
 *	\return Standard. */
struct Matrix_d* constructor_file_Matrix_d
	(FILE* data_file,           ///< The pointer to the file from which to read the data.
	 const bool check_container ///< Flag for whether the container type should be checked.
	);

/** \brief Check the difference between entries in the input \ref Matrix_i\*s.
 *	\return The `true` if inputs differ; `false` otherwise. */
bool diff_Matrix_i
	(const struct Matrix_i*const a, ///< Input 0.
	 const struct Matrix_i*const b  ///< Input 1.
	);

/** \brief Check the relative difference between entries in the input \ref Matrix_d\*s up to the input tolerance.
 *	\return The `true` if inputs differ; `false` otherwise. */
bool diff_Matrix_d
	(const struct Matrix_d*const a, ///< Input 0.
	 const struct Matrix_d*const b, ///< Input 1.
	 const double tol               ///< The tolerance.
	);

/** \brief `const` version of \ref diff_Matrix_d.
 *	\return See brief. */
bool diff_const_Matrix_d
	(const struct const_Matrix_d*const a, ///< Input 0.
	 const struct const_Matrix_d*const b, ///< Input 1.
	 const double tol                     ///< The tolerance.
	);

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
