// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_support_vector_h__INCLUDED
#define DPG__test_support_vector_h__INCLUDED
/**	\file
 *	Provides support functions for testing relating to the containers defined in \ref vector.h.
 */

#include <stdio.h>
#include <stdbool.h>

/** \brief Constructor for a \ref Vector_i\* from data in the input file of the given name.
 *	\return Standard. */
struct Vector_i* constructor_file_name_Vector_i
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief Constructor for a \ref Vector_i\* from the current line in the input file.
 *	\return Standard. */
struct Vector_i* constructor_file_Vector_i
	(FILE* data_file,           ///< The pointer to the file from which to read the data.
	 const bool check_container ///< Flag for whether the container type should be checked.
	);

/** \brief Check the difference between entries in the input \ref Vector_i\*s.
 *	\return The `true` if inputs differ; `false` otherwise. */
bool diff_Vector_i
	(const struct Vector_i*const a, ///< Input 0.
	 const struct Vector_i*const b  ///< Input 1.
	);

/// \brief Print the difference of the input \ref Vector_i\*s.
void print_diff_Vector_i
	(const struct Vector_i*const a, ///< Input 0.
	 const struct Vector_i*const b  ///< Input 1.
	);

#endif // DPG__test_support_vector_h__INCLUDED
