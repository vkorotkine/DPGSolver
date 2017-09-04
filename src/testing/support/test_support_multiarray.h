// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_support_multiarray_h__INCLUDED
#define DPG__test_support_multiarray_h__INCLUDED
/**	\file
 *	Provides support functions for testing relating to the containers defined in \ref multiarray.h.
 */

#include <stdio.h>
#include <stdbool.h>

/** \brief Constructor for a \ref Multiarray_Vector_i\* from data in the input file of the given name.
 *	\return Standard. */
struct Multiarray_Vector_i* constructor_file_name_Multiarray_Vector_i
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/// \brief Check that the container type is that which is expected.
void check_container_type
	(FILE* data_file,                ///< The file containing the data.
	 const char*const container_type ///< The container type.
	);

/** \brief Check the difference between entries in the input \ref Multiarray_Vector_i\*s.
 *	\return The `true` if inputs differ; `false` otherwise. */
bool diff_Multiarray_Vector_i
	(const struct Multiarray_Vector_i*const a, ///< Input 0.
	 const struct Multiarray_Vector_i*const b  ///< Input 1.
	);

/// \brief Print the difference of the input \ref Multiarray_Vector_i\*s.
void print_diff_Multiarray_Vector_i
	(const struct Multiarray_Vector_i*const a, ///< Input 0.
	 const struct Multiarray_Vector_i*const b  ///< Input 1.
	);

#endif // DPG__test_support_multiarray_h__INCLUDED
