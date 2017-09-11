// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_support_multiarray_h__INCLUDED
#define DPG__test_support_multiarray_h__INCLUDED
/** \file
 *  \brief Provides support functions for testing relating to the containers defined in \ref multiarray.h.
 */

#include <stdio.h>
#include <stdbool.h>

struct const_Multiarray_d;

// Constructor functions ******************************************************************************************** //

/** \brief Constructor for a \ref Multiarray_d\* from data in the input file of the given name.
 *  \return Standard. */
struct Multiarray_d* constructor_file_name_Multiarray_d
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief Constructor for a \ref Multiarray_Vector_i\* from data in the input file of the given name.
 *  \return Standard. */
struct Multiarray_Vector_i* constructor_file_name_Multiarray_Vector_i
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief Constructor for a \ref Multiarray_d\* from the current line in the input file.
 *  \return Standard. */
struct Multiarray_d* constructor_file_Multiarray_d
	(FILE* data_file,           ///< The pointer to the file from which to read the data.
	 const bool check_container ///< Flag for whether the container type should be checked.
	);

// Difference functions ********************************************************************************************* //

/** \brief Check the difference between entries in the input \ref Multiarray_Vector_i\*s.
 *  \return The `true` if inputs differ; `false` otherwise. */
bool diff_Multiarray_Vector_i
	(const struct Multiarray_Vector_i*const a, ///< Input 0.
	 const struct Multiarray_Vector_i*const b  ///< Input 1.
	);

/** \brief Check the relative difference between entries in the input \ref Multiarray_d\*s up to the input tolerance.
 *  \return The `true` if inputs differ; `false` otherwise. */
bool diff_Multiarray_d
	(const struct Multiarray_d*const a, ///< Input 0.
	 const struct Multiarray_d*const b, ///< Input 1.
	 const double tol                   ///< The tolerance.
	);

/** \brief `const` version of \ref diff_Multiarray_d.
 *  \return See brief. */
bool diff_const_Multiarray_d
	(const struct const_Multiarray_d*const a, ///< Defined for \ref diff_Multiarray_d.
	 const struct const_Multiarray_d*const b, ///< Defined for \ref diff_Multiarray_d.
	 const double tol                         ///< Defined for \ref diff_Multiarray_d.
	);

// Printing functions *********************************************************************************************** //

/// \brief Print the difference of the input \ref Multiarray_Vector_i\*s.
void print_diff_Multiarray_Vector_i
	(const struct Multiarray_Vector_i*const a, ///< Input 0.
	 const struct Multiarray_Vector_i*const b  ///< Input 1.
	);

/// \brief Print the relative difference of the input \ref Multiarray_d\*s, outputting 0 if less than the tolerance.
void print_diff_Multiarray_d
	(const struct Multiarray_d*const a, ///< Input 0.
	 const struct Multiarray_d*const b, ///< Input 1.
	 const double tol                   ///< The tolerance.
	);

/// \brief `const` version of \ref print_diff_Multiarray_d.
void print_diff_const_Multiarray_d
	(const struct const_Multiarray_d*const a, ///< Input 0.
	 const struct const_Multiarray_d*const b, ///< Input 1.
	 const double tol                         ///< The tolerance.
	);

#endif // DPG__test_support_multiarray_h__INCLUDED
