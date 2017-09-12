// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_support_cubature_h__INCLUDED
#define DPG__test_support_cubature_h__INCLUDED
/** \file
 *  \brief Provides support functions for testing relating to the \ref Cubature container.
 */

#include <stdbool.h>

// Constructor functions ******************************************************************************************** //

/** \brief Constructor for a \ref Cubature\* from data in the input file of the given name.
 *  \return Standard. */
struct Cubature* constructor_file_name_Cubature
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief `const` version of \ref constructor_file_name_Cubature.
 *  \return Standard. */
const struct const_Cubature* constructor_file_name_const_Cubature
	(const char*const var_name,      ///< Defined for \ref constructor_file_name_Cubature.
	 const char*const file_name_full ///< Defined for \ref constructor_file_name_Cubature.
	);

// Difference functions ********************************************************************************************* //

/** \brief Check the difference between members of the input \ref Cubature\*s.
 *  \return The `true` if inputs differ; `false` otherwise. */
bool diff_const_Cubature
	(const struct const_Cubature*const a, ///< Input 0.
	 const struct const_Cubature*const b, ///< Input 1.
	 const double tol                     ///< The tolerance.
	);

// Printing functions *********************************************************************************************** //

/// \brief Print the relative difference of the \ref Cubature members, outputting 0 if less than the tolerance.
void print_diff_const_Cubature
	(const struct const_Cubature*const a, ///< Input 0.
	 const struct const_Cubature*const b, ///< Input 1.
	 const double tol                     ///< The tolerance.
	);

#endif // DPG__test_support_cubature_h__INCLUDED
