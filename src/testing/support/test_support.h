// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_support_h__INCLUDED
#define DPG__test_support_h__INCLUDED
/**	\file
 *	\brief Provides base support functions for testing.
 */

#include <stdbool.h>

/// \brief Update the value of pass based on the old/new values.
void update_pass
	(bool*const pass,     ///< The old value for `pass` (Modified in the function if necessary).
	 const bool pass_new  ///< The new value for `pass`.
	);

/**	\brief Compare the input `int`s and optionally print them if they differ.
 *	\return `true` if tests passed. */
bool compare_i
	(const int a,              ///< Input 0.
	 const int b,              ///< Input 1.
	 const bool print_enabled, ///< Flag for whether printing is enabled.
	 const char*const var_name ///< The variable name.
	);

/**	\brief Compare the input `bool`s and optionally print them if they differ.
 *	\return `true` if tests passed. */
bool compare_b
	(const bool a,             ///< Input 0.
	 const bool b,             ///< Input 1.
	 const bool print_enabled, ///< Flag for whether printing is enabled.
	 const char*const var_name ///< The variable name.
	);

/**	\brief Set the variable_name to be printed as part of one of the comparison functions based on a container member.
 *	\return See brief. */
char* set_print_name_container_member
	(const char*const name_container, ///< The name of the container.
	 int ind_container,               ///< The container index.
	 const char*const name_member     ///< The name of the member.
	);

#endif // DPG__test_support_h__INCLUDED
