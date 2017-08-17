// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Vector_h__INCLUDED
#define DPG__Vector_h__INCLUDED
/**	\file
 *	\brief Provides Vector_\* containers and related functions.
 *
 *	Potentially relevant comments may be found in containers.h.
 *
 *	Vectors are 1D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>

/// \brief Vector (`unsigned int`).
struct Vector_ui {
	size_t extents[1];

	bool    owns_data;
	unsigned int* data;
};

/// \brief Vector (`const unsigned int`).
struct const_Vector_ui {
	const size_t extents[1];

	const bool    owns_data;
	const unsigned int* data;
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructs a default \ref Vector_ui\*\*.
struct Vector_ui** constructor_default_Vector_ui_2
	(const size_t n_dest ///< The number of \ref Vector_ui\* components.
	);

/// \brief Constructs an empty \ref Vector_ui\*.
struct Vector_ui* constructor_empty_Vector_ui
	(const size_t ext_0 ///< The value of extents[0].
	);

/// \brief Copy constructor for a \ref Vector_ui\* from a `const unsigned int*`.
struct Vector_ui* constructor_copy_Vector_ui_ui
	(const size_t ext_0,               ///< The value of extents[0].
	 const unsigned int*const data_src ///< The source data.
	);

/// \brief Move constructor for a \ref Vector_ui\* from a `unsigned int*`.
struct Vector_ui* constructor_move_Vector_ui_ui
	(const size_t ext_0,     ///< The value of extents[0].
	 const bool owns_data,   ///< Standard.
	 unsigned int*const data ///< Standard.
	);

/// \brief Move constructor for a \ref const_Vector_ui\* from a `const unsigned int*`.
struct const_Vector_ui* constructor_move_const_Vector_ui_ui
	(const size_t ext_0,           ///< The value of extents[0].
	 const bool owns_data,         ///< Standard.
	 const unsigned int*const data ///< Standard.
	);

/// \brief Move constructor for a `const` \ref const_Vector_ui `*const`.
void const_constructor_move_Vector_ui
	(const struct const_Vector_ui*const* dest, ///< Destination.
	 struct Vector_ui* src                     ///< Source.
	);

/// \brief Destructs a \ref Vector_ui\*.
void destructor_Vector_ui
	(struct Vector_ui* a ///< Standard.
	);

/// \brief Destructs a \ref Vector_ui\*\*.
void destructor_Vector_ui_2
	(struct Vector_ui** a, ///< Standard.
	 const size_t n_src    ///< The number of \ref Vector_ui\* components.
	);

// Helper functions ************************************************************************************************* //

/** \brief Reorder a \ref Vector_ui based on the provided ordering.
 *	\warning This is not currently done in place.
 */
void reorder_Vector_ui
	(struct Vector_ui*const a,         ///< Standard.
	 const unsigned int*const ordering ///< The ordering.
	);

/// \brief Reserve space for a \ref Vector_ui\*.
void reserve_Vector_ui
	(struct Vector_ui*const a, ///< Standard.
	 const size_t ext_0        ///< New value for extents[0].
	);

/// \brief Set all data entries to zero.
void set_to_zero_Vector_ui
	(struct Vector_ui*const a ///< Standard.
	);

/// \brief Sort the data of the \ref Vector_ui\*.
void sort_Vector_ui
	(struct Vector_ui* a ///< Standard.
	);

/** \brief See return.
 *	\return The sum of the components of the \ref Vector_ui\*.
 */
unsigned int sum_Vector_ui
	(struct Vector_ui* a ///< Standard.
	);

/** \brief See return.
 *	\return `bool` indicating whether the \ref Vector_ui\* inputs are equal.
 */
bool check_equal_Vector_ui
	(const struct Vector_ui*const a, ///< Input 1.
	 const struct Vector_ui*const b  ///< Input 2.
	);

/** \brief Comparison function for std::qsort between \ref Vector_ui\*\* `a` and `b`.
 *	\return The lexicographical comparison of `a` and `b`.
 *
 *	\note Input Vectors must be have sorted data.
 */
int cmp_Vector_ui
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

/// \brief Copy the data of the `src` to that of the `dest` \ref Vector_ui\*.
void copy_data_Vector_ui_Vector_ui
	(const struct Vector_ui*const src, ///< Standard.
	 struct Vector_ui*const dest       ///< Standard.
	);

// Printing functions *********************************************************************************************** //

/// \brief Print a \ref Vector_ui\* to the terminal.
void print_Vector_ui
	(const struct Vector_ui*const a ///< Standard.
	);

/// \brief Print a \ref const_Vector_ui\* to the terminal.
void print_const_Vector_ui
	(const struct const_Vector_ui*const a ///< Standard.
	);

#endif // DPG__Vector_h__INCLUDED
