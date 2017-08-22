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

/// \brief Vector (`int`).
struct Vector_i {
	ptrdiff_t extents[1];

	bool owns_data;
	int* data;
};

/// \brief Vector (`const int`).
struct const_Vector_i {
	const ptrdiff_t extents[1];

	const bool owns_data;
	const int* data;
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructs a default \ref Vector_i\*\*.
struct Vector_i** constructor_default_Vector_i_2
	(const ptrdiff_t n_dest ///< The number of \ref Vector_i\* components.
	);

/// \brief Constructs an empty \ref Vector_i\*.
struct Vector_i* constructor_empty_Vector_i
	(const ptrdiff_t ext_0 ///< The value of extents[0].
	);

/// \brief Copy constructor for a \ref Vector_i\* from a `Vector_i*`.
struct Vector_i* constructor_copy_Vector_i
	(const struct Vector_i*const src ///< The source data.
	);

/// \brief Copy constructor for a \ref Vector_i\* from a `const int*`.
struct Vector_i* constructor_copy_Vector_i_i
	(const ptrdiff_t ext_0,   ///< The value of extents[0].
	 const int*const data_src ///< The source data.
	);

/// \brief Move constructor for a \ref Vector_i\* from a `int*`.
struct Vector_i* constructor_move_Vector_i_i
	(const ptrdiff_t ext_0, ///< The value of extents[0].
	 const bool owns_data,  ///< Standard.
	 int*const data         ///< Standard.
	);

/// \brief Move constructor for a \ref const_Vector_i\* from a `const int*`.
struct const_Vector_i* constructor_move_const_Vector_i_i
	(const ptrdiff_t ext_0, ///< The value of extents[0].
	 const bool owns_data,  ///< Standard.
	 const int*const data   ///< Standard.
	);

/// \brief Move constructor for a `const` \ref const_Vector_i `*const`.
void const_constructor_move_Vector_i
	(const struct const_Vector_i*const* dest, ///< Destination.
	 struct Vector_i* src                     ///< Source.
	);

/// \brief Destructs a \ref Vector_i\*.
void destructor_Vector_i
	(struct Vector_i* a ///< Standard.
	);

/// \brief Destructs a \ref Vector_i\*\*.
void destructor_Vector_i_2
	(struct Vector_i** a,   ///< Standard.
	 const ptrdiff_t n_src, ///< The number of \ref Vector_i\* components.
	 const bool owns_data   ///< Standard.
	);

// Helper functions ************************************************************************************************* //

/** \brief Reorder a \ref Vector_i based on the provided ordering.
 *	\warning This is not currently done in place.
 */
void reorder_Vector_i
	(struct Vector_i*const a, ///< Standard.
	 const int*const ordering ///< The ordering.
	);

/// \brief Reserve space for a \ref Vector_i\*.
void reserve_Vector_i
	(struct Vector_i*const a, ///< Standard.
	 const ptrdiff_t ext_0    ///< New value for extents[0].
	);

/// \brief Set all data entries to zero.
void set_to_zero_Vector_i
	(struct Vector_i*const a ///< Standard.
	);

/// \brief Set data entries to those of the src data.
void set_to_data_Vector_i
	(struct Vector_i*const a, ///< Standard.
	 const int*const data_src ///< The source data.
	);

/// \brief Sort the data of the \ref Vector_i\*.
void sort_Vector_i
	(struct Vector_i* a ///< Standard.
	);

/** \brief See return.
 *	\return The sum of the components of the \ref Vector_i\*.
 */
int sum_Vector_i
	(struct Vector_i* a ///< Standard.
	);

/** \brief See return.
 *	\return `bool` indicating whether the \ref Vector_i\* inputs are equal.
 */
bool check_equal_Vector_i
	(const struct Vector_i*const a, ///< Input 1.
	 const struct Vector_i*const b  ///< Input 2.
	);

/** \brief Comparison function for std::qsort between \ref Vector_i\*\* `a` and `b`.
 *	\return The lexicographical comparison of `a` and `b`.
 *
 *	\note Input Vectors must be have sorted data.
 */
int cmp_Vector_i
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

/// \brief Copy the data of the `src` to that of the `dest` \ref Vector_i\*.
void copy_data_Vector_i_Vector_i
	(const struct Vector_i*const src, ///< Standard.
	 struct Vector_i*const dest       ///< Standard.
	);

// Printing functions *********************************************************************************************** //

/// \brief Print a \ref Vector_i\* to the terminal.
void print_Vector_i
	(const struct Vector_i*const a ///< Standard.
	);

/// \brief Print a \ref const_Vector_i\* to the terminal.
void print_const_Vector_i
	(const struct const_Vector_i*const a ///< Standard.
	);

#endif // DPG__Vector_h__INCLUDED
