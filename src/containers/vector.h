// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Vector_h__INCLUDED
#define DPG__Vector_h__INCLUDED
/**	\file
 *	\brief Provides Vector_\* containers and related functions.
 *
 *	Potentially relevant comments may be found in \ref multiarray.h.
 *
 *	Vectors are 1D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>

struct const_Matrix_d;

/// \brief Vector (`int`).
struct Vector_i {
	ptrdiff_t ext_0; ///< Defined in \ref Matrix_d.

	bool owns_data; ///< Defined in \ref Multiarray_d.
	int* data;      ///< Defined in \ref Multiarray_d.
};

/// \brief Vector (`const int`).
struct const_Vector_i {
	const ptrdiff_t ext_0; ///< Defined in \ref Vector_i.

	const bool owns_data; ///< Defined in \ref Vector_i.
	const int*const data; ///< Defined in \ref Vector_i.
};

/// \brief Vector (`double`).
struct Vector_d {
	ptrdiff_t ext_0; ///< Defined in \ref Matrix_d.

	bool owns_data; ///< Defined in \ref Multiarray_d.
	double* data;   ///< Defined in \ref Multiarray_d.
};

/// \brief Vector (`const int`).
struct const_Vector_d {
	const ptrdiff_t ext_0; ///< Defined in \ref Vector_i.

	const bool owns_data;    ///< Defined in \ref Vector_i.
	const double*const data; ///< Defined in \ref Vector_i.
};

// Constructor/Destructor functions ********************************************************************************* //

/** \brief Constructs a default \ref Vector_i\*\*.
 *	\return Standard. */
struct Vector_i** constructor_default_Vector_i_2
	(const ptrdiff_t n_dest ///< The number of \ref Vector_i\* components.
	);

/** \brief Constructs an empty \ref Vector_i\*.
 *	\return Standard. */
struct Vector_i* constructor_empty_Vector_i
	(const ptrdiff_t ext_0 ///< The value of ext_0.
	);

/** \brief Copy constructor for a \ref Vector_i\* from a `Vector_i*`.
 *	\return Standard. */
struct Vector_i* constructor_copy_Vector_i
	(const struct Vector_i*const src ///< The source data.
	);

/** \brief Copy constructor for a \ref Vector_i\* from a `const int*`.
 *	\return Standard. */
struct Vector_i* constructor_copy_Vector_i_i
	(const ptrdiff_t ext_0,   ///< The value of ext_0.
	 const int*const data_src ///< The source data.
	);

/** \brief Move constructor for a \ref Vector_i\* from a `int*`.
 *	\return Standard. */
struct Vector_i* constructor_move_Vector_i_i
	(const ptrdiff_t ext_0, ///< The value of ext_0.
	 const bool owns_data,  ///< Standard.
	 int*const data         ///< Standard.
	);

/** \brief Move constructor for a \ref const_Vector_i\* from a `const int*`.
 *	\return Standard. */
struct const_Vector_i* constructor_move_const_Vector_i_i
	(const ptrdiff_t ext_0, ///< The value of ext_0.
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

/** \brief Constructs an empty \ref Vector_d\*.
 *	\return Standard. */
struct Vector_d* constructor_empty_Vector_d
	(const ptrdiff_t ext_0 ///< The value of ext_0.
	);

/** \brief Constructs a \ref Vector_d\* as the sum of the rows/columns of the input Matrix in the specified direction.
 *	\return Standard. */
struct Vector_d* constructor_sum_Vector_d_const_Matrix_d
	(const char sum_dir,                   ///< The direction in which to sum the entries. Options: 'R'ow, 'C'olumn.
	 const struct const_Matrix_d*const src ///< The source matrix.
	);

/// \brief Destructs a \ref Vector_d\*.
void destructor_Vector_d
	(struct Vector_d* a ///< Standard.
	);

// Helper functions ************************************************************************************************* //

/** \brief Reorder a \ref Vector_i based on the provided ordering.
 *	\warning This is not currently done in place.
 */
void reorder_Vector_i
	(struct Vector_i*const a, ///< Standard.
	 const int*const ordering ///< The ordering.
	);

/// \brief Resize a \ref Vector_i\*.
void resize_Vector_i
	(struct Vector_i*const a, ///< Standard.
	 const ptrdiff_t ext_0    ///< New value for ext_0.
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

/**	\brief See return.
 *	\return `bool` indicating whether the \ref Vector_i\* inputs are equal. */
bool check_equal_Vector_i
	(const struct Vector_i*const a, ///< Input 1.
	 const struct Vector_i*const b  ///< Input 2.
	);

/**	\brief See return.
 *	\return `bool` indicating whether the \ref Vector_i\* and `int*` inputs are equal. */
bool check_equal_Vector_i_i
	(const struct Vector_i*const a, ///< The vector input.
	 const int* data_b              ///< The data input.
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

/** \brief Push a value to the back of a \ref Vector_i\*.
 *	The vector may optionally be sorted and/or accept only unique values. */
void push_back_Vector_i
	(struct Vector_i*const src, ///< The source vector.
	 const int val,             ///< The value to be added.
	 const bool sorted,         ///< Flag for whether the vector should be sorted.
	 const bool unique          ///< Flag for whether the value should only be added if it is unique.
	);

/** \brief Check if the input value is present in the optionally sorted src \ref Vector_i\*.
 *	\return `true` if found. */
bool find_val_Vector_i
	(const struct const_Vector_i*const src, ///< The input vector.
	 const int val,                         ///< The value to find.
	 const bool sorted                      ///< Flag indicating whether the vector is sorted.
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

/// \brief Print a \ref Vector_d\* to the terminal displaying entries below the tolerance as 0.0.
void print_Vector_d
	(const struct Vector_d*const a, ///< Standard.
	 const double tol               ///< The tolerance.
	);

#endif // DPG__Vector_h__INCLUDED
