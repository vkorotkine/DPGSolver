// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__containers_h__INCLUDED
#define DPG__containers_h__INCLUDED
/**	\file
 *	\brief Provides standard containers and related functions.
 *
 *	\section s1 Functions
 *
 *	\subsection s1_1 Naming Convention
 *
 *	Function names are chosen according to the following template: `constructor_{0}_(1)_{2}_{3}_{4}` where elements in
 *	curly braces {} are required and those in round brackets () are optional.
 *		- {0} : Type of constructor
 *			- Options: move
 *		- (1) : Optional `const` specifier
 *		- {2} : Type of container to be returned
 *			- Options: Multiarray_d, Matrix_d
 *		- {3} : Level of dereferencing of the returned container object
 *		- (4) : Type of input from which the container is constructed
 *
 *	\subsection s1_2 Variadic Arguments
 *
 *	In the interest of greater generic programming, variadic functions are used for the constructors such that a
 *	variable number of `extent` values may be passed for variable order Multiarrays. This is similar to the
 *	implementation of printf. A detailed explanation of this procedure can be found in this
 *	[Wikipedia article on stdarg.h][stdarg.h].
 *
 *	\section s2 General
 *
 *	The Multiarray struct is intended to be used as a higher-dimensional matrix where move constructors are used to
 *	form matrix structs for appropriate sub-blocks. As the data is stored contiguously in memory, the multi-array may
 *	also be acted on over multiple dimensions at once.
 *
 *	Constructors for const_Multiarray_* structs are not provided directly. To assign to these variables, an lvalue cast
 *	to the appropriate non-const type should be used.
 *
 *	\subsection s2_1 Example
 *
 *	To assign to `const struct const_Multiarray_d*const var1`, use the relevant constructor to set a variable:
 *	`struct Multiarray_d* var2`, then cast the lvalue `*(struct Multiarray_d**)&var1 = var2`.
 *
 *	\warning This procedure exhibits undefined behaviour (relating to changing a const-qualified type) unless the memory
 *	for the struct is **dynamically allocated**. See [this][SO_dyn_const_struct] SO answer for detailed explanation of
 *	procedures adopted here.
 *
 *	<!-- References: -->
 *	[SO_dyn_const_struct]: https://stackoverflow.com/questions/2219001/how-to-initialize-const-members-of-structs-on-the-heap
 *	[stdarg.h]: https://en.wikipedia.org/wiki/Stdarg.h
 */

#include <stddef.h>
#include <stdbool.h>
#include <stdarg.h>

struct Multiarray_d;
struct const_Multiarray_d;

struct Matrix_d;
struct const_Matrix_d;

// Multiarray_* ***************************************************************************************************** //

/// \brief Supports dense multi-dimensional (double) arrays.
struct Multiarray_d {
	char layout; ///< The layout may be 'R'ow or 'C'olumn major.

	size_t  order;   ///< Number of dimensions.
	size_t* extents; ///< Size of arrays in each dimension.

	bool    owns_data; /**< Flag for whether the data should be freed in the destructor. This would be false if a move
	                        constructor was used. */
	double* data;      ///< The data.
};

/// \brief \c const version of \ref Multiarray_d.
struct const_Multiarray_d {
	const char layout;

	const size_t       order;
	const size_t*const extents;

	const bool         owns_data;
	const double*const data;
};

struct Multiarray_d* constructor_move_Multiarray_d_1_d (const char layout, double*const data, const size_t order, ...);

void destructor_Multiarray_d_1 (struct Multiarray_d* A);


/// \brief Set `extents` for a `Multiarray_*`.
size_t* set_extents
	(const size_t order, ///< Defined in \ref Multiarray_d.
	 va_list ap          ///< List of variadic arguments.
	);

/// \brief `size` is the product of the `extents`.
size_t compute_size
	(const size_t order,
	 const size_t *const extents
	);

// Matrix_* ********************************************************************************************************* //

/// \brief 2-dimensional version of \ref Multiarray_d.
struct Matrix_d {
	char layout;

	size_t extents[2];

	bool    owns_data;
	double* data;
};

/// \brief const version of \ref Matrix_d.
struct const_Matrix_d {
	const char layout;

	const size_t extents[2];

	const bool         owns_data;
	const double*const data;
};

/// \brief Constructs an empty \ref Matrix_d.
struct Matrix_d* constructor_empty_Matrix_d
	(const char layout,   ///< Defined in \ref Matrix_d.
	 const size_t n_rows, ///< Value of extents[0].
	 const size_t n_cols  ///< Value of extents[1].
	);

/** \brief Get pointer to row of row-major \ref Matrix_d.
 *	\return Pointer to the first entry of the row.
 */
double* get_row_Matrix_d
	(const size_t row,        ///< Desired row.
	 const struct Matrix_d* a ///< Matrix.
	);

/// \brief Print a \ref Matrix_d to the terminal.
void print_Matrix_d
	(const struct Matrix_d*const a /// Standard.
	);

// Vector_* ********************************************************************************************************* //

/// \brief 1-dimensional Multiarray.
struct Vector_ui {
	size_t extents[1];

	bool    owns_data;
	unsigned int* data;
};

/// \brief Constructs an empty \ref Vector_ui.
struct Vector_ui* constructor_empty_Vector_ui
	(const size_t n_rows ///< The value of extents[0].
	);

/// \brief Print a \ref Vector_ui to the terminal.
void print_Vector_ui
	(const struct Vector_ui*const a ///< Standard.
	);

#endif // DPG__containers_h__INCLUDED
