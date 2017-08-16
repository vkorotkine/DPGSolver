// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Multiarray_h__INCLUDED
#define DPG__Multiarray_h__INCLUDED
/**	\file
 *	\brief Provides Multiarray_\* containers and related functions.
 *
 *	\section s1_Multi General
 *
 *	\subsection s11_Multi Standard Data Types
 *
 *	For standard data types, the Multiarray container is intended to be used as a higher-dimensional matrix where move
 *	constructors are used to form matrix containers for appropriate sub-blocks. As the data is stored contiguously in
 *	memory, the Multiarray may also be acted on over multiple dimensions at once.
 *
 *	\subsection s12_Multi Defined Types
 *
 *	Two specializations of the Multiarray exist: Matrix (2D Multiarray) and Vector (1D Multiarray). These containers are
 *	used when the data is most intuitively considered to be of the given form (e.g. Mathematical operators are
 *	Matrices).
 *
 *	Further, functionality is provided for storing Multiarrays of these specialized types (i.e. Multiarrays of Matrices
 *	and Vectors). When comparing with multiply dereferenced specialized containers (e.g. struct Matrix**), this results
 *	in the fundamental advantage of containers carrying around their size information, easing the burden on the
 *	developer. This implies that containers with multiple levels of dereferencing should never be used.
 *
 *	In order to avoid unintentionally overwriting data which should be constant, `const` versions of the containers
 *	are also provided where relevant. A complication which arises as a result of declaring objects `const` is that it is
 *	not possible to define them! To overcome this difficulty, lvalue casts are used to set the data.
 *
 *	\warning This procedure exhibits undefined behaviour (relating to changing a `const`-qualified type) unless the
 *	memory for the struct is **dynamically allocated**. See [this SO answer][SO_dyn_const_struct] for a detailed
 *	explanation of the problem and possible approaches used to overcome it.
 *
 *	\warning It is **required** that `const` and non-`const` versions of the containers have identical memory layout
 *	         such that casts can be used to convert between them.
 *
 *	Noting that containers should be dynamically allocated, we have the implication that containers with no
 *	dereferencing should never be used. Taken together with the upper limitation on the level of dereferencing, it is
 *	thus required that **containers have exactly one level of dereferencing**. The only place in which an additional
 *	level of dereferencing is permitted is as part of the constructor for a Multiarray.
 *
 *	\section s2_Multi Functions
 *
 *	\subsection s21_Multi Naming Convention
 *
 *	Function names are chosen according to the following template: `constructor_{0}_(1)_{2}_(3)_(4)` where elements in
 *	curly braces {} are required and those in round brackets () are optional.
 *		- {0} : Type of constructor
 *			- default: reserve memory only for the container itself (not for the data).
 *			- empty:   reserve memory storage for the data but do not set it.
 *			- copy:    reserve memory storage for the data and set it to a copy of the input data.
 *			- move:    move data to the container being constructed.
 *		- (1) : Optional `const` specifier
 *		- {2} : Type of container to be returned
 *			- Multiarray_\*: d, Vector_ui
 *		- (3) : Level of dereferencing if not equal to 1.
 *		- (4) : Type of input from which the container is constructed
 *
 *	\subsection s22_Multi Constructors for `const` Containers
 *
 *	The `const_constructor_move_...` functions are used to *define* the `const` equivalent of the container. The
 *	functions move **all** container members through lvalue and rvalue casts.
 *
 *	\note This means that a single destructor call should be made for both the `src` and `dest` variables.
 *
 *	\subsection s23_Multi Variadic Arguments
 *
 *	In the interest of greater generic programming, variadic functions are used for the constructors such that a
 *	variable number of `extent` values may be passed for variable order Multiarrays. This is similar to the
 *	implementation of printf. A detailed explanation of this procedure can be found in this
 *	[Wikipedia article on stdarg.h][stdarg.h].
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

/// \brief Multiarray (`double`).
struct Multiarray_d {
	char layout; ///< The layout may be 'R'ow or 'C'olumn major.

	size_t  order;   ///< Number of dimensions.
	size_t* extents; ///< Size of arrays in each dimension.

	bool    owns_data; /**< Flag for whether the data should be freed in the destructor. This would be false if a move
	                        constructor was used. */
	double* data;      ///< The data.
};

/// \brief Multiarray (`const double`).
struct const_Multiarray_d {
	const char layout;

	const size_t       order;
	const size_t*const extents;

	const bool         owns_data;
	const double*const data;
};

/// \brief Multiarray (`Vector_ui*`).
struct Multiarray_Vector_ui {
	char layout;

	size_t  order;
	size_t* extents;

	bool    owns_data;

	struct Vector_ui** data;
};

/// \brief Multiarray (`const Vector_ui*`).
struct const_Multiarray_Vector_ui {
	const char layout;

	const size_t  order;
	const size_t*const  extents;

	const bool    owns_data;

	const struct const_Vector_ui*const*const data;
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Move constructor for a \ref Multiarray_d\* from a `double*`.
struct Multiarray_d* constructor_move_Multiarray_d_d
	(const char layout,  ///< Defined in \ref Multiarray_d.
	 double*const data,  ///< Defined in \ref Multiarray_d.
	 const size_t order, ///< Defined in \ref Multiarray_d.
	 ...                 ///< Variadic arguments holding the extents.
	);

/** \brief Constructs an empty \ref Multiarray_Vector_ui\*.
 *	\note The layout is set to row-major by default as the data cannot be used directly as for the standard datatypes.
 */
struct Multiarray_Vector_ui* constructor_empty_Multiarray_Vector_ui
	(const size_t order, ///< Defined in \ref Multiarray_d.
	 ...                 ///< Variadic arguments.
	);

/// \brief Constructs a \ref Multiarray_Vector_ui\* and sets the values of its \ref Vector_ui\* components.
struct Multiarray_Vector_ui* constructor_copy_Multiarray_Vector_ui
	(const unsigned int* data_V,     ///< Defined in \ref set_Multiarray_Vector_ui_ui.
	 const unsigned int*const ext_V, ///< Defined in \ref set_Multiarray_Vector_ui_ui.
	 const size_t order,             ///< Defined in \ref Multiarray_d.
	 ...                             ///< Variadic arguments holding the extents of the Multiarray.
	);

/// \brief Move constructor for a `const` \ref const_Multiarray_Vector_ui `*const`.
void const_constructor_move_Multiarray_Vector_ui
	(const struct const_Multiarray_Vector_ui*const* dest, ///< Destination.
	 struct Multiarray_Vector_ui* src                     ///< Source.
	);

/// \brief Destructs a \ref Multiarray_d\*.
void destructor_Multiarray_d
	(struct Multiarray_d* a ///< Standard.
	);

/// \brief Destructs a \ref Multiarray_Vector_ui\*.
void destructor_Multiarray_Vector_ui
	(struct Multiarray_Vector_ui* a ///< Standard.
	);

// Helper functions ************************************************************************************************* //

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

/// \brief Set the values of the \ref Multiarray_Vector_ui based on the input `unsigned int*` data.
void set_Multiarray_Vector_ui_ui
	(struct Multiarray_Vector_ui* a, ///< Standard.
	 const unsigned int*data_V,      ///< Input data for the Vectors.
	 const unsigned int*const ext_V  ///< Input extent[0] for the Vectors.
	);

/** \brief Sort the data of the \ref Multiarray_Vector_ui\*.
 *	\return Optionally return indices.
 */
struct Vector_ui* sort_Multiarray_Vector_ui
	(struct Multiarray_Vector_ui* a, ///< Standard.
	 const bool return_indices       ///< Flag for whether the indices should also be returned.
	);

// Printing functions *********************************************************************************************** //

/// \brief Print a \ref Multiarray_Vector_ui\* to the terminal.
void print_Multiarray_Vector_ui
	(const struct Multiarray_Vector_ui*const a ///< Standard.
	);

/// \brief Print a \ref const_Multiarray_Vector_ui\* to the terminal.
void print_const_Multiarray_Vector_ui
	(const struct const_Multiarray_Vector_ui*const a ///< Standard.
	);

#endif // DPG__Multiarray_h__INCLUDED
