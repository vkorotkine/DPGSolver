// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Matrix_h__INCLUDED
#define DPG__Matrix_h__INCLUDED
/**	\file
 *	\brief Provides Matrix_\* containers and related functions.
 *
 *	Potentially relevant comments may be found in \ref multiarray.h.
 *
 *	Matrices are 2D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>

struct const_Vector_i;

/// \brief Matrix (`double`).
struct Matrix_d {
	char layout; ///< Defined in \ref Multiarray_d.

	ptrdiff_t ext_0, ///< Size of array in the 1st dimension.
	          ext_1; ///< Size of array in the 2nd dimension.

	bool owns_data; ///< Defined in \ref Multiarray_d.
	double* data;   ///< Defined in \ref Multiarray_d.
};

/// \brief Matrix (`const double`).
struct const_Matrix_d {
	const char layout; ///< Defined in \ref Matrix_d.

	const ptrdiff_t ext_0, ///< Defined in \ref Matrix_d.
	                ext_1; ///< Defined in \ref Matrix_d.

	const bool owns_data;    ///< Defined in \ref Matrix_d.
	const double*const data; ///< Defined in \ref Matrix_d.
};

/// \brief Matrix (`int`).
struct Matrix_i {
	char layout; ///< Defined in \ref Matrix_d.

	ptrdiff_t ext_0, ///< Defined in \ref Matrix_d.
	          ext_1; ///< Defined in \ref Matrix_d.

	bool owns_data; ///< Defined in \ref Matrix_d.
	int* data;      ///< Defined in \ref Matrix_d.
};

/// \brief Matrix (`const int`).
struct const_Matrix_i {
	const char layout; ///< Defined in \ref Matrix_d.

	const ptrdiff_t ext_0, ///< Defined in \ref Matrix_d.
	                ext_1; ///< Defined in \ref Matrix_d.

	const bool owns_data; ///< Defined in \ref Matrix_d.
	const int*const data; ///< Defined in \ref Matrix_d.
};

// Constructor/Destructor functions ********************************************************************************* //

/**	\brief Constructor for a default \ref Matrix_d\*.
 *	\return Standard. */
struct Matrix_d* constructor_default_Matrix_d ();

/** \brief Constructs an empty \ref Matrix_d\*.
 *	\return Standard. */
struct Matrix_d* constructor_empty_Matrix_d
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1  ///< Standard.
	);

/** \brief Copy constructor for a \ref Matrix_i\* from a `const int*`.
 *	\return Standard. */
struct Matrix_i* constructor_copy_Matrix_i_i
	(const char layout,       ///< Standard.
	 const ptrdiff_t ext_0,   ///< Standard.
	 const ptrdiff_t ext_1,   ///< Standard.
	 const int*const data_src ///< The source data.
	);

/** \brief Copy constructor for a \ref Matrix_d\* from a `const double*`.
 *	\return Standard. */
struct Matrix_d* constructor_copy_Matrix_d_d
	(const char layout,          ///< Standard.
	 const ptrdiff_t ext_0,      ///< Standard.
	 const ptrdiff_t ext_1,      ///< Standard.
	 const double*const data_src ///< The source data.
	);

/** \brief Copy constructor for a \ref Matrix_d\* from a \ref Matrix_d\*.
 *	\return Standard. */
struct Matrix_d* constructor_copy_Matrix_d
	(struct Matrix_d* src /// The source matrix.
	);

/** \brief Copy constructor for a \ref const_Matrix_d\* from a partial number of rows/columns of another.
 *	\return Standard. */
const struct const_Matrix_d* constructor_copy_extract_const_Matrix_d
	(const struct const_Matrix_d*const src,    ///< The source Matrix.
	 const struct const_Vector_i*const indices ///< The indices of the rows/columns to copy.
	);

/** \brief Move Constructor for a \ref const_Matrix_d\*.
 *	\return Standard. */
struct const_Matrix_d* constructor_move_const_Matrix_d_Matrix_d
	(struct Matrix_d*const src ///< Source Matrix.
	);

/** \brief Copy constructor for a `const` \ref const_Matrix_d\* from a `const` \ref const_Matrix_d\*.
 *	\return Standard. */
void const_constructor_copy_Matrix_d
	(const struct const_Matrix_d*const* dest, ///< Destination.
	 const struct const_Matrix_d*const src    ///< Source.
	);

/// \brief Move Constructor for a `const` \ref const_Matrix_d `*const`.
void const_constructor_move_Matrix_d
	(const struct const_Matrix_d*const* dest, ///< Destination.
	 struct Matrix_d* src                     ///< Source.
	);

/// \brief Destructs a \ref Matrix_d\*.
void destructor_Matrix_d
	(struct Matrix_d* a ///< Standard.
	);

/** \brief Constructs an empty \ref Matrix_i\*.
 *	\return Standard. */
struct Matrix_i* constructor_empty_Matrix_i
	(const char layout,  ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1  ///< Standard.
	);

/** \brief Move Constructor for a \ref const_Matrix_i\*.
 *	\return Standard. */
struct const_Matrix_i* constructor_move_const_Matrix_i_Matrix_i
	(struct Matrix_i*const src ///< Source Matrix.
	);

/** \brief Move Constructor for a `const` \ref const_Matrix_i `*const`.
 *	\return Standard. */
void const_constructor_move_Matrix_i
	(const struct const_Matrix_i*const* dest, ///< Destination.
	 struct Matrix_i* src                     ///< Source.
	);

/// \brief Destructs a \ref Matrix_i\*.
void destructor_Matrix_i
	(struct Matrix_i* a ///< Standard.
	);

// Helper functions ************************************************************************************************* //

/** \brief Get pointer to row of row-major \ref Matrix_d\*.
 *	\return Pointer to the first entry of the row.
 */
double* get_row_Matrix_d
	(const ptrdiff_t row,     ///< Desired row.
	 const struct Matrix_d* a ///< Matrix.
	);

/** \brief Get pointer to row of row-major \ref const_Matrix_d\*.
 *	\return Pointer to the first entry of the row.
 */
const double* get_row_const_Matrix_d
	(const ptrdiff_t row,           ///< Desired row.
	 const struct const_Matrix_d* a ///< Matrix.
	);

/** \brief Get pointer to row of row-major \ref Matrix_i\*.
 *	\return Pointer to the first entry of the row.
 */
int* get_row_Matrix_i
	(const ptrdiff_t row,     ///< Desired row.
	 const struct Matrix_i* a ///< Matrix.
	);

/** \brief Compute the norm of the specified row of the input \ref Matrix_d.
 *	\return See brief. */
double compute_norm_Matrix_d_row
	(const ptrdiff_t row,           ///< The row.
	 const struct Matrix_d*const a, ///< The input matrix.
	 const char*const norm_type     ///< The norm type.
	);

/** \brief Get value of the (row,col) entry of a \ref Matrix_i\*.
 *	\return See brief. */
int get_val_Matrix_i
	(const ptrdiff_t row,          ///< The row.
	 const ptrdiff_t col,          ///< The column.
	 const struct Matrix_i*const a ///< Standard.
	);

/** \brief Get value of the (row,col) entry of a \ref const_Matrix_i\*.
 *	\return See brief. */
int get_val_const_Matrix_i
	(const ptrdiff_t row,                ///< The row.
	 const ptrdiff_t col,                ///< The column.
	 const struct const_Matrix_i*const a ///< Standard.
	);

/// \brief Set the values of the destination row to that of the source data.
void set_row_Matrix_d
	(const ptrdiff_t row,         ///< The destination row.
	 const struct Matrix_d* dest, ///< The destination Matrix.
	 const double*const data_src  ///< The source data.
	);

/// \brief Transpose the \ref Matrix_d\* optionally leaving the values of the extents unchanged in `mem_only = true`.
void transpose_Matrix_d
	(struct Matrix_d* a, ///< Matrix to be transposed.
	 const bool mem_only ///< Flag for whether only the memory should be transposed (with ext_0/ext_1 unchanged).
	);

// Printing functions *********************************************************************************************** //

/// \brief Print a \ref Matrix_d\* to the terminal displaying entries below the tolerance as 0.0.
void print_Matrix_d
	(const struct Matrix_d*const a, ///< Standard.
	 const double tol               ///< The tolerance.
	);

/// \brief Print a \ref const_Matrix_d\* as in \ref print_Matrix_d.
void print_const_Matrix_d
	(const struct const_Matrix_d*const a, ///< Standard.
	 const double tol                     ///< The tolerance.
	);

/// \brief Print a \ref Matrix_i\* to the terminal.
void print_Matrix_i
	(const struct Matrix_i*const a ///< Standard.
	);

/// \brief Print a \ref const_Matrix_i\* to the terminal.
void print_const_Matrix_i
	(const struct const_Matrix_i*const a ///< Standard.
	);


#endif // DPG__Matrix_h__INCLUDED
