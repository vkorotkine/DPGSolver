// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__matrix_h__INCLUDED
#define DPG__matrix_h__INCLUDED
/**	\file
 *	\brief Provides Matrix_\* containers and related functions.
 *
 *	Potentially relevant comments may be found in \ref multiarray.h.
 *
 *	Matrices are 2D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>

#include "matrix_constructors.h"
#include "matrix_math.h"
#include "matrix_print.h"

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

// Interface functions ********************************************************************************************** //

/** \brief Get pointer to row of row-major \ref Matrix_d\*.
 *  \return Pointer to the first entry of the row. */
double* get_row_Matrix_d
	(const ptrdiff_t row,     ///< Desired row.
	 const struct Matrix_d* a ///< Matrix.
	);

/** \brief Get pointer to column of col-major \ref Matrix_d\*.
 *  \return Pointer to the first entry of the column. */
double* get_col_Matrix_d
	(const ptrdiff_t col,     ///< Desired column.
	 const struct Matrix_d* a ///< Matrix.
	);

/** \brief Get pointer to row of row-major \ref const_Matrix_d\*.
 *  \return Pointer to the first entry of the row. */
const double* get_row_const_Matrix_d
	(const ptrdiff_t row,           ///< Desired row.
	 const struct const_Matrix_d* a ///< Matrix.
	);

/** \brief Get pointer to column of col-major \ref const_Matrix_d\*.
 *  \return Pointer to the first entry of the column. */
const double* get_col_const_Matrix_d
	(const ptrdiff_t col,                ///< Desired column.
	 const struct const_Matrix_d*const a ///< Matrix.
	);

/**	\brief Get pointer to row of row-major \ref Matrix_i\*.
 *	\return Pointer to the first entry of the row.
 */
int* get_row_Matrix_i
	(const ptrdiff_t row,     ///< Desired row.
	 const struct Matrix_i* a ///< Matrix.
	);

/**	\brief Get value of the (row,col) entry of a \ref Matrix_i\*.
 *	\return See brief. */
int get_val_Matrix_i
	(const ptrdiff_t row,          ///< The row.
	 const ptrdiff_t col,          ///< The column.
	 const struct Matrix_i*const a ///< Standard.
	);

/**	\brief Get value of the (row,col) entry of a \ref const_Matrix_i\*.
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

#endif // DPG__matrix_h__INCLUDED
