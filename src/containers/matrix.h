/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */

#ifndef DPG__matrix_h__INCLUDED
#define DPG__matrix_h__INCLUDED
/** \file
 *  \brief Provides Matrix_\* containers and related functions.
 *
 *  Potentially relevant comments may be found in \ref multiarray.h.
 *
 *  Matrices are 2D Multiarrays.
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

/** \brief Matrix in 'C'ompressed 'S'parse 'R'ow storage format (`double`).
 *
 *  Further details regarding the CSR storage format can be found on [this intel mkl page][mkl_csr].
 *
 *  <!-- References: -->
 *  [mkl_csr]: https://software.intel.com/en-us/mkl-developer-reference-fortran-sparse-blas-csr-matrix-storage-format
 */
struct Matrix_CSR_d {
	ptrdiff_t ext_0, ///< Defined in \ref Matrix_d.
	          ext_1; ///< Defined in \ref Matrix_d.

	ptrdiff_t* row_index; ///< Indices of the entries at the start of each row (`pointerB` on the intel mkl page).
	ptrdiff_t* columns;   ///< Indices of the column corresponding to each entry in the matrix.

	bool owns_data; ///< Defined in \ref Matrix_d.
	double* data;   ///< Defined in \ref Matrix_d.
};

/// \brief `const` version of \ref Matrix_CSR_d.
struct const_Matrix_CSR_d {
	const ptrdiff_t ext_0, ///< Defined in \ref Matrix_CSR_d.
	                ext_1; ///< Defined in \ref Matrix_CSR_d.

	const ptrdiff_t*const row_index; ///< Defined in \ref Matrix_CSR_d.
	const ptrdiff_t*const columns;   ///< Defined in \ref Matrix_CSR_d.

	const bool owns_data;    ///< Defined in \ref Matrix_CSR_d.
	const double*const data; ///< Defined in \ref Matrix_CSR_d.
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

/// \brief Swap the layout.
void swap_layout
	(char*const layout ///< Pointer to the layout variable.
	);

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

/** \brief Get pointer to column of col-major \ref Matrix_i\*.
 *  \return Pointer to the first entry of the column. */
int* get_col_Matrix_i
	(const ptrdiff_t col,     ///< Desired column.
	 const struct Matrix_i* a ///< Matrix.
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

/** \brief Get pointer to a slice (row or column) of a \ref Matrix_d\*.
 *  \return Pointer to the first entry of the slice. */
double* get_slice_Matrix_d
	(const ptrdiff_t slice,   ///< Index of the desired slice.
	 const struct Matrix_d* a ///< Matrix.
	);

/** \brief `const` version of \ref get_slice_Matrix_d.
 *  \return See brief. */
const double* get_slice_const_Matrix_d
	(const ptrdiff_t slice,         ///< Defined for \ref get_slice_Matrix_d.
	 const struct const_Matrix_d* a ///< Defined for \ref get_slice_Matrix_d.
	);

/** \brief Get pointer to row of row-major \ref Matrix_i\*.
 *  \return Pointer to the first entry of the row. */
int* get_row_Matrix_i
	(const ptrdiff_t row,     ///< Desired row.
	 const struct Matrix_i* a ///< Matrix.
	);

/** \brief Get pointer to row of row-major \ref const_Matrix_i\*.
 *  \return Pointer to the first entry of the row. */
const int* get_row_const_Matrix_i
	(const ptrdiff_t row,           ///< Desired row.
	 const struct const_Matrix_i* a ///< Matrix.
	);

/** \brief Get value of the (row,col) entry of a \ref Matrix_i\*.
 *  \return See brief. */
int get_val_Matrix_i
	(const ptrdiff_t row,          ///< The row.
	 const ptrdiff_t col,          ///< The column.
	 const struct Matrix_i*const a ///< Standard.
	);

/** \brief Get value of the (row,col) entry of a \ref const_Matrix_i\*.
 *  \return See brief. */
int get_val_const_Matrix_i
	(const ptrdiff_t row,                ///< The row.
	 const ptrdiff_t col,                ///< The column.
	 const struct const_Matrix_i*const a ///< Standard.
	);

/// \brief Set the values of the destination row to that of the source data.
void set_row_Matrix_d
	(const ptrdiff_t row,        ///< The destination row.
	 struct Matrix_d* dest,      ///< The destination Matrix.
	 const double*const data_src ///< The source data.
	);

/// \brief Set the values of the destination row to that of the source data.
void set_row_Matrix_i
	(const ptrdiff_t row,     ///< The destination row.
	 struct Matrix_i* dest,   ///< The destination Matrix.
	 const int*const data_src ///< The source data.
	);

/// \brief Set the values of the destination column to that of the source data.
void set_col_Matrix_i
	(const ptrdiff_t col,     ///< The destination row.
	 struct Matrix_i* dest,   ///< The destination Matrix.
	 const int*const data_src ///< The source data.
	);

/// \brief Set the values of the destination column to the input value.
void set_col_to_val_Matrix_i
	(const ptrdiff_t col,   ///< The destination row.
	 struct Matrix_i* dest, ///< The destination Matrix.
	 const int data_src     ///< The source data.
	);

/// \brief Set all data entries to the input value.
void set_to_value_Matrix_d
	(struct Matrix_d*const a, ///< Standard.
	 const double val         ///< The value.
	);

/// \brief Set all data entries to the input value.
void set_to_value_Matrix_i
	(struct Matrix_i*const a, ///< Standard.
	 const int val            ///< The value.
	);

/** \brief Compute the opposite layout.
 *  \return See brief. */
char compute_opposite_layout
	(const char layout_i ///< The input layout.
	);

/** \brief See return.
 *  \return The index of a Matrix corresponding to the given row/column input. */
ptrdiff_t compute_index_Matrix
	(const ptrdiff_t i,     ///< The index in the `ext_0` direction.
	 const ptrdiff_t j,     ///< The index in the `ext_1` direction.
	 const ptrdiff_t ext_0, ///< Defined in \ref Matrix_d.
	 const ptrdiff_t ext_1, ///< Defined in \ref Matrix_d.
	 const char layout      ///< Defined in \ref Matrix_d.
	);

#endif // DPG__matrix_h__INCLUDED
