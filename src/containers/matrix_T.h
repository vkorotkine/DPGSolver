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
/** \file
 *  \brief Provides templated Matrix_\* containers and related functions.
 *
 *  Potentially relevant comments may be found in \ref multiarray.h.
 *
 *  Matrices are 2D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>

/// \brief Templated Matrix.
struct Matrix_T {
	char layout; ///< Memory layout. Options: 'R'ow/'C'olumn-major.

	ptrdiff_t ext_0, ///< Size of array in the 1st dimension.
	          ext_1; ///< Size of array in the 2nd dimension.

	bool owns_data; ///< Flag for whether the data should be freed when calling the destructor.
	Type* data;     ///< The data.
};

/// \brief `const` version of \ref Matrix_T.
struct const_Matrix_T { ///\{
	const char layout;

	const ptrdiff_t ext_0,
	                ext_1;

	const bool owns_data;
	const Type*const data;
}; ///\}

/** \brief Templated Matrix in 'C'ompressed 'S'parse 'R'ow storage format.
 *
 *  Further details regarding the CSR storage format can be found on [this intel mkl page][mkl_csr].
 *
 *  <!-- References: -->
 *  [mkl_csr]: https://software.intel.com/en-us/mkl-developer-reference-fortran-sparse-blas-csr-matrix-storage-format
 */
struct Matrix_CSR_T {
	ptrdiff_t ext_0, ///< Defined in \ref Matrix_T.
	          ext_1; ///< Defined in \ref Matrix_T.

	ptrdiff_t* row_index; ///< Indices of the entries at the start of each row (`pointerB` on the intel mkl page).
	ptrdiff_t* columns;   ///< Indices of the column corresponding to each entry in the matrix.

	bool owns_data; ///< Defined in \ref Matrix_T.
	Type* data;     ///< Defined in \ref Matrix_T.
};

/// \brief `const` version of \ref Matrix_CSR_T.
struct const_Matrix_CSR_T { ///\{
	const ptrdiff_t ext_0,
	                ext_1;

	const ptrdiff_t*const row_index;
	const ptrdiff_t*const columns;

	const bool owns_data;
	const Type*const data;
}; ///\}

// Interface functions ********************************************************************************************** //

/** \brief Pointer to value setting functions.
 *
 *  \param dest The destination.
 *  \param src  The source.
 */
typedef void (*set_value_fptr_T)
	(Type*const dest,
	 const Type src
	);

/** \brief Get pointer to row of row-major \ref Matrix_T\*.
 *  \return Pointer to the first entry of the row. */
Type* get_row_Matrix_T
	(const ptrdiff_t row,     ///< Desired row.
	 const struct Matrix_T* a ///< Matrix.
	);

/** \brief Get pointer to row of row-major \ref const_Matrix_T\*.
 *  \return Pointer to the first entry of the row. */
const Type* get_row_const_Matrix_T
	(const ptrdiff_t row,           ///< Desired row.
	 const struct const_Matrix_T* a ///< Matrix.
	);

/** \brief Get pointer to column of col-major \ref Matrix_T\*.
 *  \return Pointer to the first entry of the column. */
Type* get_col_Matrix_T
	(const ptrdiff_t col,     ///< Desired column.
	 const struct Matrix_T* a ///< Matrix.
	);

/** \brief Get pointer to column of col-major \ref const_Matrix_T\*.
 *  \return Pointer to the first entry of the column. */
const Type* get_col_const_Matrix_T
	(const ptrdiff_t col,                ///< Desired column.
	 const struct const_Matrix_T*const a ///< Matrix.
	);

/** \brief Get pointer to a slice (row or column) of a \ref Matrix_T\*.
 *  \return Pointer to the first entry of the slice. */
Type* get_slice_Matrix_T
	(const ptrdiff_t slice,   ///< Index of the desired slice.
	 const struct Matrix_T* a ///< Matrix.
	);

/** \brief `const` version of \ref get_slice_Matrix_T.
 *  \return See brief. */
const Type* get_slice_const_Matrix_T
	(const ptrdiff_t slice,         ///< See brief.
	 const struct const_Matrix_T* a ///< See brief.
	);

/** \brief Get value of the (row,col) entry of a \ref Matrix_T\*.
 *  \return See brief. */
Type get_val_Matrix_T
	(const ptrdiff_t row,          ///< The row.
	 const ptrdiff_t col,          ///< The column.
	 const struct Matrix_T*const a ///< Standard.
	);

/** \brief Get value of the (row,col) entry of a \ref const_Matrix_T\*.
 *  \return See brief. */
Type get_val_const_Matrix_T
	(const ptrdiff_t row,                ///< The row.
	 const ptrdiff_t col,                ///< The column.
	 const struct const_Matrix_T*const a ///< Standard.
	);

/// \brief Set the values of the destination row to that of the source data.
void set_row_Matrix_T
	(const ptrdiff_t row,      ///< The destination row.
	 struct Matrix_T* dest,    ///< The destination Matrix.
	 const Type*const data_src ///< The source data.
	);

/// \brief Set the values of the destination column to that of the source data.
void set_col_Matrix_T
	(const ptrdiff_t col,      ///< The destination row.
	 struct Matrix_T* dest,    ///< The destination Matrix.
	 const Type*const data_src ///< The source data.
	);

/// \brief Set the values of the destination column to the input value.
void set_col_to_val_Matrix_T
	(const ptrdiff_t col,   ///< The destination row.
	 struct Matrix_T* dest, ///< The destination Matrix.
	 const Type data_src    ///< The source data.
	);

/// \brief Set all data entries to the input value.
void set_to_value_Matrix_T
	(struct Matrix_T*const a, ///< Standard.
	 const Type val         ///< The value.
	);

/// \brief Set a sub-block of a \ref Matrix_T to the entries of the input matrix.
void set_block_Matrix_T
	(struct Matrix_T* a,                 ///< The large matrix.
	 const struct const_Matrix_T* a_sub, ///< The matrix holding the values to set in the sub-block.
	 const ptrdiff_t row0,               ///< The index of the first row where the sub-block should be placed.
	 const ptrdiff_t col0,               ///< The index of the first column where the sub-block should be placed.
	 const char set_type                 ///< The type of setting to use. Options: 'i'nsert, 'a'dd.
	);
#if TYPE_RC == TYPE_COMPLEX
/// \brief Version of \ref set_block_Matrix_T with a real input.
void set_block_Matrix_T_R
	(struct Matrix_T* a,                 ///< See brief.
	 const struct const_Matrix_R* a_sub, ///< See brief.
	 const ptrdiff_t row0,               ///< See brief.
	 const ptrdiff_t col0,               ///< See brief.
	 const char set_type                 ///< See brief.
	);

/** \brief Set a sub-block of a real \ref Matrix_T to the imaginary entries of the input complex \ref const_Matrix_T
 *         matrix after after division by \ref CX_STEP.
 *  This function is identical to \ref set_block_Matrix_T excluding the setting portion.
 */
void set_block_Matrix_R_cmplx_step
	(struct Matrix_R* a,                 ///< The large matrix.
	 const struct const_Matrix_C* a_sub, ///< The matrix holding the values to set in the sub-block.
	 const ptrdiff_t row0,               ///< The index of the first row where the sub-block should be placed.
	 const ptrdiff_t col0,               ///< The index of the first column where the sub-block should be placed.
	 const char set_type                 ///< The type of setting to use. Options: 'i'nsert, 'a'dd.
	);
#endif

/// \brief Version of \ref set_value_fptr_T inserting values.
void set_value_insert_T
	(Type*const dest, ///< See brief.
	 const Type src   ///< See brief.
	);

/// \brief Version of \ref set_value_fptr_T adding values.
void set_value_add_T
	(Type*const dest, ///< See brief.
	 const Type src   ///< See brief.
	);

