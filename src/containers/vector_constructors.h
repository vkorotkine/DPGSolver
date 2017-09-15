// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__vector_constructors_h__INCLUDED
#define DPG__vector_constructors_h__INCLUDED
/** \file
 *  \brief Provides Vector_\* constructors and destructors.
 */

#include <stddef.h>
#include <stdbool.h>

struct Matrix_d;
struct Multiarray_d;
struct const_Vector_d;
struct const_Matrix_d;
struct const_Multiarray_d;

// Default constructors ********************************************************************************************* //

/** \brief Constructs a default \ref Vector_d\*.
 *  \return Standard. */
struct Vector_d* constructor_default_Vector_d ();

/** \brief Constructor for a default \ref const_Vector_d\*.
 *  \return Standard. */
const struct const_Vector_d* constructor_default_const_Vector_d ();

/** \brief Constructs a default \ref Vector_i\*.
 *  \return Standard. */
struct Vector_i* constructor_default_Vector_i ();

/** \brief Constructs a default \ref Vector_i\*\*.
 *  \return Standard. */
struct Vector_i** constructor_default_Vector_i_2
	(const ptrdiff_t n_dest ///< The number of \ref Vector_i\* components.
	);

// Empty constructors *********************************************************************************************** //

/** \brief Constructs an empty \ref Vector_d\*.
 *  \return Standard. */
struct Vector_d* constructor_empty_Vector_d
	(const ptrdiff_t ext_0 ///< Defined in \ref Vector_d.
	);

/** \brief Constructs an empty \ref Vector_i\*.
 *  \return Standard. */
struct Vector_i* constructor_empty_Vector_i
	(const ptrdiff_t ext_0 ///< The value of ext_0.
	);

// Copy constructors ************************************************************************************************ //

/** \brief Copy constructor for a \ref Vector_i\* from a `Vector_i*`.
 *  \return Standard. */
struct Vector_i* constructor_copy_Vector_i
	(const struct Vector_i*const src ///< The source data.
	);

/** \brief Copy constructor for a \ref Vector_d\* from a `const double*`.
 *  \return Standard. */
struct Vector_d* constructor_copy_Vector_d_d
	(const ptrdiff_t ext_0,      ///< The value of ext_0.
	 const double*const data_src ///< The source data.
	);

/** \brief Copy constructor for a \ref Vector_i\* from a `const int*`.
 *  \return Standard. */
struct Vector_i* constructor_copy_Vector_i_i
	(const ptrdiff_t ext_0,   ///< The value of ext_0.
	 const int*const data_src ///< The source data.
	);

// Move constructors ************************************************************************************************ //

/** \brief Move constructor for a \ref Vector_i\* from a `int*`.
 *  \return Standard. */
struct Vector_i* constructor_move_Vector_i_i
	(const ptrdiff_t ext_0, ///< The value of ext_0.
	 const bool owns_data,  ///< Standard.
	 int*const data         ///< Standard.
	);

/** \brief Move constructor for a \ref Vector_d\* from a `double*`.
 *  \return Standard. */
struct Vector_d* constructor_move_Vector_d_d
	(const ptrdiff_t ext_0, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 double*const data      ///< Standard.
	);

/** \brief Move constructor for a \ref const_Vector_i\* from a `const int*`.
 *  \return Standard. */
struct const_Vector_i* constructor_move_const_Vector_i_i
	(const ptrdiff_t ext_0, ///< The value of ext_0.
	 const bool owns_data,  ///< Standard.
	 const int*const data   ///< Standard.
	);

/** \brief Move constructor for a \ref Vector_d\* from a \ref Matrix_d\*.
 *  \return Standard. */
struct Vector_d* constructor_move_Vector_d_Matrix_d
	(struct Matrix_d* src ///< The source matrix.
	);

/// \brief Move constructor for a `const` \ref const_Vector_d `*const`.
void const_constructor_move_Vector_d
	(const struct const_Vector_d*const* dest, ///< Destination.
	 struct Vector_d* src                     ///< Source.
	);

/// \brief Move constructor for a `const` \ref const_Vector_i `*const`.
void const_constructor_move_Vector_i
	(const struct const_Vector_i*const* dest, ///< Destination.
	 struct Vector_i* src                     ///< Source.
	);

// Set constructors ************************************************************************************************* //

/** \brief Constructor for a \ref Vector_d\* from a sub range of a \ref Multiarray_d\*.
 *  \return Standard. */
struct Vector_d* constructor_set_Vector_d_Multiarray_d
	(struct Multiarray_d* src,         ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/** \brief `const` version of \ref constructor_set_Vector_d_Multiarray_d.
 *  \return Standard. */
const struct const_Vector_d* constructor_set_const_Vector_d_Multiarray_d
	(const struct const_Multiarray_d* src, ///< Defined for \ref constructor_set_Vector_d_Multiarray_d.
	 const ptrdiff_t*const sub_indices     ///< Defined for \ref constructor_set_Vector_d_Multiarray_d.
	);

// Special constructors ********************************************************************************************* //

/** \brief Constructs a \ref Vector_d\* as the sum of the rows/columns of the input Matrix in the specified direction.
 *  \return Standard. */
struct Vector_d* constructor_sum_Vector_d_const_Matrix_d
	(const char sum_dir,                   ///< The direction in which to sum the entries. Options: 'R'ow, 'C'olumn.
	 const struct const_Matrix_d*const src ///< The source matrix.
	);

/** \brief `const` version of \ref constructor_sum_Vector_d_const_Matrix_d.
 *  \return Standard. */
const struct const_Vector_d* constructor_sum_const_Vector_d_const_Matrix_d
	(const char sum_dir,                   ///< Defined for \ref constructor_sum_Vector_d_const_Matrix_d.
	 const struct const_Matrix_d*const src ///< Defined for \ref constructor_sum_Vector_d_const_Matrix_d.
	);

/** \brief Constructor for a \ref Vector_d\* from a matrix-vector multiplication.
 *  \return Standard. */
struct Vector_d* constructor_mv_Vector_d
	(const char trans_a_i,                ///< Defined for \ref mv_d.
	 const double alpha,                  ///< Defined for \ref mv_d.
	 const double beta,                   ///< Defined for \ref mv_d.
	 const struct const_Matrix_d*const a, ///< Defined for \ref mv_d.
	 const struct const_Vector_d*const b  ///< Defined for \ref mv_d.
	);

/** \brief `const` version of \ref constructor_mv_Vector_d.
 *  \return Standard. */
const struct const_Vector_d* constructor_mv_const_Vector_d
	(const char trans_a_i,                ///< Defined for \ref mv_d.
	 const double alpha,                  ///< Defined for \ref mv_d.
	 const double beta,                   ///< Defined for \ref mv_d.
	 const struct const_Matrix_d*const a, ///< Defined for \ref mv_d.
	 const struct const_Vector_d*const b  ///< Defined for \ref mv_d.
	);

/** \brief Constructor for a \ref Vector_d\* from the solution of a linear system using [LAPACKE_dsgesv][dsgesv].
 *  \return Standard.
 *
 *  <!-- References: -->
 *  [dsgesv]: https://software.intel.com/en-us/mkl-developer-reference-c-gesv
 */
struct Vector_d* constructor_sgesv_Vector_d
	(struct Matrix_d* A_i, ///< The LHS input matrix A.
	 struct Vector_d* B_i  ///< The RHS input vector B.
	);

/** \brief `const` version of \ref constructor_sgesv_Vector_d.
 *  \return Standard. */
const struct const_Vector_d* constructor_sgesv_const_Vector_d
	(const struct const_Matrix_d* A_i, ///< The LHS input matrix A.
	 const struct const_Vector_d* B_i  ///< The RHS input vector B.
	);

/// \brief Set a \ref Vector_d\* from a sub range of a \ref Matrix_d\*.
void set_Vector_from_Matrix_d
	(struct Vector_d* dest,            ///< The destination.
	 struct Matrix_d* src,             ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Vector_from_Matrix_d.
void set_const_Vector_from_Matrix_d
	(const struct const_Vector_d* dest, ///< Defined for \ref set_Vector_from_Matrix_d.
	 const struct const_Matrix_d* src,  ///< Defined for \ref set_Vector_from_Matrix_d.
	 const ptrdiff_t*const sub_indices  ///< Defined for \ref set_Vector_from_Matrix_d.
	);

/// \brief Set a \ref Vector_d\* from a sub range of a \ref Multiarray_d\*.
void set_Vector_from_Multiarray_d
	(struct Vector_d* dest,            ///< The destination.
	 struct Multiarray_d* src,         ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Vector_from_Multiarray_d.
void set_const_Vector_from_Multiarray_d
	(const struct const_Vector_d* dest,    ///< Defined for \ref set_Vector_from_Multiarray_d.
	 const struct const_Multiarray_d* src, ///< Defined for \ref set_Vector_from_Multiarray_d.
	 const ptrdiff_t*const sub_indices     ///< Defined for \ref set_Vector_from_Multiarray_d.
	);

// Destructors ****************************************************************************************************** //

/// \brief Destructs a \ref Vector_d\*.
void destructor_Vector_d
	(struct Vector_d* a ///< Standard.
	);

/// \brief Destructs a \ref const_Vector_d\*.
void destructor_const_Vector_d
	(const struct const_Vector_d* a ///< Standard.
	);

/// \brief Destructs a \ref Vector_i\*.
void destructor_Vector_i
	(struct Vector_i* a ///< Standard.
	);

/// \brief Destructs a \ref const_Vector_i\*.
void destructor_const_Vector_i
	(const struct const_Vector_i* a ///< Standard.
	);

/// \brief Destructs a \ref Vector_i\*\*.
void destructor_Vector_i_2
	(struct Vector_i** a,   ///< Standard.
	 const ptrdiff_t n_src, ///< The number of \ref Vector_i\* components.
	 const bool owns_data   ///< Standard.
	);

#endif // DPG__vector_constructors_h__INCLUDED
