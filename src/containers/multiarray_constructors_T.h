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
 *  \brief Provides Multiarray_\* constructors and destructors.
 */

#include <stddef.h>
#include <stdbool.h>

struct Vector_T;
struct Matrix_T;
struct Multiarray_R;
struct Multiarray_Vector_T;
struct Multiarray_Matrix_T;
struct const_Vector_i;
struct const_Vector_T;
struct const_Matrix_R;
struct const_Matrix_T;
struct const_Multiarray_R;
struct const_Multiarray_T;
struct const_Multiarray_Vector_T;
struct const_Multiarray_Matrix_T;

// Helper functions ************************************************************************************************* //

/**	\brief Allocated and set the `extents` for a `Multiarray_*`.
 *	\return See brief. */
ptrdiff_t* allocate_and_set_extents
	(const int order,                ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

// Default constructors ********************************************************************************************* //

/** \brief Constructor for a default \ref Multiarray_T\*.
 *  \return Standard. */
struct Multiarray_T* constructor_default_Multiarray_T ();

/** \brief `const` version of \ref constructor_default_Multiarray_T.
 *  \return Standard. */
const struct const_Multiarray_T* constructor_default_const_Multiarray_T ();

/** \brief Constructor for a default \ref Multiarray_Matrix_T\*.
 *  \return Standard. */
struct Multiarray_Matrix_T* constructor_default_Multiarray_Matrix_T ();

/** \brief `const` version of \ref constructor_default_Multiarray_Matrix_T.
 *  \return Standard. */
const struct const_Multiarray_Matrix_T* constructor_default_const_Multiarray_Matrix_T ();

// Empty constructors *********************************************************************************************** //

/** \brief Constructor for an empty \ref Multiarray_T\*.
 *  \return Standard. */
struct Multiarray_T* constructor_empty_Multiarray_T
	(const char layout,              ///< Defined in \ref Multiarray_T.
	 const int order,                ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

/** \brief Constructor for an empty \ref Multiarray_T\* with extents having been previously dynamically allocated.
 *  \return Standard. */
struct Multiarray_T* constructor_empty_Multiarray_T_dyn_extents
	(const char layout,            ///< Defined in \ref Multiarray_T.
	 const int order,              ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const extents ///< Defined in \ref Multiarray_T.
	);

/** \brief Constructor for an empty \ref Multiarray_Vector_T\*.
 *  \return Standard. */
struct Multiarray_Vector_T* constructor_empty_Multiarray_Vector_T
	(const bool alloc_V,             ///< Flag for whether memory should be reserved for the individual Vectors.
	 const int order,                ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

/** \brief `const` version of \ref constructor_empty_Multiarray_Vector_T.
 *  \return Standard. */
const struct const_Multiarray_Vector_T* constructor_empty_const_Multiarray_Vector_T
	(const bool alloc_V,             ///< Defined for \ref constructor_empty_Multiarray_Vector_T.
	 const int order,                ///< Defined for \ref constructor_empty_Multiarray_Vector_T.
	 const ptrdiff_t*const extents_i ///< Defined for \ref constructor_empty_Multiarray_Vector_T.
	);

/** \brief Constructor for an empty \ref Multiarray_Vector_T\* where a \ref const_Vector_T holds the order/extents.
 *  \return Standard. */
const struct const_Multiarray_Vector_T* constructor_empty_const_Multiarray_Vector_T_V
	(const bool alloc_V,                           ///< Defined for \ref constructor_empty_Multiarray_Vector_T.
	 const struct const_Vector_i*const extents_i_V ///< Vector holding the order and extents.
	);

/** \brief Constructor for an empty \ref Multiarray_Vector_T\*.
 *  \return Standard. */
struct Multiarray_Vector_T* constructor_empty_Multiarray_Vector_T
	(const bool alloc_V,             ///< Flag for whether memory should be reserved for the individual Vectors.
	 const int order,                ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

/** \brief `const` version of \ref constructor_empty_Multiarray_Vector_T.
 *  \return Standard. */
const struct const_Multiarray_Vector_T* constructor_empty_const_Multiarray_Vector_T
	(const bool alloc_V,             ///< Defined for \ref constructor_empty_Multiarray_Vector_T.
	 const int order,                ///< Defined for \ref constructor_empty_Multiarray_Vector_T.
	 const ptrdiff_t*const extents_i ///< Defined for \ref constructor_empty_Multiarray_Vector_T.
	);

/** \brief Constructor for an empty \ref Multiarray_Matrix_T\*.
 *  \return Standard. */
struct Multiarray_Matrix_T* constructor_empty_Multiarray_Matrix_T
	(const bool alloc_M,             ///< Flag for whether memory should be reserved for the individual Matrices.
	 const int order,                ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

/** \brief `const` version of \ref constructor_empty_Multiarray_Matrix_T.
 *  \return Standard. */
const struct const_Multiarray_Matrix_T* constructor_empty_const_Multiarray_Matrix_T
	(const bool alloc_M,             ///< Defined for \ref constructor_empty_Multiarray_Matrix_T.
	 const int order,                ///< Defined for \ref constructor_empty_Multiarray_Matrix_T.
	 const ptrdiff_t*const extents_i ///< Defined for \ref constructor_empty_Multiarray_Matrix_T.
	);

/** \brief Constructor for an empty \ref Multiarray_Matrix_T\* with extents input as a \ref const_Vector_T\*.
 *  \return Standard. */
struct Multiarray_Matrix_T* constructor_empty_Multiarray_Matrix_T_V
	(const bool alloc_M,                           ///< Defined for \ref constructor_empty_Multiarray_Matrix_T.
	 const struct const_Vector_i*const extents_i_V ///< The input extents in vector format.
	);

/** \brief `const` version of \ref constructor_empty_Multiarray_Matrix_T_V.
 *  \return Standard. */
const struct const_Multiarray_Matrix_T* constructor_empty_const_Multiarray_Matrix_T_V
	(const bool alloc_M,                           ///< Defined for \ref constructor_empty_Multiarray_Matrix_T.
	 const struct const_Vector_i*const extents_i_V ///< The input extents in vector format.
	);

// Zero constructors ************************************************************************************************ //

/** \brief Same as \ref constructor_empty_Multiarray_T but with data calloc'ed.
 *  \return Standard. */
struct Multiarray_T* constructor_zero_Multiarray_T
	(const char layout,              ///< Defined in \ref Multiarray_T.
	 const int order,                ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

/** \brief Same as \ref constructor_empty_Multiarray_T_dyn_extents but with data calloc'ed.
 *  \return Standard. */
struct Multiarray_T* constructor_zero_Multiarray_T_dyn_extents
	(const char layout,              ///< Defined in \ref Multiarray_T.
	 const int order,                ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const extents_i ///< Defined in \ref Multiarray_T.
	);

// Copy constructors ************************************************************************************************ //

/** \brief Constructs a \ref Multiarray_Vector_T\* and sets the values of its \ref Vector_T\* components from the
 *         intput `Type*` data.
 *  \return Standard. */
struct Multiarray_Vector_T* constructor_copy_Multiarray_Vector_T_T
	(const Type* data_V,              ///< The `Type` data.
	 const int*const ext_V,          /**< The ext_0 values of each \ref Vector_T in the multiarray. Note: using
	                                      `int` and not `ptrdiff_t` as this is sometimes initialized from
	                                      \ref Vector_T::data. */
	 const int order,                ///< Defined in \ref Multiarray_T.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

/** \brief Copy constructor for a \ref Multiarray_T\* from a \ref Multiarray_T\*.
 *  \return Standard. */
struct Multiarray_T* constructor_copy_Multiarray_T
	(struct Multiarray_T* src ///< Source.
	);

/** \brief `const` version of \ref constructor_copy_Multiarray_T.
 *  \return Standard. */
const struct const_Multiarray_T* constructor_copy_const_Multiarray_T
	(const struct const_Multiarray_T*const src ///< Defined for \ref constructor_copy_Multiarray_T.
	);

/** \brief Copy constructor for a `const` \ref const_Multiarray_T\* from a `const` \ref const_Multiarray_T\*.
 *  \return Standard. */
void const_constructor_copy_Multiarray_T
	(const struct const_Multiarray_T*const* dest, ///< Destination.
	 const struct const_Multiarray_T*const src    ///< Source.
	);
#if TYPE_RC == TYPE_COMPLEX
/** \brief Copy constructor for a \ref Multiarray_T\* from a \ref Multiarray_R\*.
 *  \return Standard. */
struct Multiarray_T* constructor_copy_Multiarray_T_Multiarray_R
	(struct Multiarray_R* src ///< Source.
	);

/** \brief `const` version of \ref constructor_copy_Multiarray_T_Multiarray_R.
 *  \return Standard. */
const struct const_Multiarray_T* constructor_copy_const_Multiarray_T_Multiarray_R
	(const struct const_Multiarray_R* src ///< See brief.
	);
#endif
// Move constructors ************************************************************************************************ //

/** \brief Move constructor for a \ref Multiarray_T\* from a `Type*`.
 *  \return Standard. */
struct Multiarray_T* constructor_move_Multiarray_T_T
	(const char layout,               ///< Standard.
	 const int order,                 ///< Standard.
	 const ptrdiff_t*const extents_i, ///< The input extents.
	 const bool owns_data,            ///< Standard.
	 Type*const data                ///< Standard.
	);

/** \brief `const` version of \ref constructor_move_Multiarray_T_T.
 *  \return Standard. */
const struct const_Multiarray_T* constructor_move_const_Multiarray_T_T
	(const char layout,               ///< Defined for \ref constructor_move_Multiarray_T_T.
	 const int order,                 ///< Defined for \ref constructor_move_Multiarray_T_T.
	 const ptrdiff_t*const extents_i, ///< Defined for \ref constructor_move_Multiarray_T_T.
	 const bool owns_data,            ///< Defined for \ref constructor_move_Multiarray_T_T.
	 const Type*const data          ///< Defined for \ref constructor_move_Multiarray_T_T.
	);

/** \brief Move constructor for a \ref Multiarray_Vector_T\* with the input extents having been previously
 *         dynamically allocated.
 *  \return See brief. */
struct Multiarray_Vector_T* constructor_move_Multiarray_Vector_T_dyn_extents
	(const int order,            ///< Standard.
	 ptrdiff_t*const extents,    ///< Standard.
	 const bool owns_data,       ///< Standard.
	 struct Vector_T**const data ///< Standard.
	);

/** \brief Move constructor for a \ref Multiarray_Matrix_T\* with the input extents having been previously
 *         dynamically allocated.
 *  \return See brief. */
struct Multiarray_Matrix_T* constructor_move_Multiarray_Matrix_T_dyn_extents
	(const int order,            ///< Standard.
	 ptrdiff_t*const extents,    ///< Standard.
	 const bool owns_data,       ///< Standard.
	 struct Matrix_T**const data ///< Standard.
	);

/** \brief Move constructor for a \ref Multiarray_T\* with the input extents having been previously dynamically
 *         allocated.
 *  \return See brief. */
struct Multiarray_T* constructor_move_Multiarray_T_dyn_extents
	(const char layout,       ///< Standard.
	 const int order,         ///< Standard.
	 ptrdiff_t*const extents, ///< Standard.
	 const bool owns_data,    ///< Standard.
	 Type*const data           ///< Standard.
	);

/** `const` version of \ref constructor_move_Multiarray_T_dyn_extents.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_move_const_Multiarray_T_dyn_extents
	(const char layout,       ///< Defined for \ref constructor_move_Multiarray_T_dyn_extents.
	 const int order,         ///< Defined for \ref constructor_move_Multiarray_T_dyn_extents.
	 ptrdiff_t*const extents, ///< Defined for \ref constructor_move_Multiarray_T_dyn_extents.
	 const bool owns_data,    ///< Defined for \ref constructor_move_Multiarray_T_dyn_extents.
	 const Type*const data  ///< Defined for \ref constructor_move_Multiarray_T_dyn_extents.
	);

/** \brief Move constructor for a \ref Multiarray_T\* from a \ref Matrix_T\*.
 *  \return Standard. */
struct Multiarray_T* constructor_move_Multiarray_T_Matrix_T
	(struct Matrix_T* src ///< The source matrix.
	);

/** \brief `const` version of \ref constructor_move_Multiarray_T_Matrix_T.
 *  \return Standard. */
const struct const_Multiarray_T* constructor_move_const_Multiarray_T_Matrix_T
	(const struct const_Matrix_T* src ///< Defined for \ref constructor_move_Multiarray_T_Matrix_T.
	);

/// \brief Move Constructor for a `const` \ref const_Multiarray_T `*const`.
void const_constructor_move_Multiarray_T
	(const struct const_Multiarray_T*const* dest, ///< Destination.
	 struct Multiarray_T* src                     ///< Source.
	);

/// \brief `const` version of \ref const_constructor_move_Multiarray_T.
void const_constructor_move_const_Multiarray_T
	(const struct const_Multiarray_T*const* dest, ///< Defined for \ref const_constructor_move_Multiarray_T.
	 const struct const_Multiarray_T* src         ///< Defined for \ref const_constructor_move_Multiarray_T.
	);

/// \brief Move constructor for a `const` \ref const_Multiarray_Vector_T `*const`.
void const_constructor_move_Multiarray_Vector_T
	(const struct const_Multiarray_Vector_T*const* dest, ///< Destination.
	 struct Multiarray_Vector_T* src                     ///< Source.
	);

/// \brief Move constructor for a `const` \ref const_Multiarray_Matrix_T `*const`.
void const_constructor_move_Multiarray_Matrix_T
	(const struct const_Multiarray_Matrix_T*const* dest, ///< Destination.
	 struct Multiarray_Matrix_T* src                     ///< Source.
	);

// Special constructors ********************************************************************************************* //

/** \brief Constructor for a \ref const_Multiarray_T\* **of order 2** from the matrix-vector multiplications of the
 *         matrices of `A` with the vector `b`.
 *  \return See brief.
 */
const struct const_Multiarray_T* constructor_MaM1_V_const_Multiarray_T
	(const char layout,                              ///< The layout of the output Multiarray.
	 const char trans_a,                             ///< Defined for \ref mv_d.
	 const Real alpha,                               ///< Defined for \ref mv_d.
	 const Real beta,                                ///< Defined for \ref mv_d.
	 const struct const_Multiarray_Matrix_T*const A, ///< Defined for \ref mv_d.
	 const struct const_Vector_T* b                  ///< Defined for \ref mv_d.
	);

/// \brief Set a \ref Multiarray_Matrix_T\* from a sub range of a \ref Multiarray_Matrix_T\*.
void set_Multiarray_Matrix_from_Multiarray_Matrix_T
	(struct Multiarray_Matrix_T* dest, ///< The destination.
	 struct Multiarray_Matrix_T* src,  ///< The source.
	 const int order_o,                ///< The order of the output (destination).
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Multiarray_Matrix_from_Multiarray_Matrix_T.
void set_const_Multiarray_Matrix_from_Multiarray_Matrix_T
	(const struct const_Multiarray_Matrix_T* dest, ///< Defined for mutable version.
	 const struct const_Multiarray_Matrix_T* src,  ///< Defined for mutable version.
	 const int order_o,                            ///< Defined for mutable version.
	 const ptrdiff_t*const sub_indices             ///< Defined for mutable version.
	);

/** \brief Constructor for a \ref Multiarray_T\* using a matrix-matrix multiplication, interpreting the input
 *         multiarray as a matrix with the appropriate extents.
 *  \return The result of the mm function call with the same number of columns as the reinterpreted input.
 *
 *  The first extent **must** be equal to `ext_1` of the `a` matrix.
 *  See comments in \ref constructor_mm_NN1C_Matrix_T for the preset matrix-matrix multiplication parameters.
 */
struct Multiarray_T* constructor_mm_NN1C_Multiarray_T
	(const struct const_Matrix_R*const a,    ///< Defined for \ref mm_d.
	 const struct const_Multiarray_T*const b ///< Input `b` in multiarray format.
	);

/** \brief `const` version of \ref constructor_mm_NN1C_Multiarray_T.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_mm_NN1C_const_Multiarray_T
	(const struct const_Matrix_R*const a,    ///< Defined for \ref  constructor_mm_NN1C_Multiarray_T.
	 const struct const_Multiarray_T*const b ///< Defined for \ref  constructor_mm_NN1C_Multiarray_T.
	);

/** \brief Constructor for a \ref const_Multiarray_T\* using a matrix-matrix multiplication, interpreting the input
 *         multiarray as a matrix with the appropriate extents.
 *  \return The result of the mm function call with the same number of columns as the reinterpreted input.
 *
 *  The first extent **must** be equal to `ext_1` of the `a` matrix.
 *  See comments in \ref constructor_mm_NN1C_Matrix_T for the preset matrix-matrix multiplication parameters, excluding
 *  the layout.
 */
/// \todo Remove if unused.
/*const struct const_Multiarray_T* constructor_mm_NN1_const_Multiarray_T
	(const struct const_Matrix_T*const a,     ///< Defined for \ref mm_d.
	 const struct const_Multiarray_T*const b, ///< Input `b` in multiarray format.
	 const char layout_c                      ///< The desired layout for the output.
	);*/

/** \brief Constructor for a \ref const_Multiarray_T\* by applying sub-operator matrices along each direction.
 *  \return Standard.
 *
 *  See comments in \ref constructor_mm_NN1C_Matrix_T for the preset matrix-matrix multiplication parameters.
 *
 *  The value of `b->extents[0]` must be equal to the product of `ext_1` of the matrices in `a_tp`. If the input `b` has
 *  an order greater than 1, the output of this function can be interpreted as the sub-operators being applied to each
 *  column individually.
 */
const struct const_Multiarray_T* constructor_mm_tp_NN1C_const_Multiarray_T
	(const struct const_Multiarray_Matrix_T* a_tp, ///< The tensor-product sub-operators.
	 const struct const_Multiarray_T* b            ///< The input multiarray.
	);

// Destructors ****************************************************************************************************** //

/// \brief Destructs a \ref Multiarray_T\*.
void destructor_Multiarray_T
	(struct Multiarray_T* a ///< Standard.
	);

/// \brief Destructs a \ref const_Multiarray_T\*.
void destructor_const_Multiarray_T
	(const struct const_Multiarray_T* a ///< Standard.
	);

/// \brief Destructs a \ref Multiarray_Vector_T\*.
void destructor_Multiarray_Vector_T
	(struct Multiarray_Vector_T* a ///< Standard.
	);

/// \brief `const` version of \ref destructor_Multiarray_Vector_T.
void destructor_const_Multiarray_Vector_T
	(const struct const_Multiarray_Vector_T* a ///< Defined for \ref destructor_Multiarray_Vector_T.
	);

/// \brief Destructs a \ref Multiarray_Matrix_T\*.
void destructor_Multiarray_Matrix_T
	(struct Multiarray_Matrix_T* a ///< Standard.
	);

/// \brief `const` version of \ref destructor_Multiarray_Matrix_T.
void destructor_const_Multiarray_Matrix_T
	(const struct const_Multiarray_Matrix_T* a ///< Defined for \ref destructor_Multiarray_Matrix_T.
	);

/// \brief Destructor for a \ref Multiarray_Matrix_T\*[2].
void destructor_const_Multiarray2_Matrix_T
	(const struct const_Multiarray_Matrix_T* a[2] ///< Standard.
	);
