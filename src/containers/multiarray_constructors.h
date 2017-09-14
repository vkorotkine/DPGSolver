// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__multiarray_constructors_h__INCLUDED
#define DPG__multiarray_constructors_h__INCLUDED
/** \file
 *  \brief Provides Multiarray_\* constructors and destructors.
 *
 *  \section s1_Multi_con Naming Convention
 *
 *  Function names are chosen according to the following template: `constructor_{0}_(1)_{2}_(3)_(4)` where elements
 *  in curly braces `{}` are required and those in round brackets `()` are optional.
 *  	- {0} : Type of constructor
 *  		- default: reserve memory only for the container itself (not for the data).
 *  		- empty:   reserve memory storage for the data but do not set it.
 *  		- copy:    reserve memory storage for the data and set it to a copy of the input data.
 *  		- move:    move data to the container being constructed.
 *  	- (1) : Optional `const` specifier
 *  	- {2} : Type of container to be returned
 *  		- Multiarray_\*: d, Vector_i
 *  	- (3) : Level of dereferencing if not equal to 1.
 *  	- (4) : Type of input from which the container is constructed
 *
 *  \section s2_Multi_con Const Versions of Multiarray Containers
 *
 *  In order to avoid unintentionally overwriting data which should be constant, `const` versions of the containers
 *  are also provided where relevant. A complication which arises as a result of declaring objects `const` is that
 *  it is not possible to define them! To overcome this difficulty, casts are used to set the data.
 *
 *  \warning This procedure exhibits undefined behaviour (relating to changing a `const`-qualified type) unless the
 *  memory for the struct is **dynamically allocated**. See [this SO answer][SO_dyn_const_struct] for a detailed
 *  explanation of the problem and possible approaches used to overcome it.
 *
 *  \warning It is **required** that `const` and non-`const` versions of the containers have identical memory layout
 *           such that casts can be used to convert between them.
 *
 *  Noting that containers should be dynamically allocated, we have the implication that containers with no
 *  dereferencing should never be used. Taken together with the upper limitation on the level of dereferencing, it
 *  is thus required that **containers have exactly one level of dereferencing**. The only place in which an
 *  additional level of dereferencing is permitted is as a list of pointers to containers in a multiarray.
 *
 *  \subsection s21_Multi_con Constructors for `const` Containers
 *
 *  The `const_constructor_move_*` functions are used to *define* the `const` equivalent of the container. The
 *  functions move **all** container members through casts.
 *
 *  \note This means that a single destructor call should be made for both the `src` and `dest` variables.
 *
 *  \section s3_Multi_con Compound Literals
 *
 *  In the interest of greater generic programming, compound literals should be passed as multiarray constructor
 *  arguments when `order > 2`. Example implementations can be found in questions and answers of
 *  [this][SO_compound1] and [this][SO_compound2] SO discussion.
 *
 *  <!-- References: -->
 *  [SO_dyn_const_struct]: https://stackoverflow.com/q/2219001/5983549
 *  [SO_compound1]: https://stackoverflow.com/q/17374792/5983549
 *  [SO_compound2]: https://stackoverflow.com/a/11268460/5983549
 */

#include <stddef.h>
#include <stdbool.h>

struct Vector_i;
struct Matrix_d;
struct Multiarray_Matrix_d;
struct const_Matrix_d;
struct const_Multiarray_d;
struct const_Multiarray_Vector_i;
struct const_Multiarray_Matrix_d;

// Default constructors ********************************************************************************************* //

/** \brief Constructor for a default \ref Multiarray_d\*.
 *  \return Standard. */
struct Multiarray_d* constructor_default_Multiarray_d ();

/** \brief Constructor for a default \ref Multiarray_Matrix_d\*.
 *  \return Standard. */
struct Multiarray_Matrix_d* constructor_default_Multiarray_Matrix_d ();

/** \brief `const` version of \ref constructor_default_Multiarray_Matrix_d.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_default_const_Multiarray_Matrix_d ();

// Empty constructors *********************************************************************************************** //

/** \brief Constructor for an empty \ref Multiarray_d\*.
 *  \return Standard.
 */
struct Multiarray_d* constructor_empty_Multiarray_d
	(const int order,                ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

/** \brief Constructor for an empty \ref Multiarray_Vector_i\*.
 *  \return Standard. */
struct Multiarray_Vector_i* constructor_empty_Multiarray_Vector_i
	(const bool alloc_V,             ///< Flag for whether memory should be reserved for the individual Vectors.
	 const int order,                ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

/** \brief Constructor for an empty \ref Multiarray_Matrix_d\*.
 *  \return Standard. */
struct Multiarray_Matrix_d* constructor_empty_Multiarray_Matrix_d
	(const bool alloc_M,             ///< Flag for whether memory should be reserved for the individual Matrices.
	 const int order,                ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

// Copy constructors ************************************************************************************************ //

/** \brief Constructs a \ref Multiarray_Vector_i\* and sets the values of its \ref Vector_i\* components from the
 *         intput `int*` data.
 *  \return Standard. */
struct Multiarray_Vector_i* constructor_copy_Multiarray_Vector_i_i
	(const int* data_V,              ///< The `int` data.
	 const int*const ext_V,          /**< The ext_0 values of each \ref Vector_i in the multiarray. Note: using
	                                      `int` and not `ptrdiff_t` as this is sometimes initialized from
	                                      \ref Vector_i::data. */
	 const int order,                ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

/** \brief Copy constructor for a `const` \ref const_Multiarray_d\* from a `const` \ref const_Multiarray_d\*.
 *  \return Standard. */
void const_constructor_copy_Multiarray_d
	(const struct const_Multiarray_d*const* dest, ///< Destination.
	 const struct const_Multiarray_d*const src    ///< Source.
	);

// Move constructors ************************************************************************************************ //

/** \brief Move constructor for a \ref Multiarray_d\* from a `double*`.
 *  \return Standard. */
struct Multiarray_d* constructor_move_Multiarray_d_d
	(const char layout,               ///< Standard.
	 const int order,                 ///< Standard.
	 const ptrdiff_t*const extents_i, ///< The input extents.
	 const bool owns_data,            ///< Standard.
	 double*const data                ///< Standard.
	);

/** \brief Move constructor for a \ref Multiarray_Vector_i\* with the input extents having been previously
 *         dynamically allocated.
 *  \return See brief. */
struct Multiarray_Vector_i* constructor_move_Multiarray_Vector_i_dyn_extents
	(const int order,            ///< Standard.
	 ptrdiff_t*const extents,    ///< Standard.
	 const bool owns_data,       ///< Standard.
	 struct Vector_i**const data ///< Standard.
	);

/** \brief Move constructor for a \ref Multiarray_Matrix_d\* with the input extents having been previously
 *         dynamically allocated.
 *  \return See brief. */
struct Multiarray_Matrix_d* constructor_move_Multiarray_Matrix_d_dyn_extents
	(const int order,            ///< Standard.
	 ptrdiff_t*const extents,    ///< Standard.
	 const bool owns_data,       ///< Standard.
	 struct Matrix_d**const data ///< Standard.
	);

/** \brief Move constructor for a \ref Multiarray_d\* with the input extents having been previously dynamically
 *         allocated.
 *  \return See brief. */
struct Multiarray_d* constructor_move_Multiarray_d_dyn_extents
	(const char layout,       ///< Standard.
	 const int order,         ///< Standard.
	 ptrdiff_t*const extents, ///< Standard.
	 const bool owns_data,    ///< Standard.
	 double*const data        ///< Standard.
	);

/** \brief Move constructor for a \ref Multiarray_d\* from a \ref Matrix_d\*.
 *  \return Standard. */
struct Multiarray_d* constructor_move_Multiarray_d_Matrix_d
	(struct Matrix_d* src ///< The source matrix.
	);

/// \brief Move Constructor for a `const` \ref const_Multiarray_d `*const`.
void const_constructor_move_Multiarray_d
	(const struct const_Multiarray_d*const* dest, ///< Destination.
	 struct Multiarray_d* src                     ///< Source.
	);

/// \brief Move constructor for a `const` \ref const_Multiarray_Vector_i `*const`.
void const_constructor_move_Multiarray_Vector_i
	(const struct const_Multiarray_Vector_i*const* dest, ///< Destination.
	 struct Multiarray_Vector_i* src                     ///< Source.
	);

// Special constructors ********************************************************************************************* //

/// \brief Set a \ref Multiarray_Matrix_d\* from a sub range of a \ref Multiarray_Matrix_d\*.
void set_Multiarray_Matrix_from_Multiarray_Matrix_d
	(struct Multiarray_Matrix_d* dest, ///< The destination.
	 struct Multiarray_Matrix_d* src,  ///< The source.
	 const int order_o,                ///< The order of the output (destination).
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Multiarray_Matrix_from_Multiarray_Matrix_d.
void set_const_Multiarray_Matrix_from_Multiarray_Matrix_d
	(const struct const_Multiarray_Matrix_d* dest, ///< Defined for mutable version.
	 const struct const_Multiarray_Matrix_d* src,  ///< Defined for mutable version.
	 const int order_o,                            ///< Defined for mutable version.
	 const ptrdiff_t*const sub_indices             ///< Defined for mutable version.
	);

/// \brief Set a \ref Matrix_d\* from an entry of a \ref Multiarray_Matrix_d\*.
void set_Matrix_from_Multiarray_Matrix_d
	(struct Matrix_d* dest,            ///< The destination.
	 struct Multiarray_Matrix_d* src,  ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Matrix_from_Multiarray_Matrix_d.
void set_const_Matrix_from_Multiarray_Matrix_d
	(const struct const_Matrix_d* dest,           ///< Defined for \ref set_Matrix_from_Multiarray_Matrix_d.
	 const struct const_Multiarray_Matrix_d* src, ///< Defined for \ref set_Matrix_from_Multiarray_Matrix_d.
	 const ptrdiff_t*const sub_indices            ///< Defined for \ref set_Matrix_from_Multiarray_Matrix_d.
	);

// Destructors ****************************************************************************************************** //

/// \brief Destructs a \ref Multiarray_d\*.
void destructor_Multiarray_d
	(struct Multiarray_d* a ///< Standard.
	);

/// \brief Destructs a \ref const_Multiarray_d\*.
void destructor_const_Multiarray_d
	(const struct const_Multiarray_d* a ///< Standard.
	);

/// \brief Destructs a \ref Multiarray_Vector_i\*.
void destructor_Multiarray_Vector_i
	(struct Multiarray_Vector_i* a ///< Standard.
	);

/// \brief `const` version of \ref destructor_Multiarray_Vector_i.
void destructor_const_Multiarray_Vector_i
	(const struct const_Multiarray_Vector_i* a ///< Defined for \ref destructor_Multiarray_Vector_i.
	);

/// \brief Destructs a \ref Multiarray_Matrix_d\*.
void destructor_Multiarray_Matrix_d
	(struct Multiarray_Matrix_d* a ///< Standard.
	);

/// \brief `const` version of \ref destructor_Multiarray_Matrix_d.
void destructor_const_Multiarray_Matrix_d
	(const struct const_Multiarray_Matrix_d* a ///< Defined for \ref destructor_Multiarray_Matrix_d.
	);

#endif // DPG__multiarray_constructors_h__INCLUDED
