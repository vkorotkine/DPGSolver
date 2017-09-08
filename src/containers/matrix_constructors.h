// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Matrix_h__INCLUDED
#define DPG__Matrix_h__INCLUDED
/**	\file
 *	\brief Provides Matrix_\* constructors and destructors.
 *
 *	Matrices are 2D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>

struct const_Vector_i;
struct const_Matrix_i;

// Default constructors ********************************************************************************************* //

/**	\brief Constructor for a default \ref Matrix_d\*.
 *	\return Standard. */
struct Matrix_d* constructor_default_Matrix_d ();

// Empty constructors *********************************************************************************************** //

/**	\brief Constructs an empty \ref Matrix_d\*.
 *	\return Standard. */
struct Matrix_d* constructor_empty_Matrix_d
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1  ///< Standard.
	);

/**	\brief Constructs an empty \ref Matrix_i\*.
 *	\return Standard. */
struct Matrix_i* constructor_empty_Matrix_i
	(const char layout,  ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1  ///< Standard.
	);

// Copy constructors ************************************************************************************************ //

/**	\brief Copy constructor for a \ref Matrix_d\* from a \ref Matrix_d\*.
 *	\return Standard. */
struct Matrix_d* constructor_copy_Matrix_d
	(struct Matrix_d* src /// The source matrix.
	);

/**	\brief Copy constructor for a \ref Matrix_i\* from a `const int*`.
 *	\return Standard. */
struct Matrix_i* constructor_copy_Matrix_i_i
	(const char layout,       ///< Standard.
	 const ptrdiff_t ext_0,   ///< Standard.
	 const ptrdiff_t ext_1,   ///< Standard.
	 const int*const data_src ///< The source data.
	);

/**	\brief Copy constructor for a \ref Matrix_d\* from a `const double*`.
 *	\return Standard. */
struct Matrix_d* constructor_copy_Matrix_d_d
	(const char layout,          ///< Standard.
	 const ptrdiff_t ext_0,      ///< Standard.
	 const ptrdiff_t ext_1,      ///< Standard.
	 const double*const data_src ///< The source data.
	);

/**	\brief Copy constructor for a \ref const_Matrix_d\* from a partial number of rows/columns of another.
 *	\return Standard. */
const struct const_Matrix_d* constructor_copy_extract_const_Matrix_d
	(const struct const_Matrix_d*const src,    ///< The source Matrix.
	 const struct const_Vector_i*const indices ///< The indices of the rows/columns to copy.
	);

/**	\brief Copy constructor for a `const` \ref const_Matrix_d\* from a `const` \ref const_Matrix_d\*.
 *	\return Standard. */
void const_constructor_copy_Matrix_d
	(const struct const_Matrix_d*const* dest, ///< Destination.
	 const struct const_Matrix_d*const src    ///< Source.
	);

// Move constructors ************************************************************************************************ //

/**	\brief Move constructor for a \ref Matrix_d\* from a `double*`.
 *	\return Standard. */
struct Matrix_d* constructor_move_Matrix_d_d
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 double*const data      ///< Standard.
	);

/**	\brief Move constructor for a \ref Matrix_i\* from a `int*`.
 *	\return Standard. */
struct Matrix_i* constructor_move_Matrix_i_i
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 int*const data         ///< Standard.
	);

///	\brief Move Constructor for a `const` \ref const_Matrix_d `*const`.
void const_constructor_move_Matrix_d
	(const struct const_Matrix_d*const* dest, ///< Destination.
	 struct Matrix_d* src                     ///< Source.
	);

/**	\brief Move Constructor for a `const` \ref const_Matrix_i `*const`.
 *	\return Standard. */
void const_constructor_move_Matrix_i
	(const struct const_Matrix_i*const* dest, ///< Destination.
	 struct Matrix_i* src                     ///< Source.
	);

// Destructors ****************************************************************************************************** //

///	\brief Destructs a \ref Matrix_d\*.
void destructor_Matrix_d
	(struct Matrix_d* a ///< Standard.
	);

///	\brief Destructs a \ref Matrix_i\*.
void destructor_Matrix_i
	(struct Matrix_i* a ///< Standard.
	);

#endif // DPG__Matrix_h__INCLUDED
