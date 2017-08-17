// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Matrix_h__INCLUDED
#define DPG__Matrix_h__INCLUDED
/**	\file
 *	\brief Provides Matrix_\* containers and related functions.
 *
 *	Potentially relevant comments may be found in containers.h.
 *
 *	Matrices are 2D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>

/// \brief Matrix (`double`).
struct Matrix_d {
	char layout;

	size_t extents[2];

	bool    owns_data;
	double* data;
};

/// \brief Matrix (`const double`).
struct const_Matrix_d {
	const char layout;

	const size_t extents[2];

	const bool         owns_data;
	const double*const data;
};

/// \brief Matrix (`unsigned int`).
struct Matrix_ui {
	char layout;

	size_t extents[2];

	bool    owns_data;
	unsigned int* data;
};

/// \brief Matrix (`const unsigned int`).
struct const_Matrix_ui {
	const char layout;

	const size_t extents[2];

	const bool    owns_data;
	const unsigned int* data;
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructs an empty \ref Matrix_d\*.
struct Matrix_d* constructor_empty_Matrix_d
	(const char layout,  ///< Standard.
	 const size_t ext_0, ///< Standard.
	 const size_t ext_1  ///< Standard.
	);

/// \brief Move Constructor for a \ref const_Matrix_d\*.
struct const_Matrix_d* constructor_move_const_Matrix_d_Matrix_d
	(struct Matrix_d*const src ///< Source Matrix.
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

/// \brief Constructs an empty \ref Matrix_ui\*.
struct Matrix_ui* constructor_empty_Matrix_ui
	(const char layout,  ///< Standard.
	 const size_t ext_0, ///< Standard.
	 const size_t ext_1  ///< Standard.
	);

/// \brief Move Constructor for a \ref const_Matrix_ui\*.
struct const_Matrix_ui* constructor_move_const_Matrix_ui_Matrix_ui
	(struct Matrix_ui*const src ///< Source Matrix.
	);

/// \brief Move Constructor for a `const` \ref const_Matrix_ui `*const`.
void const_constructor_move_Matrix_ui
	(const struct const_Matrix_ui*const* dest, ///< Destination.
	 struct Matrix_ui* src                     ///< Source.
	);

/// \brief Destructs a \ref Matrix_ui\*.
void destructor_Matrix_ui
	(struct Matrix_ui* a ///< Standard.
	);

// Helper functions ************************************************************************************************* //

/** \brief Get pointer to row of row-major \ref Matrix_d\*.
 *	\return Pointer to the first entry of the row.
 */
double* get_row_Matrix_d
	(const size_t row,        ///< Desired row.
	 const struct Matrix_d* a ///< Matrix.
	);

/** \brief Get pointer to row of row-major \ref const_Matrix_d\*.
 *	\return Pointer to the first entry of the row.
 */
const double* get_row_const_Matrix_d
	(const size_t row,        ///< Desired row.
	 const struct const_Matrix_d* a ///< Matrix.
	);

/** \brief Get pointer to row of row-major \ref Matrix_ui\*.
 *	\return Pointer to the first entry of the row.
 */
unsigned int* get_row_Matrix_ui
	(const size_t row,         ///< Desired row.
	 const struct Matrix_ui* a ///< Matrix.
	);

/** \brief Get value of the (row,col) entry of a \ref Matrix_ui\*.
 *	\return See brief. */
unsigned int get_val_Matrix_ui
	(const size_t row,              ///< The row.
	 const size_t col,              ///< The column.
	 const struct Matrix_ui*const a ///< Standard.
	);

/** \brief Get value of the (row,col) entry of a \ref const_Matrix_ui\*.
 *	\return See brief. */
unsigned int get_val_const_Matrix_ui
	(const size_t row,              ///< The row.
	 const size_t col,              ///< The column.
	 const struct const_Matrix_ui*const a ///< Standard.
	);

// Printing functions *********************************************************************************************** //

/// \brief Print a \ref Matrix_d\* to the terminal.
void print_Matrix_d
	(const struct Matrix_d*const a ///< Standard.
	);

/// \brief Print a \ref const_Matrix_d\* to the terminal.
void print_const_Matrix_d
	(const struct const_Matrix_d*const a ///< Standard.
	);

/// \brief Print a \ref Matrix_ui\* to the terminal.
void print_Matrix_ui
	(const struct Matrix_ui*const a ///< Standard.
	);

/// \brief Print a \ref const_Matrix_ui\* to the terminal.
void print_const_Matrix_ui
	(const struct const_Matrix_ui*const a ///< Standard.
	);


#endif // DPG__Matrix_h__INCLUDED
