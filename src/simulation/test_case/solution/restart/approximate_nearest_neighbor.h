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

#ifndef DPG__approximate_nearest_neighbor_h__INCLUDED
#define DPG__approximate_nearest_neighbor_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions relating to approximate nearest neighbor computation.
 *
 *  \note The only reason for the reimplementation of this functionality was a lack of c libraries with the desired
 *        features; it was determined to be simpler to do the implementation than to maintain a c wrapper for an
 *        existing c++ library. If this project is ever converted to c++, please consider directly using an available
 *        library which would almost certainly be more efficient and more extensively tested.
 */

#include <stddef.h>
#include <float.h>
#include <limits.h>

#include "definitions_core.h"

struct const_Vector_i;
struct const_Matrix_d;
struct Matrix_d;

/// Use `int` type if true and `long long int` otherwise.
#define USE_SINGLE false

#define Real double
#define REAL_MAX DBL_MAX

/** \{ \name Templating-related definitions.
 *
 *  The division by 4 in the maximum indices is required because:
 *  - The shift results in all numbers potentially being increased (maximally) by a factor of 2;
 *  - Something related to the subtraction requires an additional factor of 2?
 */
#if USE_SINGLE == true
	#define Index     int
	#define INDEX_MIN 0
	#define INDEX_MAX (INT_MAX/4)
	#define R_INDEX_MIN ((Real)0)
	#define R_INDEX_MAX ((Real)(INT_MAX/4))

	/** The accuracy tolerance for the approximate nearest neighbor search. Chosen based on the limitations of the
	 *  conversion to `int`. */
	#define EPS_ANN 1e-9
#else
	#define Index     long long
	#define INDEX_MIN 0
	#define INDEX_MAX (LLONG_MAX/4)
	#define R_INDEX_MIN ((Real)0)
	#define R_INDEX_MAX ((Real)(LLONG_MAX/4))

	#define EPS_ANN 1e-16
#endif
/// \}

/// \brief Container for 'A'pproximate 'N'earest 'N'eighbor input data.
struct Input_ANN {
	const struct const_Matrix_d* nodes_b; ///< Background nodes.
	const struct const_Matrix_d* nodes_s; ///< Search nodes.
};

/// \brief Container for a node represented in binary.
struct Node_ANN {
	Index xyz[DIM]; ///< The coordinates with type `Index`.
	int index;      ///< The index.
};

/// \brief Container for 'S'hift-'S'huffle-'S'ort 'A'pproximate 'N'earest 'N'eighbor information.
struct SSS_ANN {
	/** Matrix holding the values of the minimum and maximum node coordinates in each of the xyz coordinate
	 *  directions; the last two rows are used to hold the average ((min+max)/2) and the difference (max-min). */
	const struct const_Matrix_d* xyz_min_max;

	ptrdiff_t ext_b; ///< The number of entries in \ref SSS_ANN::b.
	ptrdiff_t ext_s; ///< The number of entries in \ref SSS_ANN::s.

	const struct Node_ANN* b; ///< The background nodes.
	const struct Node_ANN* s; ///< The search nodes.
};

/// \brief Container for output data from the constructor for a sorted list of nodes.
struct Nodes_Sorted_ANN {
	const struct const_Matrix_d* nodes;   ///< The sorted nodes.
	const struct const_Vector_i* indices; ///< The indices of the sorted nodes in relation to the input nodes.
};

/** \brief Constructor for the \ref const_Vector_T\* (`int`) of 'a'pproximate 'n'earest 'n'eighbor indices for the input
 *         background and search node lists.
 *  \return See brief.
 *
 *  As the current implementation relies on the conversion of input scalar coordinates to integer values, it is expected
 *  that this function should begin to give incorrect results when the length scales of the input nodes vary by more
 *  than ~1e9.
 */
const struct const_Vector_i* constructor_ann_indices
	(const struct Input_ANN*const ann_i ///< Standard.
	);

/** \brief Version of \ref constructor_ann_indices but using a \ref SSS_ANN input.
 *  \return See brief. */
const struct const_Vector_i* constructor_ann_indices_from_sss
	(const struct SSS_ANN*const sss ///< Standard.
	);

/// \brief Destructor for a \ref Input_ANN container.
void destructor_Input_ANN
	(struct Input_ANN*const ann_info ///< Standard.
	);

/** \brief Constructor for \ref SSS_ANN::xyz_min_max.
 *  \return See brief. */
void constructor_SSS_ANN_xyz
	(const struct Input_ANN*const ann_i, ///< Standard.
	 struct SSS_ANN*const sss            ///< Standard.
	);

/// \brief Destructor for \ref SSS_ANN::xyz_min_max.
void destructor_SSS_ANN_xyz
	(const struct SSS_ANN*const sss ///< Standard.
	);

/** \brief Constructor for the background portion of the input \ref SSS_ANN container.
 *  \return See brief. */
void constructor_SSS_ANN_b
	(const struct Input_ANN*const ann_i, ///< Standard.
	 struct SSS_ANN*const sss            ///< Standard.
	);

/// \brief Destructor for the background portion of the input \ref SSS_ANN container.
void destructor_SSS_ANN_b
	(const struct SSS_ANN*const sss ///< Standard.
	);

/** \brief Constructor for the search portion of the input \ref SSS_ANN container.
 *  \return See brief. */
void constructor_SSS_ANN_s
	(const struct Input_ANN*const ann_i, ///< Standard.
	 struct SSS_ANN*const sss            ///< Standard.
	);

/// \brief Destructor for the search portion of the input \ref SSS_ANN container.
void destructor_SSS_ANN_s
	(const struct SSS_ANN*const sss ///< Standard.
	);

/// \brief Sort the input list of nodes using the 'A'pproximate 'N'earest 'N'eighbor shuffle comparison.
void sort_nodes_ANN
	(const ptrdiff_t n_n,        ///< The length of the array of nodes.
	 struct Node_ANN*const nodes ///< The array of nodes.
	);

/** \brief Constructor for a \ref Nodes_Sorted_ANN container.
 *  \return See brief. */
const struct Nodes_Sorted_ANN* constructor_Nodes_Sorted_ANN
	(const struct const_Matrix_d*const nodes_i ///< Input unsorted nodes.
	);

/** \brief Version of \ref constructor_Nodes_Sorted_ANN transposing the input nodes to be row-major and then back if the
 *         input is not row-major.
 *  \return See brief. */
const struct Nodes_Sorted_ANN* constructor_Nodes_Sorted_ANN_with_trans
	(struct Matrix_d*const nodes_i ///< See brief.
	);

/// \brief Destructor for a \ref Nodes_Sorted_ANN container.
void destructor_Nodes_Sorted_ANN
	(const struct Nodes_Sorted_ANN*const nsa ///< Standard.
	);

#endif // DPG__approximate_nearest_neighbor_h__INCLUDED
