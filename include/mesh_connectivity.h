// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__mesh_connectivity_h__INCLUDED
#define DPG__mesh_connectivity_h__INCLUDED
/**	\file
 *	Provides the interface to mesh connectivity containers and functions.
 */

#include "mesh_readers.h"
#include "Intrusive.h"

#include "Multiarray.h"
#include "Vector.h"

/// \brief The number of 'M'aster and 'S'lave entities (always two but defined to avoid use of the magic number).
#define N_MS 2

/// \brief Container for local connectivity related information.
struct Conn_info {
	// Available from mesh_data:
	const unsigned int d;                 ///< The dimension.
	struct Vector_ui* elem_per_dim;       ///< The number of elements of each dimension.
	struct const_Vector_ui* volume_types; ///< Pointer to the first volume entry in \ref Mesh_Data::elem_types.
	struct Vector_ui* v_n_lf;             ///< The number of local faces for each volume.

	// Computed here:
	struct Multiarray_Vector_ui* f_ve;     ///< Global face to vertex correspondence.
	struct Vector_ui*            ind_f_ve; ///< Indices of \ref f_ve after sorting.
};

/// \brief Holds data relating to the mesh connectivity.
struct Mesh_Connectivity {
	const struct const_Multiarray_Vector_ui*const v_to_v;  ///< Volume to volume connectivity.

	/** Volume to local face connectvity.
	 *	Redundant self-reference entries are replaced with the number corresponding to the appropriate boundary
	 *	condition. */
	const struct const_Multiarray_Vector_ui*const v_to_lf;
};

/// \brief Set up the mesh connectivity.
struct Mesh_Connectivity* mesh_connect
	(const struct Mesh_Data*const mesh_data,     ///< \ref Mesh_Data.
	 const struct const_Intrusive_List* elements ///< The base \ref Element list.
	);

/// \brief Destructor for \ref Mesh_Connectivity.
void destructor_Mesh_Connectivity
	(struct Mesh_Connectivity* mesh_conn ///< Standard.
	);

/** \brief See return.
 *	\return The index of the first volume. */
size_t get_first_volume_index
	(const struct Vector_ui*const elem_per_dim, ///< Defined in \ref Conn_info.
	 const unsigned int d                       ///< Defined in \ref Conn_info.
	);

/// \brief Set the face node numbers and sort them.
void set_f_node_nums
	(struct Vector_ui**const f_node_nums,         ///< The face node numbers.
	 const struct const_Vector_ui*const node_nums ///< The entry of \ref Mesh_Data::node_nums for the current face.
	);

#endif // DPG__mesh_connectivity_h__INCLUDED
