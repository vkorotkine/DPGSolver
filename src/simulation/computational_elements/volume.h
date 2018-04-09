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

#ifndef DPG__volume_h__INCLUDED
#define DPG__volume_h__INCLUDED
/** \file
 *  \brief Provides the interface for the base \ref Volume container and associated functions.
 *
 *  A \ref Volume is a `d` dimensional finite element.
 */

#include <stdbool.h>

#include "definitions_elements.h"
#include "intrusive.h"

struct Mesh;
struct Simulation;
struct const_Vector_i;
struct const_Multiarray_Vector_i;

/// \brief Container for data relating to the base Volumes.
struct Volume {
	struct Intrusive_Link lnk; ///< The \ref Intrusive_Link.

	const int index; ///< The index of the volume.

	const struct const_Element*const element; ///< Pointer to the associated \ref const_Element.

	/// Boundary conditions of each of the volume faces. Equal to \ref BC_INVALID for non-boundary faces.
	const struct const_Vector_i* bc_faces;

	/// Boundary conditions of each of the volume edges. Equal to \ref BC_INVALID for non-boundary edges.
	const struct const_Vector_i* bc_edges;

	const bool boundary, ///< Flag for whether the volume is on a domain boundary.
	           curved;   ///< Flag for whether the volume is curved.

	const struct const_Multiarray_d*const xyz_ve; ///< The xyz coordinates of the volume vertices.

	double h; ///< Measure of volume size (maximum edge length).

	const struct Face*const faces[NFMAX][NSUBFMAX]; ///< Array of pointers to the neighbouring \ref Face containers.
};

// Constructor/Destructor functions ********************************************************************************* //

/** \brief Constructs the base \ref Volume \ref Intrusive_List.
 *  \return Standard. */
struct Intrusive_List* constructor_Volumes
	(struct Simulation*const sim, ///< The \ref Simulation.
	 const struct Mesh*const mesh ///< The \ref Mesh.
	);

/// \brief Destructs the base \ref Volume \ref Intrusive_List.
void destructor_Volumes
	(struct Intrusive_List* volumes ///< Standard.
	);

/// \brief Destructor for a \ref Volume.
void destructor_Volume
	(struct Volume* vol ///< Pointer to the volume.
	);

// Helper functions ************************************************************************************************* //

/** \brief Check if a sufficient number of vertices (2) satisfy the condition on a domain boundary.
 *
 *  When two vertices satisfy the condition, this implies that a volume edge satisfies the condition and the
 *  qualification is passed to the associated volume.
 *
 *  This function can be used to check:
 *	- faces of volumes: `ve_inds` = volume vertex indices, `f_ve` = \ref Volume::element;
 *	- edges of faces:   `ve_inds` = face vertex indices,   `f_ve` = \ref Face::element.
 *
 *  Currently used to check:
 *	- curved (curved_only = `true`);
 *	- boundary (curved only = `false`).
 *
 *  \return See brief. */
bool check_ve_condition
	(const struct const_Multiarray_Vector_i*const f_ve,  ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,          ///< The vertex indices for the volume.
	 const struct const_Vector_i*const ve_condition,     ///< The vertex condition to check for.
	 const struct const_Multiarray_Vector_i*const ve_bc, ///< Defined in \ref Mesh_Vertices.
	 const bool curved_only                              /**< Flag indicating whether only curved boundary
	                                                      *   conditions should be considered. */
	);

/// \brief Update the \ref Element pointed to by each of the \ref Volume\*s.
void update_volumes_element
	(struct Intrusive_List* volumes,             ///< The list of \ref Volume\*s to be updated.
	 const struct const_Intrusive_List* elements ///< The list of elements from which to replace the current ones.
	);

/** \brief Copy constructor for a \ref Volume.
 *  \return Standard.
 *  The pointers for \ref Volume::faces **are all set to NULL**.
 *  The pointer for \ref Volume::element is set according to the value of the `independent_elements` flag.
 */
struct Volume* constructor_copy_Volume
	(const struct Volume*const vol_i,   ///< The input \ref Volume.
	 const struct Simulation*const sim, ///< \ref Simulation.
	 const bool independent_elements    ///< Flag for whether
	);

/** \brief Compute and return the \ref Volume::h based on the element type and the vertex coordinates.
 *  \return See brief. */
double compute_h_volume
	(const struct Volume*const vol ///< \ref Volume.
	);

#endif // DPG__volume_h__INCLUDED
