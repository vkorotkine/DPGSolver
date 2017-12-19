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

#ifndef DPG__face_h__INCLUDED
#define DPG__face_h__INCLUDED
/** \file
 *  \brief Provides the interface for the base \ref Face container and associated functions.
 *
 *  A \ref Face is a `d-1` dimensional finite element which is found on the face of a \ref Volume.
 */

#include <stdbool.h>

#include "intrusive.h"

struct Mesh;
struct Simulation;

/// \brief Container for data relating to the base Faces.
struct Face {
	struct Intrusive_Link lnk; ///< \ref Intrusive_Link.

	const int index; ///< The index of the face.

	const bool boundary, ///< Flag for whether the face is on a domain boundary.
	           curved;   ///< Flag for whether the face is curved.

	const int bc;        ///< The boundary condition associated with the face (if relevant).

	const struct const_Element*const element; ///< Pointer to the associated \ref const_Element.

	/** \brief Container for information relating to the neighbouring \ref Volume on either side of the \ref Face.
	 *
	 *  The information for the first index `neigh_info[0]` relates to the Volume whose outward normal vector on the
	 *  current face coincides with that stored as part of the Face; this is generally referred to as the left
	 *  volume.
	 */
	struct Neigh_Info {
		/// Local face index in relation to the neighbouring volume.
		int ind_lf;

		/** Local face h-refinement index.
		 *  Range: [0:NFREFMAX).
		 *  In the case of a conforming mesh (h-refinement disabled), ind_href = 0. */
		int ind_href;

		/** Local sub-face h-refinement index.
		 *  Range: [0:NFSUBFMAX).
		 *
		 *  Example:
		 *
		 *  For a Face of type QUAD, there are nine possible values for ind_sref, related to the allowed
		 *  h-refinements:
		 *  	None       (1 -> 1): ind_href = [0]
		 *  	Isotropic  (1 -> 4): ind_href = [1,4]
		 *  	Horizontal (1 -> 2): ind_href = [5,6]
		 *  	Vertical   (1 -> 2): ind_href = [7,8]
		 *
		 *  	ind_href = 0 : Indsfh = 0 (This is a conforming Face)
		 *  	           1 : Indsfh = 0 (This is the first  sub-face of the h-refined macro Face)
		 *  	           2 : Indsfh = 1 (This is the second sub-face of the h-refined macro Face)
		 *  	           3 : Indsfh = 2 (This is the third  sub-face of the h-refined macro Face)
		 *  	           4 : Indsfh = 3 (This is the fourth sub-face of the h-refined macro Face)
		 *  	           5 : Indsfh = 0
		 *  	           6 : Indsfh = 1
		 *  	           7 : Indsfh = 0
		 *  	           8 : Indsfh = 1
		 */
		int ind_sref;

		/** The local face ordering index.
		 *  Specifies the index for the ordering between Faces as seen from adjacent Volumes. When
		 *  interpolating values from Volumes to Faces, it is implicitly assumed that the Face is in the
		 *  reference configuration of the current Volume. When these interpolated values must be used in
		 *  relation to the neighbouring Face, their orientation is likely to be incorrect and the reordering
		 *  necessary to be seen correctly is required. This variable provides the index to the ordering array
		 *  which transfers the current ordering to that required as if the Face were seen from the opposite
		 *  Volume.
		 */
		int ind_ord;

		/// Pointer to the neighbouring \ref Volume. The second pointer is `NULL` for a boundary face.
		struct Volume *volume;
	} neigh_info[2]; ///< \ref Neigh_Info.
};

// Constructor/Destructor functions ********************************************************************************* //

/** \brief Constructs the base \ref Face \ref Intrusive_List.
 *  \return Standard. */
struct Intrusive_List* constructor_Faces
	(struct Simulation*const sim, ///< The \ref Simulation.
	 const struct Mesh*const mesh ///< The \ref Mesh.
	);

/// \brief Destructs the base \ref Face \ref Intrusive_List.
void destructor_Faces
	(struct Intrusive_List* faces ///< Standard.
	);

/// \brief Cast from \ref Face\* to `const` \ref Face `*const`.
void const_cast_Face
	(const struct Face*const* dest, ///< Destination.
	 const struct Face*const src    ///< Source.
	);

// Helper functions ************************************************************************************************* //

/** \brief Call \ref get_face_element_index_by_ind_lf using parameters from the input \ref Face.
 *  \return See brief. */
int get_face_element_index
	(const struct Face*const face ///< The input face.
	);

/** \brief Compute the side index of the current face. Options: 0 (input left volume), 1 (input right volume).
 *  \return See brief. */
int compute_side_index_face
	(const struct Face* face, ///< \ref Face.
	 const struct Volume* vol ///< The pointer to the neighbouring \ref Volume.
	);

/** \brief Copy constructor for a \ref Face.
 *  \return Standard.
 *  The pointers for \ref Face::Neigh_Info::volume **are all set to NULL**.
 *  The pointer for \ref Face::element is set according to the value of the `independent_elements` flag.
 */
struct Face* constructor_copy_Face
	(const struct Face*const face_i,    ///< The input \ref Face.
	 const struct Simulation*const sim, ///< \ref Simulation.
	 const bool independent_elements    ///< Flag for whether
	);

#endif // DPG__face_h__INCLUDED
