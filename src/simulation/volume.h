// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__volume_h__INCLUDED
#define DPG__volume_h__INCLUDED
/**	\file
 *	\brief Provides the interface for the base \ref Volume container and associated functions.
 *
 *	A \ref Volume is a `d` dimensional finite element.
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

	const bool boundary, ///< Flag for whether the volume is on a domain boundary.
	           curved;   ///< Flag for whether the volume is curved.

	const struct const_Matrix_d*const xyz_ve; ///< The xyz coordinates of the volume vertices.

/// \todo Potentially make this a Multiarray.
	/** The geometry coefficients of the volume in the \ref Simulation::basis_geom. For each of the supported
	 * \ref Simulation::domain_type options, geom_coef represents:
	 *	- DOM_STRAIGHT: the projection of xyz_ve into the geometry basis of order 1.
	 *	- DOM_CURVED:
	 *		- straight volumes: [See DOM_STRAIGHT];
	 *		- boundary volumes: the coefficients of the *blended* face geometry of order k_g.
	 *	- DOM_PARAMETRIC: the coefficients of the mapped geometry of order k_g.
	 */
	const struct const_Matrix_d*const geom_coef;

	const struct Face*const faces[NFMAX][NSUBFMAX]; ///< Array of pointers to the neighbouring \ref Face containers.

	const struct const_Element*const element; ///< Pointer to the associated \ref const_Element.
};

// Constructor/Destructor functions ********************************************************************************* //

/** \brief Constructs the base \ref Volume \ref Intrusive_List.
 *	\return Standard. */
struct Intrusive_List* constructor_Volumes
	(struct Simulation*const sim, ///< The \ref Simulation.
	 const struct Mesh*const mesh ///< The \ref Mesh.
	);

/// \brief Destructs the base \ref Volume \ref Intrusive_List.
void destructor_Volumes
	(struct Intrusive_List* volumes ///< Standard.
	);

/// \brief Destructor for an individual \ref Volume.
void destructor_Volume
	(struct Volume* volume ///< Standard.
	);

// Helper functions ************************************************************************************************* //

/** \brief Check if a sufficient number of vertices (2) satisfy the condition on a domain boundary.
 *
 *	When two vertices satisfy the condition, this implies that a volume edge satisfies the condition and the
 *	qualification is passed to the associated volume.
 *
 *	This function can be used to check:
 *		- faces of volumes: `ve_inds` = volume vertex indices, `f_ve` = \ref Volume::element;
 *		- edges of faces:   `ve_inds` = face vertex indices,   `f_ve` = \ref Face::element.
 *
 *	Currently used to check:
 *		- curved (curved_only = `true`);
 *		- boundary (curved only = `false`).
 *
 *	\return See brief. */
bool check_ve_condition
	(const struct const_Multiarray_Vector_i*const f_ve,  ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,          ///< The vertex indices for the volume.
	 const struct const_Vector_i*const ve_condition,     ///< The vertex condition to check for.
	 const struct const_Multiarray_Vector_i*const ve_bc, ///< Defined in \ref Mesh_Vertices.
	 const bool curved_only                              /**< Flag indicating whether only curved boundary conditions
	                                                      *   should be considered. */
	);

#endif // DPG__volume_h__INCLUDED
