// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Volume_h__INCLUDED
#define DPG__Volume_h__INCLUDED
/**	\file
 *	Provides the interface for the base \ref Volume container and associated functions.
 *
 *	A \ref Volume is a `d` dimensional finite element.
 */

#include <stdbool.h>

#include "Intrusive.h"
#include "Matrix.h"

#include "Simulation.h"
#include "Mesh.h"

#include "constants_elements.h"

/// \brief Container for data relating to the base Volumes.
struct Volume {
	struct Intrusive_Link lnk; ///< The \ref Intrusive_Link.

	const bool boundary, ///< Flag for whether the volume is on a domain boundary.
	           curved;   ///< Flag for whether the volume is curved.

	const struct const_Matrix_d*const xyz_ve; ///< The xyz coordinates of the volume vertices.

	const struct Face*const faces[NFMAX][NSUBFMAX]; ///< Array of pointers to the neighbouring \ref Face containers.

	const struct const_Element*const element; ///< Pointer to the associated \ref const_Element.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructs the base \ref Volume \ref Intrusive_List.
struct Intrusive_List* constructor_Volume_List
	(struct Simulation*const sim, ///< The \ref Simulation.
	 const struct Mesh*const mesh ///< The \ref Mesh.
	);

/// \brief Destructs the base \ref Volume \ref Intrusive_List.
void destructor_Volumes
	(struct Intrusive_List* Volumes ///< Standard.
	);

// Helper functions ************************************************************************************************* //

#endif // DPG__Volume_h__INCLUDED
