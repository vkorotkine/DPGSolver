// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Face_h__INCLUDED
#define DPG__Face_h__INCLUDED
/**	\file
 *	Provides the interface for the base \ref Face container and associated functions.
 */

#include <stdbool.h>

#include "Intrusive.h"

#include "Simulation.h"
#include "Mesh.h"

/// \brief Container for data relating to the base Volumes.
struct Face {
	struct Intrusive_Link lnk; ///< \ref Intrusive_Link.

	const bool boundary, ///< Flag for whether the volume is on a domain boundary.
	           curved;   ///< Flag for whether the volume is curved.

	const struct const_Element*const element; ///< Pointer to the associated \ref const_Element.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructs the base \ref Volume \ref Intrusive_List.
struct Intrusive_List* constructor_Face_List
	(const struct Simulation*const sim, ///< The \ref Simulation.
	 const struct Mesh*const mesh       ///< The \ref Mesh.
	);

/// \brief Destructs the base \ref Face \ref Intrusive_List.
void destructor_Faces
	(struct Intrusive_List* Faces ///< Standard.
	);

// Helper functions ************************************************************************************************* //

#endif // DPG__Face_h__INCLUDED
