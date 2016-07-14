// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>

#include "database.h"
#include "functions.h"

/*
 *	Purpose:
 *		Free remaining allocated memory.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	Memory freed:
 *
 *
 *	References:
 *
 */

void memory_free(void)
{
	// DB Parameters

		// Initialization
		free(DB.TestCase), free(DB.MeshType), free(DB.Form), free(DB.NodeType), free(DB.BasisType), free(DB.MeshFile);

		// Preprocessing

			// SetupParameters
			free(DB.Parametrization);
			array_free3_c(DB.NP,DB.NEC,DB.NodeTypeS);
			array_free3_c(DB.NP,DB.NEC,DB.NodeTypeF);
			array_free3_c(DB.NP,DB.NEC,DB.NodeTypeFrs);
			array_free3_c(DB.NP,DB.NEC,DB.NodeTypeFrc);
			array_free3_c(DB.NP,DB.NEC,DB.NodeTypeIfs);
			array_free3_c(DB.NP,DB.NEC,DB.NodeTypeIfc);
			array_free3_c(DB.NP,DB.NEC,DB.NodeTypeIvs);
			array_free3_c(DB.NP,DB.NEC,DB.NodeTypeIvc);

			free(DB.PGc);
			free(DB.PF);
			array_free3_ui(DB.NP,2,DB.SF_BE);
			array_free2_ui(DB.NP,DB.PCs);
			array_free2_ui(DB.NP,DB.PCc);
			array_free2_ui(DB.NP,DB.PJs);
			array_free2_ui(DB.NP,DB.PJc);
			array_free2_ui(DB.NP,DB.PFrs);
			array_free2_ui(DB.NP,DB.PFrc);
			array_free2_ui(DB.NP,DB.PIfs);
			array_free2_ui(DB.NP,DB.PIfc);
			array_free2_ui(DB.NP,DB.PIvs);
			array_free2_ui(DB.NP,DB.PIvc);
			free(DB.VFPartUnity);

			// SetupMesh
			free(DB.PVe), free(DB.NE), free(DB.EType), free(DB.ETags), free(DB.EToVe), free(DB.EToPrt);
			free(DB.VToV), free(DB.VToF), free(DB.VToGF), free(DB.VToBC), free(DB.GFToVe), free(DB.VC), free(DB.GFC);
			free(DB.VeXYZ);

	// ELEMENTs
	struct S_ELEMENT *ELEMENT, *ELEMENTnext;

	ELEMENT = DB.ELEMENT;
	while (ELEMENT != NULL) {
		ELEMENTnext = ELEMENT->next;
		memory_destructor_E(ELEMENT);
		ELEMENT = ELEMENTnext;
	}

	// VOLUMEs
	struct S_VOLUME *VOLUME, *VOLUMEnext;
	VOLUME = DB.VOLUME;
	while (VOLUME != NULL) {
		VOLUMEnext = VOLUME->next;
		memory_destructor_V(VOLUME);
		VOLUME = VOLUMEnext;
	}

	// FACETs
	struct S_FACET *FACET, *FACETnext;
	for (FACET = DB.FACET; FACET != NULL; ) {
		FACETnext = FACET->next;
		memory_destructor_F(FACET);
		FACET = FACETnext;
	}
}

void memory_free_children(void)
{
	/*
	 *	Purpose:
	 *		Free memory of non-leaf VOLUME and FACET children.
	 */

// freeing past the coarsed VOLUME/FACET? Investigate (ToBeDeleted)

	// VOLUMEs
	struct S_VOLUME *VOLUME, *VOLUMEc, *VOLUMEcnext;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOLUMEc = VOLUME->child0;
		while (VOLUMEc) {
			VOLUMEcnext = VOLUMEc->next;
			memory_destructor_V(VOLUMEc);
			VOLUMEc = VOLUMEcnext;
		}
		VOLUME->child0 = NULL;
	}

	// FACETs
	struct S_FACET *FACET, *FACETc, *FACETcnext;
	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		FACETc = FACET->child0;
ptrdiff_t count = 0;
		while (FACETc) {
			FACETcnext = FACETc->next;
			memory_destructor_F(FACETc);
			FACETc = FACETcnext;
printf("FACETfree: %d\n",count++);
		}
		FACET->child0 = NULL;
	}

}
